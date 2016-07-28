#include "ChemotaxisWorld.h"

//Initialize the parameter list.
std::shared_ptr<ParameterLink<bool>> ChemotaxisWorld::use_lin_gradient_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-use_lin_gradient", true, "Create a linear attractant gradient. Otherwise, it will be exponential.");
std::shared_ptr<ParameterLink<bool>> ChemotaxisWorld::clear_outputs_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-clear_outputs", true, "Reset output nodes between brain updates.");
std::shared_ptr<ParameterLink<bool>> ChemotaxisWorld::environment_variability_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-environment_variability", false, "Cause the environment to vary between replicate runs.");
std::shared_ptr<ParameterLink<bool>> ChemotaxisWorld::use_integral_sensor_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-use_integral_sensor", false, "Use the (experimental) integral sensor.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::rot_diff_coeff_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-rot_diff_coeff", 0.016, "Rotational diffusion constant.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::speed_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-speed", 1.0, "Magnitude of cell movement per tick.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::slope_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-slope", 10.0, "Rate of gradient increase. m if linear, k if exponential.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::base_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-base", 255.0, "Base concentration at x = 0.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::variability_slope_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-variability_slope", 1.0, "Slope will change by as much as +/- this number.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::variability_base_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-variability_base", 1.0, "Base will change by as much as +/- this number.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::variability_rot_diff_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-variability_rot_diff", 0.01, "Rotational diffusion constant will change by as much as +/- this number.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::variability_speed_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-variability_speed", 0.2, "Speed will change by as much as +/- this number.");
std::shared_ptr<ParameterLink<int>> ChemotaxisWorld::eval_ticks_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-eval_ticks", 5000, "Number of ticks to evaluate.");
std::shared_ptr<ParameterLink<int>> ChemotaxisWorld::brain_updates_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-brain_updates", 1, "Number of times the brain is set to update before the output is read.");

ChemotaxisWorld::ChemotaxisWorld(std::shared_ptr<ParametersTable> _PT) : //Initializer
  AbstractWorld(_PT) {
    //Grab the values from the parameter list.
    use_lin_gradient = (PT == nullptr) ? use_lin_gradient_pl->lookup() : PT->lookupBool("WORLD_CHEMOTAXIS-use_lin_gradient");
    clear_outputs = (PT == nullptr) ? clear_outputs_pl->lookup() : PT->lookupBool("WORLD_CHEMOTAXIS-clear_outputs");
    environment_variability = (PT == nullptr) ? environment_variability_pl->lookup() : PT->lookupBool("WORLD_CHEMOTAXIS-environment_variability");
    use_integral_sensor = (PT == nullptr) ? use_integral_sensor_pl->lookup() : PT->lookupBool("WORLD_CHEMOTAXIS-use_integral_sensor");
    rot_diff_coeff = (PT == nullptr) ? rot_diff_coeff_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-rot_diff_coeff");
    speed = (PT == nullptr) ? speed_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-speed");
    slope = (PT == nullptr) ? slope_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-slope");
    base = (PT == nullptr) ? base_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-base");
    variability_slope = (PT == nullptr) ? variability_slope_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-slope");
    variability_base = (PT == nullptr) ? variability_base_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-base");
    variability_rot_diff = (PT == nullptr) ? variability_rot_diff_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-rot_diff_coeff");
    variability_speed = (PT == nullptr) ? variability_speed_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-speed");
    eval_ticks = (PT == nullptr) ? eval_ticks_pl->lookup() : PT->lookupInt("WORLD_CHEMOTAXIS-eval_ticks");
    brain_updates = (PT == nullptr) ? brain_updates_pl->lookup() : PT->lookupInt("WORLD_CHEMOTAXIS-brain_updates");

    //Set up data columns.
    aveFileColumns.clear();
    aveFileColumns.push_back("x_displacement");
    aveFileColumns.push_back("y_displacement");
}

//Should return a linear concentration unless it's negative, in which case return 0.
//Will not attempt to convert negative value, and will return 0 instead.
inline double get_conc_linear(const double &x, const double &slope, const double &base) {
    double conc = (x * slope) + base;
  return ((conc > 0) ? conc : 0);
}

//Should return an exponential concentration, unless it's negative, in which case return 0.
//Will not attempt to convert negative value, and will return 0 instead.
inline double get_conc_exp(const double &x, const double &k, const double &base) {
  double conc = exp(k * x) + base;
  return ((conc > 0) ? conc : 0);
}

//Use the cell's x position, y position, angle theta, and the cell's speed
//to calculate the new position. Void because it works directly on the vector by ref.
//pos_vec has the form [x,y,theta,speed]
inline void run(std::vector<double> &pos_vec) {
  pos_vec[0] += (pos_vec[3] * cos(pos_vec[2]));
  pos_vec[1] += (pos_vec[3] * sin(pos_vec[2]));
  return;
}

//The cell tumbles; resulting in a new theta between 0 and 2*pi
//Potentially add in ability to bias later.
inline double tumble() {
  return (Random::getDouble(2 * M_PI));
}

//Causes random changes in theta which depend on the diffusion coefficient.
inline void rot_diffuse(double& theta, const double &diff_coeff) {
  theta += (Random::getNormal(0,1) * diff_coeff);
  return;
}

//Take a vector of bits and convert it into a probability of tumbling.
//The endianness of how the brain returns outputs shouldn't matter.
//The organism will 'learn' it. Will only take 8 bits for now; make generic size if needed.
//Takes vector of ints, which should be 1 or 0, multiply each position by 2**n and add
//the result to an accumulator to convert to int. Then take accumulator+1 (so bias can't be 0)
//and divide by 256 to get the resulting tumble bias.
double bit_to_prob(const bool (&bit_arr)[16]) {
  double accumulator = 0;
  for(int idx = 0; idx != 16; ++idx){
    accumulator += bit_arr[idx];
  }
  return ((accumulator+1)/17);
}

//Determine whether the nth bit of an uint32_t is 1. Will fail if idx > 31.
inline int is_set(const uint32_t &num, const int &idx){
  return ((num & (1 << idx)) >> idx);
}

void ChemotaxisWorld::runWorldSolo(std::shared_ptr<Organism> org, bool analyse, bool visualize, bool debug) {
  //Should be able to thread and run concurrently at the loss of reproducibility.
  //Starting orientation is random between 0 and 2pi.
  //Also starts the organism on the origin for now.
  std::vector<double> pos_vec{0.0, 0.0, Random::getDouble(0,2 * M_PI), speed};
  std::vector<std::vector<double>> pos_hist;
  std::vector<double> tumble_hist;
  std::vector<double> concentration_hist;
  std::vector<double> delta_hist;
  std::vector<int> memory_hist;
  bool output_array[16];


  //Pre-reserve capacity to prevent reallocations.
  pos_hist.reserve(eval_ticks+1);
  tumble_hist.reserve(eval_ticks+1);
  concentration_hist.reserve(eval_ticks+1);
  memory_hist.reserve(eval_ticks+1);
  delta_hist.reserve(eval_ticks+1);

  //Initialize conc_hist to keep things simpler when calculating delta.
  concentration_hist.push_back((use_lin_gradient) ? get_conc_linear(pos_vec[0], slope, base) : get_conc_exp(pos_vec[0], slope, base));

  //Initialize the cell.
  bool is_tumbling;
  double tumble_bias;
  double concentration;
  double delta;

  //Sanity check reset the brain. Shouldn't need to do this
  org->brain->resetBrain();

  /*
  //Add in environmental_variability, if enabled.
  if (environment_variability) {
    //Re-grab the initial values to modify.
    rot_diff_coeff = (PT == nullptr) ? rot_diff_coeff_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-rot_diff_coeff");
    speed = (PT == nullptr) ? speed_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-speed");
    slope = (PT == nullptr) ? slope_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-slope");
    base = (PT == nullptr) ? base_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-base");

    rot_diff_coeff += Random::getDouble((-variability_rot_diff), variability_rot_diff);
    speed += Random::getDouble((-variability_speed), variability_speed);
    slope += Random::getDouble((-variability_slope), variability_slope);
    base += Random::getDouble((-variability_slope), variability_slope);

    //Catch negatives and just set them to zero. You could maybe skip this, but let's do it for now.
    if (rot_diff_coeff < 0) {rot_diff_coeff = 0;}
    if (speed < 0) {speed = 0;}
    if (slope < 0) {slope = 0;}
    if (base < 0) {base = 0;}
  } // End environmental variability
  */
  //Evaluate the organism for eval_ticks ticks.
  for (int t = 0; t != eval_ticks; ++t){
    //Get concentration at location.
    if(use_lin_gradient){
      concentration = get_conc_linear(pos_vec[0], slope, base);
      concentration_hist.push_back(concentration);
    }
    else{
      concentration = get_conc_exp(pos_vec[0], slope, base);
      concentration_hist.push_back(concentration);
    }

    //Integral sensor
    /* This sensor will integrate over a timescale determined by the brain.
      It will average the values over the timescale, then subtract the average
      from the current value. */
    //Get the memory length from the first 8 bits of the brain. The length is just the number of 1s.
    /*
    int mem_length = 0;
    for (n = 0; n != 8; ++n){
      mem_length += org->brain->readOutput(n);
    }
    mem_hist.push_back(mem_length);

    //Calculate the relative change in concentration, add it to the delta history.
    //conc_hist is initialized with a value, so no need to worry about it being empty.
    //If the last element is zero, add a tiny value so things don't break. Should be pretty rare in practice.
    if(conc_history.back() == 0){conc_hist[conc_hist.size()] += 0.00001;}
    double delta = (concentration - conc_history[conc_history.back()]) / std::abs(conc_history[conc_history.back()]);
    delta_hist.push_back(delta);

    //Go through the last mem_length elements and average them.
    double accumulator = 0;
    int counter = 0;
    if (mem_length > 0) {
      auto idx_ptr = cend(delta_hist);
      while(counter != mem_length){
        --idx_ptr;
        ++counter;
        accumulator += *idx_ptr;
      }
    }
    double bits_to_add = (mem_length == 0) ? delta : accumulator;

    //Convert the delta to the unary output. A delta of 0 should produce 8 1s and 8 0s
    //Mult by two to give a bit more range
    int num_ones = std::round(avg_delta) + 8;
    if (num_ones < 0){num_ones = 0;}
    if (num_ones > 16){num_ones = 16;}

    for (n = 0; n != 16; ++n){
      if (num_ones != 0){
        org->brain->setInput(n,1)
        --num_ones;
      }
      else{
        org->brain->setInput(n,0)
      }
    }
    */

    //Let's just see how they behave with a delta sensor (1 if gradient increased, 0 otherwise);
    delta = concentration - concentration_hist[concentration_hist.size()-2];
    delta_hist.push_back(delta);
    if(delta > 0){
      org->brain->setInput(0,1);
    }
    else{
      org->brain->setInput(0,0);
    }

    //Update the brain after resetting outputs (if enabled)
    if (clear_outputs){
      org->brain->resetOutputs();
    }
    for(int runs = 0; runs != brain_updates; ++runs){
    org->brain->update();
    }
    //Read the output and convert to tumble bias.
    for (int n = 0; n != 16; ++n){
      output_array[n] = org->brain->readOutput(n);
    }
    tumble_bias = bit_to_prob(output_array);
    tumble_hist.push_back(tumble_bias);

    //Do the tumbling.
    is_tumbling = Random::P(tumble_bias);
    if (is_tumbling){
      pos_vec[2] = tumble();
      pos_hist.push_back(pos_vec);
    }
    else{ //And running
      run(pos_vec);
      pos_hist.push_back(pos_vec);
      rot_diffuse(pos_vec[2], rot_diff_coeff);
    }
  }//end eval loop

  //Score based on x movement and record some stats.
  org->score = pos_vec[0];
  org->dataMap.Append("allx_displacement", (double) pos_vec[0]);
  org->dataMap.Append("ally_displacement", (double) pos_vec[1]);

  //If in visualize mode, dump the points to file.
  if (visualize){
    //Position data
    std::ofstream out_file("chemotaxis_position_visualization_data.csv");
    for (auto sample : pos_hist){
      out_file << sample[0] << "," << sample[1] << "," << sample[2] << "," << '\n';
    }
    out_file << std::endl;
    //Tumble data
    std::ofstream out_file_tumble("chemotaxis_tumble_visualization_data.csv");
    for (auto sample : tumble_hist){
      out_file_tumble << sample << '\n';
    }
    out_file_tumble << std::endl;
    //Concentration data
    std::ofstream out_file_conc("chemotaxis_conc_visualization_data.csv");
    for (auto sample : concentration_hist){
      out_file_conc << sample << '\n';
    }
      out_file_conc << std::endl;
    //Delta data
    std::ofstream out_file_delta("chemotaxis_delta_visualization_data.csv");
    for (auto sample : delta_hist){
      out_file_delta << sample << '\n';
    }
    out_file_delta << std::endl;


  }//End visualize
}//End of RunWorldSolo fn


//The integral sensor will have 16 bits describing whether concentration is increasing or decreasing.
int ChemotaxisWorld::requiredInputs() {
    //return 16;
    return 1;
}

//The first 8 bits of output are the memory length (0-8) while the next 16 are the tumble bias.
int ChemotaxisWorld::requiredOutputs() {
  //return 24;
  return 16;
}

//Solo eval. No reason to run more than once.
int ChemotaxisWorld::maxOrgsAllowed() {
  return 1;
}

//Obviously need at least one organism.
int ChemotaxisWorld::minOrgsAllowed() {
  return 1;
}
