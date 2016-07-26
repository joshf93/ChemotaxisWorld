#include "ChemotaxisWorld.h"

//Initialize the parameter list.
std::shared_ptr<ParameterLink<bool>> ChemotaxisWorld::use_lin_gradient_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-use_lin_gradient", true, "Create a linear attractant gradient. Otherwise, it will be exponential.");
std::shared_ptr<ParameterLink<bool>> ChemotaxisWorld::clear_outputs_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-clear_outputs", true, "Reset output nodes between brain updates.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::rot_diff_coeff_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-rot_diff_coeff", 0.016, "Rotational diffusion constant.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::speed_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-speed", 1.0, "Magnitude of cell movement per tick.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::slope_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-slope", 10.0, "Rate of gradient increase. m if linear, k if exponential.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::base_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-base", 255.0, "Base concentration at x = 0.");
std::shared_ptr<ParameterLink<int>> ChemotaxisWorld::eval_ticks_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-eval_ticks", 5000, "Number of ticks to evaluate.");

ChemotaxisWorld::ChemotaxisWorld(std::shared_ptr<ParametersTable> _PT) :
  AbstractWorld(_PT) {
    //Grab the values from the parameter list.
    use_lin_gradient = (PT == nullptr) ? use_lin_gradient_pl->lookup() : PT->lookupBool("WORLD_CHEMOTAXIS-use_lin_gradient");
    clear_outputs = (PT == nullptr) ? clear_outputs_pl->lookup() : PT->lookupBool("WORLD_CHEMOTAXIS-clear_outputs");
    rot_diff_coeff = (PT == nullptr) ? rot_diff_coeff_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-rot_diff_coeff");
    speed = (PT == nullptr) ? speed_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-speed");
    slope = (PT == nullptr) ? slope_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-slope");
    base = (PT == nullptr) ? base_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-base");
    eval_ticks = (PT==nullptr) ? eval_ticks_pl->lookup() : PT->lookupInt("WORLD_CHEMOTAXIS-eval_ticks");

    //Set up data columns.
    aveFileColumns.clear();
    aveFileColumns.push_back("x_displacement");
    aveFileColumns.push_back("y_displacement");
}

//Should return a linear concentration unless it's negative, in which case return 0.
//Will not attempt to convert negative value, and will return 0 instead.
inline uint32_t get_conc_linear(const double &x, const double &slope, const double &base){
    double conc = (x * slope) + base;
  return ((conc > 0) ? conc : 0);
}

//Should return an exponential concentration, unless it's negative, in which case return 0.
//Will not attempt to convert negative value, and will return 0 instead.
inline uint32_t get_conc_exp(const double &x, const double &k, const double &base){
  double conc = exp(k * x) + base;
  return ((conc > 0) ? conc : 0);
}

//Use the cell's x position, y position, angle theta, and the cell's
//speed to calculate the new position. Void because it works directly on the vector by ref.
//pos_vec has the form [x,y,theta,speed]
inline void run(std::vector<double> &pos_vec){
  pos_vec[0] += (pos_vec[3] * cos(pos_vec[2]));
  pos_vec[1] += (pos_vec[3] * sin(pos_vec[2]));
  return;
}

//The cell tumbles; resulting in a new theta between 0 and 2*pi
//Potentially add in ability to bias later.
inline double tumble(){
  return (Random::getDouble(2 * M_PI));
}

//Causes random changes in theta which depend on the diffusion coefficient.
inline void rot_diffuse(double& theta, const double &diff_coeff){
  theta += (Random::getNormal(0,1) * diff_coeff);
  return;
}

//Take a vector of bits and convert it into a probability of tumbling.
//The endianness of how the brain returns outputs shouldn't matter.
//The organism will 'learn' it. Will only take 8 bits for now; make generic size if needed.
//Takes vector of ints, which should be 1 or 0, multiply each position by 2**n and add
//the result to an accumulator to convert to int. Then take accumulator+1 (so bias can't be 0)
//and divide by 256 to get the resulting tumble bias.
double bit_to_prob(const bool (&bit_arr)[8]){
  int accumulator = 0;
  for(int idx = 0; idx != 8; ++idx){
    accumulator += (bit_arr[idx]*(1 << idx));
  }
  return ((accumulator+1)/256);
}

//Determine whether the nth bit of an uint32_t is 1. Will fail if idx > 31.
inline int is_set(const uint32_t &num, const int &idx){
  return ((num & (1 << idx)) >> idx);
}

void ChemotaxisWorld::runWorldSolo(std::shared_ptr<Organism> org, bool analyse, bool visualize, bool debug){
  //Starting orientation should be random.
  //Also starts the organism on the origin for now.
  std::vector<double> pos_vec{0.0, 0.0, Random::getDouble(0,2 * M_PI), speed};  //Make this an option later?
  std::vector<std::vector<double>> pos_hist;
  std::vector<double> tumble_hist;
  std::vector<double> concentration_hist;
  bool output_array[8];

  //Pre-reserve to prevent reallocations.
  pos_hist.reserve(eval_ticks+1);
  tumble_hist.reserve(eval_ticks+1);
  concentration_hist.reserve(eval_ticks+1);

  //Initialize the cell.
  bool is_tumbling;
  double tumble_bias;
  uint32_t concentration;

  //Evaluate the organism for eval_ticks ticks.
  for(int t = 0; t != eval_ticks; ++t){
    //Get concentration at location.
    if(use_lin_gradient){
      concentration = get_conc_linear(pos_vec[0], slope, base);
      concentration_hist.push_back(concentration);
    }
    else{
      concentration = get_conc_exp(pos_vec[0], slope, base);
      concentration_hist.push_back(concentration);
    }

    //Put info into the brain.
    int brain_idx = 0;
    //Go through and send each bit in concentration to the brain input.
    while(brain_idx != 32){
      org->brain->setInput(brain_idx,is_set(concentration,brain_idx));
      ++brain_idx;
    }

    //Update the brain after resetting outputs (if enabled)
    if (clear_outputs){
      org->brain->resetOutputs();
    }
    org->brain->update();

    //Read the output and convert to tumble bias.
    //Take an array of 'bits', convert to int, divide (resulting ulong +1) by 256.
    //+1 to prevent tumble bias from ever being 0 ((0+1)/256 != 0).
    for (int n = 0; n != 8; ++n){
      output_array[n] = org->brain->readOutput(n);
    }
    tumble_bias = bit_to_prob(output_array);
    tumble_hist.push_back(tumble_bias);

    //Do the running and tumbling.
    is_tumbling = Random::P(tumble_bias);
    if (is_tumbling){
      pos_vec[2] = tumble();
      pos_hist.push_back(pos_vec);
    }
    else{ //Running
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
  }//End visualize
}//End of RunWorldSolo fn

//This will vary depending on what, exactly, the sensor is reading, and how it is being read.
//For noww32 bit gradient; "concentration" values will range from 0 to (2^32)-1.
//Maybe set inputs and outputs to be user modified at a later point. For now, leave them fixed.
int ChemotaxisWorld::requiredInputs() {
  return 32;
}

//Only affects tumble bias. 8 bits; tumble bias will range from 1/256 to 256/256;
int ChemotaxisWorld::requiredOutputs() {
  return 8;
}

//Solo eval. No reason to run more than once.
int ChemotaxisWorld::maxOrgsAllowed() {
  return 1;
}

//Obviously need at least one organism.
int ChemotaxisWorld::minOrgsAllowed() {
  return 1;
}
