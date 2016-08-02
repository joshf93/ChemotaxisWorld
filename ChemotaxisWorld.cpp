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
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::variability_slope_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-variability_slope", 1.0, "Slope will increase by as much as this number.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::variability_base_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-variability_base", 1.0, "Base will increase by as much as this number.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::variability_rot_diff_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-variability_rot_diff", 0.01, "Rotational diffusion constant will increase by as much as this number.");
//std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::variability_speed_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-variability_speed", 0.2, "Speed will increase by as much as this number.");
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
    //variability_speed = (PT == nullptr) ? variability_speed_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-speed");
    eval_ticks = (PT == nullptr) ? eval_ticks_pl->lookup() : PT->lookupInt("WORLD_CHEMOTAXIS-eval_ticks");
    brain_updates = (PT == nullptr) ? brain_updates_pl->lookup() : PT->lookupInt("WORLD_CHEMOTAXIS-brain_updates");

    //Set up data columns.
    aveFileColumns.clear();
    aveFileColumns.push_back("x_displacement");
    aveFileColumns.push_back("y_displacement");

    //Print some stuff so we know what's going on.
    std::cout << "#########################" << "\n";
    std::cout << "Running ChemotaxisWorld." << "\n";
    std::cout << "Cycles to simulate is: " << eval_ticks << "\n";
    std::cout << "Brain updates per cycle is: " << brain_updates << "\n";
    std::cout << "Organism speed is: " << speed << "\n";
    std::cout << "Rotational diffusion coefficient is: " << rot_diff_coeff << "\n";
    std::cout << "Using a variable environment: " << environment_variability << "\n";
    std::cout << "Using a linear gradient: " << use_lin_gradient << "\n";
    std::cout << "Slope is: " << slope << "\n";
    std::cout << "#########################" << std::endl;

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

//Unused. Left here just in case.
//Determine whether the nth bit of an uint32_t is 1. Will fail if idx > 31.
//inline int is_set(const uint32_t &num, const int &idx){
//  return ((num & (1 << idx)) >> idx);
//}


void ChemotaxisWorld::runWorldSolo(std::shared_ptr<Organism> org, bool analyse, bool visualize, bool debug) {
  //Should be able to thread and run concurrently at the loss of reproducibility.
  //Starting orientation is random between 0 and 2pi.
  //Also starts the organism on the origin for now.
  std::vector<double> pos_vec{0.0, 0.0, Random::getDouble(0,2 * M_PI), speed};
  std::vector<std::vector<double>> pos_hist;
  std::vector<double> tumble_hist;
  std::vector<double> concentration_hist;
  std::vector<double> delta_hist;
  std::vector<int> multiplier_hist;
  std::vector<int> ones_hist;

  //Pre-reserve capacity to prevent reallocations.
  pos_hist.reserve(eval_ticks+2);
  tumble_hist.reserve(eval_ticks+2);
  concentration_hist.reserve(eval_ticks+2);
  multiplier_hist.reserve(eval_ticks+2);
  delta_hist.reserve(eval_ticks+2);

  //Initialize conc_hist to keep things simpler when calculating delta.
  double initial_conc = (use_lin_gradient) ? get_conc_linear(pos_vec[0], slope, base) : get_conc_exp(pos_vec[0], slope, base);
  concentration_hist.push_back(initial_conc);

  //Initialize the cell.
  bool is_tumbling;
  double tumble_bias;
  double concentration;
  double delta;
  double accumulator;
  int multiplier;
  int num_ones;

  //Sanity check reset the brain. Shouldn't need to do this
  org->brain->resetBrain();

  //Add in environmental_variability, if enabled.
  if (environment_variability) {
    //Re-grab the initial values to modify.
    rot_diff_coeff = (PT == nullptr) ? rot_diff_coeff_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-rot_diff_coeff");
    //speed = (PT == nullptr) ? speed_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-speed");
    slope = (PT == nullptr) ? slope_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-slope");
    base = (PT == nullptr) ? base_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-base");

    //Add in the specified amount of noise.
    rot_diff_coeff += Random::getDouble(variability_rot_diff);
    //speed += Random::getDouble(variability_speed); This doesn't work; for some reason it will change even if var_speed is 0.
    slope += Random::getDouble(variability_slope);
    base += Random::getDouble(variability_base);

    //Catch negatives and just set them to zero. You could maybe skip this, but let's do it for now
    //in case someone sets negative values.
    if (rot_diff_coeff < 0) {rot_diff_coeff = 0;}
    if (speed < 0) {speed = 0;}
    if (slope < 0) {slope = 0;}
    if (base < 0) {base = 0;}
  } // End environmental variability

  //Evaluate the organism for eval_ticks ticks.
  for (int t = 0; t != eval_ticks; ++t) {
    //Get concentration at location using the appropriate method.
    concentration = (use_lin_gradient) ? get_conc_linear(pos_vec[0], slope, base) : get_conc_exp(pos_vec[0], slope, base);

    //Read in the multiplier
    multiplier = 0;
    for (int n = 16; n != 24; ++n) {
      multiplier += org->brain->readOutput(n);
    }
    //Because of how brains work, this number could be almost any positive int. Cap it at 31 to
    //prevent undefined behavior.
    //This is strange because I was pretty sure that values were OR'd, so they could only be binary.
    //Turns out that isn't the case, they seem to add. I have gotten values of 22-24 before
    //and presumably ones shifting over 31 were just selected against; they would 'wrap around'
    //and shift by multiplier%32 on x86.
    if (multiplier > 31) {
      multiplier = 31;
    }
    //In case logic gates can output negative numbers.
    if (multiplier < 0) {
      multiplier = 0;
    }
    multiplier_hist.push_back(multiplier);

    //Delta is the relative concentration change between now and the last sample.
    delta = (concentration_hist.back() == 0) ? (concentration/0.000001) : (concentration - concentration_hist.back()) / std::abs(concentration_hist.back());
    concentration_hist.push_back(concentration); //Put here to avoid annoying side effects.
    delta_hist.push_back(delta);
    delta *= (1 << multiplier); //Lshift by multiplier to allow org to control sensitivity to change.

    //Convert the delta to the unary output. A delta of 0 should produce 8 1s and 8 0s
    num_ones = std::round(delta) + 8;
    ones_hist.push_back(num_ones);
    //Cap at 0 and 16 ones.
    if (num_ones < 0) {
      num_ones = 0;
    }
    if (num_ones > 16) {
      num_ones = 16;
    }

    //Write the input bits
    for (int n = 0; n != 16; ++n) {
      if (num_ones != 0){
        org->brain->setInput(n,1);
        --num_ones;
      }
      else {
        org->brain->setInput(n,0);
      }
    }

    //Update the brain after resetting outputs (if enabled)
    if (clear_outputs) {
      org->brain->resetOutputs();
    }
    for(int runs = 0; runs != brain_updates; ++runs) {
      org->brain->update();
    }

    //Read the output and convert to tumble bias.
    accumulator = 0;
    for (int n = 0; n != 16; ++n) {
      accumulator += org->brain->readOutput(n);
    }
    if (accumulator > 16) {
      accumulator = 16;
    }
    //Ensure the accumulator isn't negative in case gates that can output -1 are used.
    tumble_bias = (accumulator <= 0) ? (1.0/17.0) : (accumulator + 1.0)/17.0;
    tumble_hist.push_back(tumble_bias);

    //Do the tumbling.
    is_tumbling = Random::P(tumble_bias);
    if (is_tumbling) {
      pos_vec[2] = tumble();
      pos_hist.push_back(pos_vec);
    }
    else { //And running
      run(pos_vec);
      rot_diffuse(pos_vec[2], rot_diff_coeff);
      pos_hist.push_back(pos_vec);
    }
  }//end eval loop

  //Finished simulating the organism, so score based on x movement and record some stats.
  org->score = pos_vec[0];
  org->dataMap.Append("allx_displacement", (double) pos_vec[0]);
  org->dataMap.Append("ally_displacement", (double) pos_vec[1]);

  //debugging
  if (pos_vec[0] > eval_ticks){
    visualize = true;
  }

  //If in visualize mode, dump the history points to file.
  if (visualize) {
    std::cout << "Rot diff coeff was: " << rot_diff_coeff << '\n';
    std::cout << "Slope was: " << slope << '\n';
    std::cout << "Speed was: " << speed << '\n';
    std::cout << "Base was: " << base << '\n';
    std::cout << std::endl;
    //Position data
    std::ofstream out_file("chemotaxis_position_visualization_data.csv");
    for (auto sample : pos_hist) {
      out_file << sample[0] << "," << sample[1] << "," << sample[2] << "," << '\n';
    }
    out_file << std::endl;
    //Tumble data
    std::ofstream out_file_tumble("chemotaxis_tumble_visualization_data.csv");
    for (auto sample : tumble_hist) {
      out_file_tumble << sample << '\n';
    }
    out_file_tumble << std::endl;
    //Concentration data
    std::ofstream out_file_conc("chemotaxis_conc_visualization_data.csv");
    for (auto sample : concentration_hist) {
      out_file_conc << sample << '\n';
    }
      out_file_conc << std::endl;
    //Delta data
    std::ofstream out_file_delta("chemotaxis_delta_visualization_data.csv");
    for (auto sample : delta_hist) {
      out_file_delta << sample << '\n';
    }
    out_file_delta << std::endl;
    //Multiplier data
    std::ofstream out_file_multiplier("chemotaxis_multiplier_visualization_data.csv");
    for (auto sample : multiplier_hist) {
      out_file_multiplier << sample << '\n';
    }
    out_file_multiplier << std::endl;
    //Ones count
    std::ofstream out_file_ones("chemotaxis_ones_visualization_data.csv");
    for(auto sample : ones_hist) {
      out_file_ones << sample << '\n';
    }
    out_file_ones << std::endl;
  }//End visualize

  if (visualize == true && pos_vec[0] > eval_ticks){
    cout << "ERROR: MOVED FURTHER THAN POSSIBLE" << std::endl;
    exit(1);
  }
}//End of RunWorldSolo fn

//The sensor will have 16 bits describing whether concentration is increasing or decreasing.
int ChemotaxisWorld::requiredInputs() {
  return 16;
}

//The first 16 bits of output are the tumble bias (1/17-17/17) while the next 8 are the sensitivity multiplier.
int ChemotaxisWorld::requiredOutputs() {
  return 24;
}

//Solo eval. No reason to run more than once.
int ChemotaxisWorld::maxOrgsAllowed() {
  return 1;
}

//Obviously need at least one organism.
int ChemotaxisWorld::minOrgsAllowed() {
  return 1;
}
