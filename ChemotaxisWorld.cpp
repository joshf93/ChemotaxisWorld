#include "ChemotaxisWorld.h"

//Initialize the parameter list.
std::shared_ptr<ParameterLink<bool>> ChemotaxisWorld::use_lin_gradient_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-use_lin_gradient", true, "Create a linear attractant gradient. Otherwise, it will be exponential.");
std::shared_ptr<ParameterLink<bool>> ChemotaxisWorld::clear_outputs_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-clear_outputs", true, "Reset output nodes between brain updates.");
std::shared_ptr<ParameterLink<bool>> ChemotaxisWorld::environment_variability_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-environment_variability", false, "Cause the environment to vary between replicate runs.");
std::shared_ptr<ParameterLink<bool>> ChemotaxisWorld::use_bit_sensor_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-use_bit_sensor", false, "Use the single bit sensor.");
std::shared_ptr<ParameterLink<bool>> ChemotaxisWorld::point_source_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-point_sources", false, "Use 'spot' sources of attractant. Still follows your choice of linear/exp gradient.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::rot_diff_coeff_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-rot_diff_coeff", 0.016, "Rotational diffusion constant.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::speed_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-speed", 1.0, "Magnitude of cell movement per tick.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::slope_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-slope", 10.0, "Rate of gradient increase. m if linear, k if exponential.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::base_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-base", 255.0, "Base concentration at x = 0.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::variability_slope_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-variability_slope", 1.0, "Slope will increase by as much as this number.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::variability_base_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-variability_base", 1.0, "Base will increase by as much as this number.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::variability_rot_diff_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-variability_rot_diff", 0.01, "Rotational diffusion constant will increase by as much as this number.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::spot_x_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-spot_x", 300.0, "The x coordinate of the attractant spot.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::spot_y_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-spot_y", 300.0, "The y coordinate of the attractant spot.");
std::shared_ptr<ParameterLink<int>> ChemotaxisWorld::eval_ticks_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-eval_ticks", 5000, "Number of ticks to evaluate.");
std::shared_ptr<ParameterLink<int>> ChemotaxisWorld::brain_updates_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-brain_updates", 1, "Number of times the brain is set to update before the output is read.");

ChemotaxisWorld::ChemotaxisWorld(std::shared_ptr<ParametersTable> _PT) : //Initializer
  AbstractWorld(_PT) {
    //Grab the values from the parameter list.
    use_lin_gradient = (PT == nullptr) ? use_lin_gradient_pl->lookup() : PT->lookupBool("WORLD_CHEMOTAXIS-use_lin_gradient");
    clear_outputs = (PT == nullptr) ? clear_outputs_pl->lookup() : PT->lookupBool("WORLD_CHEMOTAXIS-clear_outputs");
    environment_variability = (PT == nullptr) ? environment_variability_pl->lookup() : PT->lookupBool("WORLD_CHEMOTAXIS-environment_variability");
    use_bit_sensor = (PT == nullptr) ? use_bit_sensor_pl->lookup() : PT->lookupBool("WORLD_CHEMOTAXIS-use_bit_sensor");
    point_source = (PT == nullptr) ? point_source_pl->lookup() : PT->lookupBool("WORLD_CHEMOTAXIS-point_sources");
    rot_diff_coeff = (PT == nullptr) ? rot_diff_coeff_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-rot_diff_coeff");
    speed = (PT == nullptr) ? speed_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-speed");
    slope = (PT == nullptr) ? slope_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-slope");
    base = (PT == nullptr) ? base_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-base");
    variability_slope = (PT == nullptr) ? variability_slope_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-slope");
    variability_base = (PT == nullptr) ? variability_base_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-base");
    variability_rot_diff = (PT == nullptr) ? variability_rot_diff_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-rot_diff_coeff");
    spot_x = (PT == nullptr) ? spot_x_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-spot_x");
    spot_y = (PT == nullptr) ? spot_y_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-spot_y");
    eval_ticks = (PT == nullptr) ? eval_ticks_pl->lookup() : PT->lookupInt("WORLD_CHEMOTAXIS-eval_ticks");
    brain_updates = (PT == nullptr) ? brain_updates_pl->lookup() : PT->lookupInt("WORLD_CHEMOTAXIS-brain_updates");

    //Set up data columns.
    aveFileColumns.clear();
    aveFileColumns.push_back("x_displacement");
    aveFileColumns.push_back("y_displacement");
    aveFileColumns.push_back("concentration_sum");

    //Print some stuff so the user knows what is going on.
    std::cout << "#########################" << "\n";
    std::cout << "Running ChemotaxisWorld." << "\n";
    std::cout << "Cycles to simulate is: " << eval_ticks << "\n";
    std::cout << "Brain updates per cycle is: " << brain_updates << "\n";
    std::cout << "Organism speed is: " << speed << "\n";
    std::cout << "Rotational diffusion coefficient is: " << rot_diff_coeff << "\n";
    std::cout << "Using point sources: " << point_source << "\n";
    std::cout << "Using a variable environment: " << environment_variability << "\n";
    std::cout << "Using a linear gradient: " << use_lin_gradient << "\n";
    std::cout << "Slope is: " << slope << "\n";
    std::cout << "Binary sensor is: " << use_bit_sensor << "\n";
    std::cout << "#########################" << std::endl;

}

//Should return a linear concentration unless concentration would be negative,
//in which case return 0.
inline double get_conc_linear(const double &x, const double &slope, const double &base) {
    double conc = (x * slope) + base;
  return ((conc > 0) ? conc : 0);
}

//Should return an exponential concentration, unless it's negative, in which case return 0.
inline double get_conc_exp(const double &x, const double &k, const double &base) {
  double conc = exp(k * x) + base;
  return ((conc > 0) ? conc : 0);
}

//Set up a "spot" of attractant. Can have exp or linear gradient.
//Calculate the distance from the point to the organism. Then feed that into the
//Linear/exp gradients with flipped sign (closer = higher score)
//Large spot big value: linear = big base + small slope. Exp: medium (~20?) + less than one slope
//Small spot, but high value = big base + large slope. Exp: medium (~20?) + larger slope (~1)
//Large spot low value = small base + very small slope. Exp: low base + lass than one slope
//Small spot low value = small base + large slope. Exp: low base + larger slope (~1)
//X,Y are organism position, point_x, point_y are the spot position.
inline double get_conc_point(const double &x, const double &y, const double &point_x,\
   const double &point_y, const double &slope, const double &base, const bool &lin_grad){
  double distance = std::sqrt(pow((point_x-x),2) + pow((point_y-y),2));
  return((lin_grad) ? get_conc_linear(distance, -slope, base) : get_conc_exp(distance, -slope, base));
  }

//Use the cell's x position, y position, angle theta, and the cell's speed
//to calculate the new position. pos_vec has the form [x,y,theta,speed]
inline void run(std::vector<double> &pos_vec) {
  pos_vec[0] += (pos_vec[3] * cos(pos_vec[2]));
  pos_vec[1] += (pos_vec[3] * sin(pos_vec[2]));
  return;
}

//The cell tumbles; resulting in a new theta between 0 and 2*pi
//Potentially add in ability to bias later. (This is why it's a f(x))
inline double tumble() {
  return (Random::getDouble(2 * M_PI));
}

//Causes random changes in theta which depend on the diffusion coefficient.
inline void rot_diffuse(double& theta, const double &diff_coeff) {
  theta += (Random::getNormal(0,1) * diff_coeff);
  return;
}

void ChemotaxisWorld::runWorldSolo(std::shared_ptr<Organism> org, bool analyse, bool visualize, bool debug) {
  //Should be able to thread and run concurrently at the loss of reproducibility.
  //Starting orientation is random between 0 and 2pi.
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
  if(point_source){
    concentration_hist.push_back(get_conc_point(pos_vec[0], pos_vec[1], spot_x, spot_y, slope, base, use_lin_gradient));
  }
  else{
    concentration_hist.push_back((use_lin_gradient) ? \
    get_conc_linear(pos_vec[0], slope, base) : get_conc_exp(pos_vec[0], slope, base));
  }
  pos_hist.push_back(pos_vec);

  //Initialize the cell.
  bool is_tumbling;
  double tumble_bias;
  double concentration;
  double delta;
  double accumulator;
  int multiplier;
  int num_ones;

  //Sanity check reset the brain. Shouldn't need to do this, but couldn't hurt.
  org->brain->resetBrain();

  //Add in environmental_variability, if enabled.
  if (environment_variability) {
    //Re-grab the initial values to modify.
    rot_diff_coeff = (PT == nullptr) ? rot_diff_coeff_pl->lookup() :\
    PT->lookupDouble("WORLD_CHEMOTAXIS-rot_diff_coeff");
    slope = (PT == nullptr) ? slope_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-slope");
    base = (PT == nullptr) ? base_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-base");

    //Add in the specified amount of noise.
    rot_diff_coeff += Random::getDouble(variability_rot_diff);
    slope += Random::getDouble(variability_slope);
    base += Random::getDouble(variability_base);

    //Catch negatives and just set them to zero. You could maybe skip this, but let's do it for now
    //in case someone sets negative values.
    if (rot_diff_coeff < 0) {rot_diff_coeff = 0;}
    if (slope < 0) {slope = 0;}
    if (base < 0) {base = 0;}
  } // End environmental variability

  //Evaluate the organism eval_ticks times.
  for (int t = 0; t != eval_ticks; ++t) {
    //Get concentration at the org's location using the appropriate method.
    if(point_source){
      concentration = get_conc_point(pos_vec[0], pos_vec[1], spot_x, spot_y, slope, base, use_lin_gradient);
    }
    else{
      concentration = (use_lin_gradient) ? \
      get_conc_linear(pos_vec[0], slope, base) : get_conc_exp(pos_vec[0], slope, base);
    }

    //Use the bit sensor
    if (use_bit_sensor) {
      delta = concentration - concentration_hist.back();
      delta_hist.push_back(delta);
      concentration_hist.push_back(concentration);
      if (delta > 0){
        org->brain->setInput(0,0);
      }
      else{
        org->brain->setInput(0,1);
      }
    }
    //or the magnitude sensor
    else {
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
      //and shift left by multiplier%32 on x86.
      if (multiplier > 31) {
        multiplier = 31;
      }
      //In case logic gates can output negative numbers.
      if (multiplier < 0) {
        multiplier = 0;
      }
      multiplier_hist.push_back(multiplier);

      //Delta is the relative concentration change between now and the last sample.
      delta = (concentration_hist.back() == 0) ? \
      (concentration/0.000001) : (concentration - concentration_hist.back()) / std::abs(concentration_hist.back());
      concentration_hist.push_back(concentration);
      delta_hist.push_back(delta);
      //Lshift by multiplier to allow org to control sensitivity to delta.
      delta *= (1 << multiplier);

      //Convert the delta to the unary output. A delta of 0 should produce 8 1s and 8 0s
      num_ones = std::round(delta) + 8;
      //Truncate at 0 and 16 ones.
      if (num_ones < 0) {
        num_ones = 0;
      }
      if (num_ones > 16) {
        num_ones = 16;
      }
      ones_hist.push_back(num_ones);

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
    }// end magnitude sensor
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

    //Do the tumbling and running.
    is_tumbling = Random::P(tumble_bias);
    if (is_tumbling) {
      pos_vec[2] = tumble();
      pos_hist.push_back(pos_vec);
    }
    else {
      run(pos_vec);
      rot_diffuse(pos_vec[2], rot_diff_coeff);
      pos_hist.push_back(pos_vec);
    }
  }//end eval loop

  //Finished simulating the organism, so score based on x movement and record some stats.
  //Pick an eval method
  if (point_source) {
    org->score = std::accumulate(concentration_hist.cbegin(), concentration_hist.cend(), 0.0);
    org->dataMap.Append("concentration_sum", \
    (double) std::accumulate(concentration_hist.cbegin(), concentration_hist.cend(), 0.0));
    org->dataMap.Append("x_displacement", (double) pos_vec[0]);
    org->dataMap.Append("y_displacement", (double) pos_vec[1]);
    org->dataMap.setOutputBehavior("concentration_sum", DataMap::AVE);
    org->dataMap.setOutputBehavior("x_displacement", DataMap::AVE);
    org->dataMap.setOutputBehavior("y_displacement", DataMap::AVE);
  }
  else {
    org->score = pos_vec[0];
    org->dataMap.Append("concentration_sum", \
    (double) std::accumulate(concentration_hist.cbegin(), concentration_hist.cend(), 0.0));
    org->dataMap.Append("x_displacement", (double) pos_vec[0]);
    org->dataMap.Append("y_displacement", (double) pos_vec[1]);
    org->dataMap.setOutputBehavior("concentration_sum", DataMap::AVE);
    org->dataMap.setOutputBehavior("x_displacement", DataMap::AVE);
    org->dataMap.setOutputBehavior("y_displacement", DataMap::AVE);
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

}//End of RunWorldSolo fn

//The sensor will have 16 bits describing whether concentration is increasing or decreasing.
int ChemotaxisWorld::requiredInputs() {
  return 16;
}

//The first 16 bits of output are the tumble bias (1/17-17/17) while the next 8
//are the exponent of the sensitivity multiplier.
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
