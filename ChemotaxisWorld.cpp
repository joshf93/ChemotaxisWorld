#include "ChemotaxisWorld.h"

//Initialize the parameter list.
std::shared_ptr<ParameterLink<bool>> ChemotaxisWorld::use_lin_gradient_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-use_lin_gradient", true, "Create a linear attractant gradient. Otherwise, it will be exponential.");
std::shared_ptr<ParameterLink<bool>> ChemotaxisWorld::clear_outputs_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-clear_outputs", true, "Reset output nodes between brain updates.");
std::shared_ptr<ParameterLink<bool>> ChemotaxisWorld::environment_variability_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-environment_variability", false, "Cause the environment to vary between replicate runs.");
std::shared_ptr<ParameterLink<bool>> ChemotaxisWorld::use_bit_sensor_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-use_bit_sensor", false, "Use the single bit sensor.");
std::shared_ptr<ParameterLink<bool>> ChemotaxisWorld::point_source_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-point_sources", false, "Use 'spot' sources of attractant. Still follows your choice of linear/exp gradient.");
std::shared_ptr<ParameterLink<bool>> ChemotaxisWorld::matrix_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-Matrix_mode", false, "Simulate a matrix being present in the word.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::rot_diff_coeff_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-rot_diff_coeff", 0.016, "Rotational diffusion constant.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::speed_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-speed", 1.0, "Magnitude of cell movement per tick.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::slope_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-slope", 10.0, "Rate of gradient increase. m if linear, k if exponential.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::base_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-base", 255.0, "Base concentration at x = 0.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::variability_slope_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-variability_slope", 1.0, "Slope will increase by as much as this number.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::variability_base_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-variability_base", 1.0, "Base will increase by as much as this number.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::variability_rot_diff_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-variability_rot_diff", 0.01, "Rotational diffusion constant will increase by as much as this number.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::spot_x_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-spot_x", 300.0, "The x coordinate of the attractant spot.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::spot_y_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-spot_y", 300.0, "The y coordinate of the attractant spot.");
std::shared_ptr<ParameterLink<double>> ChemotaxisWorld::radius_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-radius", 2.0, "The radius of the organism.");
std::shared_ptr<ParameterLink<int>> ChemotaxisWorld::num_points_pl = Parameters::register_parameter("WORLD_CHEMOTAXIS-num_points", 300, "The number of points to insert if a matrix is present.");
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
    use_matrix = (PT == nullptr) ? matrix_pl->lookup() : PT->lookupBool("WORLD_CHEMOTAXIS-Matrix_mode");
    rot_diff_coeff = (PT == nullptr) ? rot_diff_coeff_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-rot_diff_coeff");
    speed = (PT == nullptr) ? speed_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-speed");
    slope = (PT == nullptr) ? slope_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-slope");
    base = (PT == nullptr) ? base_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-base");
    variability_slope = (PT == nullptr) ? variability_slope_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-slope");
    variability_base = (PT == nullptr) ? variability_base_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-base");
    variability_rot_diff = (PT == nullptr) ? variability_rot_diff_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-rot_diff_coeff");
    spot_x = (PT == nullptr) ? spot_x_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-spot_x");
    spot_y = (PT == nullptr) ? spot_y_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-spot_y");
    radius = (PT == nullptr) ? radius_pl->lookup() : PT->lookupDouble("WORLD_CHEMOTAXIS-radius");
    num_points = (PT == nullptr) ? num_points_pl->lookup() : PT->lookupInt("WORLD_CHEMOTAXIS-num_points");
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

double fast_atan2(const double &y, const double &x) {
    //Fast approximation of atan2, see:
    //http://www.iro.umontreal.ca/~mignotte/IFT2425/Documents/EfficientApproximationArctgFunction.pdf
    //For 100M iterations on my laptop:
    // atan2(any angle) ~ 20s
    // fast_atan2(any angle) ~13s
    // fast_atan2(small angles, same sign) ~3s, much easier on branch predictor
    // fast_atan2(small angles, random sign) ~5.5s
    //Error never exceeds 0.0038, just like the paper says. It's usually smaller
    //Than that.
    double result;
    double ratio = y/x;

    //If the ratio is outside 1, this approx can't be used. Also bow out for exact
    //zeroes.
    if (ratio > 1 || y == 0.0  || x == 0.0){
      return atan2(y, x);
    }

    //Otherwise, use eq 7 from the paper. This is the octant section.
    //Signs of x and y.
    int s_x = (x >= 0) ? 1 : -1;
    int s_y = (y >= 0) ? 1 : -1;

    if ((s_x == 1) && (fabs(ratio) < (M_PI/4))){
      //Oct I and VIII
      result = ratio * (1.0584-s_y*0.273*ratio);
    }
    else if ((s_y == 1) && (fabs(ratio) > (M_PI/4))){
      //Oct II and III
      result = (M_PI/2) - (1/ratio)*(1.0584-s_x*0.273*(1/ratio));
    }
    else if ((s_x == -1) && (fabs(ratio) < (M_PI/4))){
      //Oct IV and V
      result = s_y*M_PI + ratio*(1.0584+s_y*0.273*ratio);
    }
    else {
      //Oct VI and VII
      result = -(M_PI/2) - (1/ratio)*(1.0584+s_x*0.273*(1/ratio));
    }

    return result;
  }

  //Optimized version. Will calculate the run even if no collision occurs.
point_vector check_collision(const point_vector &pos_vec, \
    const std::vector<point_vector> &matrix_vec, const double &radius){
    //Run. If no collisions, equiv to normal run. If collisions, moves with a
    //square collision front to the nearest collision.
    point_vector move_vec {0,0,pos_vec.theta,pos_vec.magnitude};
    double cp_x;
    double cp_y;
    double cp_dist;
    double psi; //Theta converted from 0,2pi to -pi,pi for angle dist measures.
    double phi;
    double move_dist;
    double r_dist;
    //Calculate the move as if no collisions happened.
    move_vec.x = pos_vec.x + pos_vec.magnitude * cos(pos_vec.theta);
    move_vec.y = pos_vec.y + pos_vec.magnitude * sin(pos_vec.theta);

    //Tumbles can change theta to the range [0,2pi]. In order to do our next trick,
    //theta must not be greater than 2pi. This step should suffice, as rotational
    //diffusion will only take the organism so far; it shouldn't pass 4pi.
    psi = (pos_vec.theta > 2*M_PI) ? pos_vec.theta - 2*M_PI : pos_vec.theta;
    //Doing this makes angle comparisons a bit easier considering the output of
    //atan2.
    psi = (psi > M_PI) ? psi - 2*M_PI : psi;

    for (auto pt = matrix_vec.begin(); pt < matrix_vec.end(); pt++) {
      //Go through each point and check for a collision.
      //Starting from c and moving to p.
      cp_x = pt->x - pos_vec.x;
      cp_y = pt->y - pos_vec.y;

      //Filter out points that are too far away. This is an extra step, but is
      //worthwhile for the following reasons: the full suite of calcs has
      //some slow calculations (atan2 in particular). The second is that the
      //operations needed to exclude these points are very cheap; two additions, two
      //abs, and three compares/logicals, each of which are extremely quick (1-3cycles).
      //Finally, the branch predictor should handle it well, since 95% of points will
      //take the early exit opportunity.
      if(fabs(cp_x) > pos_vec.magnitude + radius || fabs(cp_y) > pos_vec.magnitude + radius) {
        continue; //Impossible to collide with it. Skip to the next point.
      }

      //Calculate dist from point to movement line.
      cp_dist = sqrt((cp_x*cp_x)+(cp_y*cp_y));
      phi =  fast_atan2(cp_y, cp_x) - psi; //Custom atan approximation
      move_dist = cos(phi) * cp_dist;
      r_dist = sin(phi) * cp_dist;
      if ((fabs(r_dist) < radius) && (cp_dist < (pos_vec.magnitude + radius)) && (fabs(phi) < M_PI/2) && (move_dist < pos_vec.magnitude)) {
        //move_vec.magnitude = move_dist;
        move_vec.x = cp_x;
        move_vec.y = cp_y;
        //Randomly orient after collision. Doesn't matter if we do it multiple times,
        //which should be relatively rare anyway.
        move_vec.theta = Random::getDouble(2 * M_PI);
        //std::cout << "COLLISION DETECTED:" << cp_x << "," << cp_y << std::endl;
      }
    }
    return (move_vec);
  }

//Use the cell's x position, y position, angle theta, and the cell's speed
//to calculate the new position. pos_vec has the form [x,y,theta,speed]
inline point_vector run(const point_vector &pos_vec) {
  point_vector new_pos {0, 0, pos_vec.theta, pos_vec.magnitude};
  new_pos.x = pos_vec.x + (pos_vec.magnitude * cos(pos_vec.theta));
  new_pos.y = pos_vec.y + (pos_vec.magnitude * sin(pos_vec.theta));
  return (new_pos);
}

//The cell tumbles; resulting in a new theta between 0 and 2*pi
//Potentially add in ability to bias later. (This is why it's a f(x))
inline double tumble() {
  return (Random::getDouble(2 * M_PI));
}

//Causes random changes in theta which depend on the diffusion coefficient.
inline double rot_diffuse(double& theta, const double &diff_coeff) {
  return (theta + (Random::getNormal(0,1) * diff_coeff));
}

void ChemotaxisWorld::runWorldSolo(std::shared_ptr<Organism> org, bool analyse, bool visualize, bool debug) {
  //Starting orientation is random between 0 and 2pi.
  point_vector pos_vec{0.0, 0.0, Random::getDouble(0,2 * M_PI), speed};
  point_vector spot_vec{spot_x, spot_y, 0, 0};
  std::vector<point_vector> pos_hist;
  std::vector<double> tumble_hist;
  std::vector<double> concentration_hist;
  std::vector<double> delta_hist;
  std::vector<int> multiplier_hist;
  std::vector<int> ones_hist;
  tumble_hist.reserve(eval_ticks+2);
  concentration_hist.reserve(eval_ticks+2);
  multiplier_hist.reserve(eval_ticks+2);
  delta_hist.reserve(eval_ticks+2);

  //Initialize conc_hist to keep things simpler when calculating delta.
  if(point_source){
    concentration_hist.push_back(get_conc_point(pos_vec.x, pos_vec.y, spot_vec.x, spot_vec.y, slope, base, use_lin_gradient));
  }
  else{
    concentration_hist.push_back((use_lin_gradient) ? \
    get_conc_linear(pos_vec.x, slope, base) : get_conc_exp(pos_vec.x, slope, base));
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

  //Matrix initialization
  //Only need to put circles in the range of positions that can be reached
  std::vector <point_vector> matrix_vec;
  double pt_lim;
  double pt_x;
  double pt_y;

  //Generate the matrix, if applicable
  if (use_matrix) {
    pt_lim = speed*eval_ticks;
    for (int n = 0; n < num_points; n++) {
      pt_x = Random::getDouble(0, pt_lim);
      pt_y = Random::getDouble(-pt_lim, pt_lim);
      matrix_vec.push_back(point_vector {pt_x,pt_y,0,0});
      }
    } //End matrix generation


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
      concentration = get_conc_point(pos_vec.x, pos_vec.y, spot_vec.x, spot_vec.y, slope, base, use_lin_gradient);
    }
    else{
      concentration = (use_lin_gradient) ? \
      get_conc_linear(pos_vec.x, slope, base) : get_conc_exp(pos_vec.x, slope, base);
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
      pos_vec.theta = tumble();
      pos_hist.push_back(pos_vec);
    }
    else {
      pos_vec = (use_matrix) ? check_collision(pos_vec, matrix_vec, radius) : run(pos_vec);
      pos_vec.theta = rot_diffuse(pos_vec.theta, rot_diff_coeff);
      pos_hist.push_back(pos_vec);
    }
    if (pos_vec.x > (speed*eval_ticks) || pos_vec.y > (speed*eval_ticks)){
      std::cout << "Went farther than possible. Bug somewhere." << "\n";
      std::cout << "Pos vec:" << pos_vec.x << "," << pos_vec.y << "," << pos_vec.theta << "," << pos_vec.magnitude << std::endl;
    }
  }//end eval loop

  //Finished simulating the organism, so score based on x movement and record some stats.
  //Pick an eval method. Look at TestWorld and BerryWorld to get an idea of how
  //the dataMaps work.
  if (point_source) {
    org->score = std::accumulate(concentration_hist.cbegin(), concentration_hist.cend(), 0.0);
    org->dataMap.Append("concentration_sum", \
    (double) std::accumulate(concentration_hist.cbegin(), concentration_hist.cend(), 0.0));
    org->dataMap.Append("x_displacement", (double) pos_vec.x);
    org->dataMap.Append("y_displacement", (double) pos_vec.y);
    org->dataMap.setOutputBehavior("concentration_sum", DataMap::AVE);
    org->dataMap.setOutputBehavior("x_displacement", DataMap::AVE);
    org->dataMap.setOutputBehavior("y_displacement", DataMap::AVE);
  }
  else {
    org->score = pos_vec.x;
    org->dataMap.Append("concentration_sum", \
    (double) std::accumulate(concentration_hist.cbegin(), concentration_hist.cend(), 0.0));
    org->dataMap.Append("x_displacement", (double) pos_vec.x);
    org->dataMap.Append("y_displacement", (double) pos_vec.y);
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
      out_file << sample.x << "," << sample.y << "," << sample.theta << "," << '\n';
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
