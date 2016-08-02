#ifndef __BasicMarkovBrainTemplate__Chemotaxis
#define __BasicMarkovBrainTemplate__Chemotaxis
#define _USE_MATH_DEFINES //I think this is okay; but if things start acting weird look here. Maybe needs
//its own header guard?

#include "../AbstractWorld.h"
#include <vector>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <fstream>

class ChemotaxisWorld: public AbstractWorld {
public:
  //Constructor and destructor.
  ChemotaxisWorld(std::shared_ptr<ParametersTable> _PT = nullptr);
  ~ChemotaxisWorld() = default;

  //Pointers to parameters
  static std::shared_ptr<ParameterLink<bool>> use_lin_gradient_pl;
  static std::shared_ptr<ParameterLink<bool>> clear_outputs_pl;
  static std::shared_ptr<ParameterLink<bool>> environment_variability_pl;
  static std::shared_ptr<ParameterLink<bool>> use_integral_sensor_pl;
  static std::shared_ptr<ParameterLink<double>> rot_diff_coeff_pl;
  static std::shared_ptr<ParameterLink<double>> speed_pl;
  static std::shared_ptr<ParameterLink<double>> slope_pl;
  static std::shared_ptr<ParameterLink<double>> base_pl;
  static std::shared_ptr<ParameterLink<double>> variability_slope_pl;
  static std::shared_ptr<ParameterLink<double>> variability_base_pl;
  static std::shared_ptr<ParameterLink<double>> variability_rot_diff_pl;
  //static std::shared_ptr<ParameterLink<double>> variability_speed_pl;
  static std::shared_ptr<ParameterLink<int>> eval_ticks_pl;
  static std::shared_ptr<ParameterLink<int>> brain_updates_pl;

  //Declare the parameters
  bool use_lin_gradient;
  bool clear_outputs;
  bool environment_variability;
  bool use_integral_sensor;
  double rot_diff_coeff;
  double speed;
  double slope; //Slope is "m" if linear, "k" if exponential.
  double base; //Starting concentration
  double variability_slope; //Magnifies or decreases the amount of variability
  double variability_base;
  double variability_rot_diff;
  //double variability_speed;
  int eval_ticks;
  int brain_updates;

  //Override AbstractWorld methods.
  virtual void runWorldSolo(std::shared_ptr<Organism> org, bool analyse, bool visualize, bool debug) override;
  virtual int requiredInputs() override;
  virtual int requiredOutputs() override;
  virtual int maxOrgsAllowed() override;
  virtual int minOrgsAllowed() override;
};
#endif
