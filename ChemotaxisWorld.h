#ifndef __BasicMarkovBrainTemplate__Chemotaxis
#define __BasicMarkovBrainTemplate__Chemotaxis
#define _USE_MATH_DEFINES //I think this is okay; but if things start acting weird look here. Maybe needs
//its own header guard?

#include "../AbstractWorld.h"

#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <sstream>

class ChemotaxisWorld: public AbstractWorld {

public:

  ChemotaxisWorld(std::shared_ptr<ParametersTable> _PT = nullptr);
  ~ChemotaxisWorld() = default;
  virtual void runWorldSolo(std::shared_ptr<Organism> org, bool analyse, bool visualize, bool debug) override;


  static std::shared_ptr<ParameterLink<bool>> use_lin_gradient_pl;
  static std::shared_ptr<ParameterLink<bool>> clear_outputs_pl;
  static std::shared_ptr<ParameterLink<double>> rot_diff_coeff_pl;
  static std::shared_ptr<ParameterLink<double>> speed_pl;
  static std::shared_ptr<ParameterLink<double>> slope_pl;
  static std::shared_ptr<ParameterLink<double>> base_pl;
  static std::shared_ptr<ParameterLink<int>> eval_ticks_pl;




  bool use_lin_gradient;
  bool clear_outputs;
  double rot_diff_coeff;
  double speed;
  double slope; //Slope is "m" if linear, "k" if exponential.
  double base; //Starting concentration
  int eval_ticks;



virtual int requiredInputs() override;
virtual int requiredOutputs() override;
virtual int maxOrgsAllowed() override;
virtual int minOrgsAllowed() override;

};
#endif
