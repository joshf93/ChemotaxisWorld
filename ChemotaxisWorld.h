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

class ChemotaxisWorld: public AbstractWorld {

public:

  ChemotaxisWorld(std::shared_ptr<ParametersTable> _PT = nullptr);
  ~ChemotaxisWorld() = default;
  virtual void runWorldSolo(std::shared_ptr<Organism> org, bool analyse, bool visualize, bool debug) override;


  bool use_lin_gradient = true;
  bool clear_outputs = true;
  double speed = 1.0;
  double slope = 1; //Slope is "m" if linear, "k" if exponential.
  double base = 1; //Starting concentration
  int eval_ticks = 1000;



virtual int requiredInputs() override;
virtual int requiredOutputs() override;
virtual int maxOrgsAllowed() override;
virtual int minOrgsAllowed() override;

};
#endif
