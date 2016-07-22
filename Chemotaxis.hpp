#ifndef __BasicMarkovBrainTemplate__Chemotaxis
#define __BasicMarkovBrainTemplate__Chemotaxis
#define _USE_MATH_DEFINES //I think this is okay; but if things start acting weird look here. Maybe needs
//its own header guard?

#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include <memory>
#include <cstdint>
#include <cmath>

class Chemotaxis: public AbstractWorld {
public:
  virtual void runWorldSolo(std::shared_ptr<Organism> org, bool analyse, bool visualize, bool debug) override;

  //Speed cannot be zero; it is the denominator calculating dy in the position function
  bool use_lin_gradient = true;
  double speed = 1;
  double slope = 0.1; //Slope is "m" if linear, "k" if exponential.
  double base = 1 //Starting concentration
  int eval_ticks = 5000;

//orientation?

virtual int requiredInputs() override;
virtual int requiredOutputs() override;
virtual int maxOrgsAllowed() override;
virtual int minOrgsAllowed() override;

};
