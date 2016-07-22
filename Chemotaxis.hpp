#ifndef __BasicMarkovBrainTemplate__Chemotaxis
#define __BasicMarkovBrainTemplate__Chemotaxis

#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include <memory>
#include <cstdint>
#include <cmath>
#include <bitset>

class Chemotaxis: public AbstractWorld {
public:
  virtual void runWorldSolo(shared_ptr<Organism> org, bool analyse, bool visualize, bool debug) override;

  //Speed cannot be zero; it is the denominator calculating dy in the position function
  bool use_lin_gradient = true;
  double speed = 1;
  double slope = 0.1; //Slope is "m" for if linear, "k" if exponential.
  int eval_ticks = 5000;

//orientation?

virtual int requiredInputs() override;
virtual int requiredOutputs() override;
virtual int maxOrgsAllowed() override;
virtual int minOrgsAllowed() override;

};
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
