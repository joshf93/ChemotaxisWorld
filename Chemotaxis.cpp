#include "Chemotaxis.hpp"
//This comment is to test git.

//Should return a linear concentration unless it's negative, in which case return 0.
//Will not attempt to convert negative value, and will return 0 instead.
uint32_t get_conc_linear(const double &x, const double &slope, const double &base){
    double conc = x * slope;
  return ((conc > 0) ? conc : 0);
}

//Should return an exponential concentration, unless it's negative, in which case return 0.
//Will not attempt to convert negative value, and will return 0 instead.
uint32_t get_conc_exp(const double &x, const double &k, const double &base){
  double conc = exp(k * x) + base;
  return ((conc > 0) ? conc : 0)
}

//Use the cell's x position, y position, angle theta, and the cell's
//speed to calculate the new position. Void because it works directly on the vector by ref.
//pos_vec has the form [x,y,theta,speed]
void run(std::vector<double> &pos_vec){
  pos_vec[0] += (pos_vec[3] * cos(pos_vec[2]));
  pos_vec[1] += (pos_vec[3] * sin(pos_vec[2]));
  return;
}

//The cell tumbles; resulting in a new theta between 0 and 2*pi
//Potentially add in ability to bias later.
inline double tumble(){
  return (Random::getDouble(2 * M_PI));
}

//Take a vector of bits and convert it into a probability of tumbling.
//The "endianness" of how the brain gives back outputs shouldn't matter.
//The organism will 'learn' it. Will only take 8 bits for now; make generic size if needed.
//Takes vector of ints, which should be 1 or 0, multiply each position by 2**n and add
//the result to an accumulator to convert to int. Then take accumulator+1 (so bias can't be 0)
//and divide by 256 to get the resulting tumble bias. If somehow this a bottleneck,
//convert pow(2,idx) to left shifts.
double bit_to_prob(const std::vector<int> &bit_vec){
  int accumulator = 0;
  for(int idx = 0; idx != 8; ++idx){
    accumulator += (bit_vec[idx]*pow(2,idx));
  }
  return ((accumulator+1)/256);
}

//Calculate total displacement, x displacement, etc. stats for analyzing organism
//After it has completed evaluation.
//Return vector format: [overall_displacement, x_displacement, ]
vector<double> end_stats(const std::vector<double> &pos_vec){


}

//Determine whether the nth bit of an uint32_t is set. Will fail if idx > 31.
inline int is_set(const uint32_t &num, const int &idx){
  if(idx > 31){
    std::cout << "Attempted to shift left by more than 31 while calculating brain input." << std::endl;
    exit(1);
  }
  return ((num & (1 << idx)) >> idx);
}


void Chemotaxis::runWorldSolo(std::shared_ptr<Organism> org, bool analyse, bool visualize, bool debug){
  //Starting orientation should be random. Make it an option later.
  //Also starts the organism on the origin for now. Make it an option later?

  std::vector<double> pos_vec(0, 0, Random::getDouble(2 * M_PI), speed);
  std::vector<std::vector<double>> pos_hist;

  //Initialize the cell to be tumbling or not; 50/50 chance hardcoded for now.
  //Make this an option later?
  bool is_tumbling = Random::P(0.5);
  //Start tumble bias at 0.5. Option later. Concentration to 0, it'll be changed first tick anyway.
  double tumble_bias = 0.5;
  uint32_t concentration = 0;


  //Evaluate the organism for eval_ticks ticks.
  for(int t = 0; t != eval_ticks; ++t){
    //Update the brain first.
    if(use_lin_gradient){
      concentration = get_conc_linear(pos_vec[0], slope, base);
    }
    else{
      concentration = get_conc_exp(pos_vec[0], slope, base);
    }

    //Our brain position index for the tick.
    int n = 0;
    //Go through and send each bit in concentration to the brain input.
    while(n != 33){
      org->brain->setInput(n,is_set(concentration,n));
      ++n;
    }

    //Update the brain after clearing outputs (if selected)
    if (clearOutputs){
      org->brain->resetOutputs;
    }
    org->brain->update();

    //Read the output and convert to tumble bias.
    //Will involve a bitset, convert to ulong, double divide (resulting long +1) by 256.
    //+1 to prevent tumble bias from ever being 0 ((0+1)//256 != 0).
    std::vector <int> output_vector;
    for (int n = 0; n != 9; ++n){
      output_vector.push_back(org->brain->readOutput(n));
    }
    tumble_bias = bit_to_prob(output_vector);


  }





}


//This will vary depending on what, exactly, the sensor is reading, and how it is being read.
//For now, let's just go with a 32 bit gradient; "concentration" values will range from 0 to (2^32)-1.
//Maybe set inputs and outputs to be user modified at a later point. For now, leave them fixed.
int Chemotaxis::requiredInputs() {
  return 32;
}

//Only affects tumble bias. Very initial value: 8 bits; tumble bias will range from 1/256 to 256/256;
int Chemotaxis::requiredOutputs() {
  return 8;
}

int Chemotaxis::MaxOrgsAllowed() {
  return 1;
}

int Chemotaxis::minOrgsAllowed() {
  return 1;
}
