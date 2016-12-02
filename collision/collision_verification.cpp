#include <vector>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>
#include <numeric>
#include <chrono>
#define _USE_MATH_DEFINES


struct point_vector {
  double x;
  double y;
  double theta;
  double magnitude;
};

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


//Optimized version
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

  move_vec.x = move_vec.magnitude * cos(move_vec.theta);
  move_vec.y = move_vec.magnitude * sin(move_vec.theta);
  psi = (pos_vec.theta > M_PI) ? (pos_vec.theta - 2*M_PI) : pos_vec.theta;

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
    if ((fabs(r_dist) < radius) && (cp_dist < (pos_vec.magnitude + radius)) && (fabs(phi) < M_PI/2) && (move_dist < move_vec.magnitude)) {
      move_vec.magnitude = move_dist;
      move_vec.x = cp_x;
      move_vec.y = cp_y;
    }
  }
  return (move_vec);
}

int main(){
  point_vector vec_p {0,0,0,10};
  double rad = 1.0;
  point_vector result;
  int temp;
  int temp2;
  std::vector<point_vector> mat_vec;// {{2,0,0,0}, {2,2,0,0}, {0,2,0,0}, {-2,2,0,0}, {-2,0,0,0}, {-2,-2,0,0}, {0,-2,0,0}, {2,-2,0,0}};
  std::vector<double> atan_vec;
  std::vector<point_vector> result_vec;
  srand(std::clock());

  //std::cout << "Enter theta to test:" << std::endl;
  //std::cin >> user_theta;
  //vec_p.theta = user_theta;
  for (int g = 0; g < 10; ++g){
    temp = (rand()%2 == 0) ? -1 : 1;
    temp2 = (rand()%2 == 0) ? -1 : 1;
    mat_vec.push_back({temp*rand()%8, temp2*rand()%8,0,0});
  }


  double t0 = std::clock();
  for (int j = 0; j < 628; ++j){
    vec_p.theta = j/100.0;
    result = check_collision(vec_p, mat_vec, rad);
    result_vec.push_back(result);

  }

  double t1 = std::clock();
  std::cout << "Time: " << (t1-t0)/1000000 << "\n";
  std::ofstream log_file("verification_collision_data.csv");
  log_file << "X,Y,T,M" << "\n";
  for (auto n : result_vec){
    log_file << n.x << "," << n.y << "," << n.theta << "," << n.magnitude << "\n";
  }
  log_file << std::endl;

  /*
  //fast_atan2 test
  long ti = std::clock();
  for (int n = 0; n < 100000000; ++n){
    temp = (rand()%2 == 0) ? -1 : 1;
    temp2 = (rand()%2 == 0) ? -1 : 1;
    result_vec.push_back(atan2(1.0*temp,double(rand()%10)+1.0)*temp2);
  }
  long tf = std::clock();
  std::cout << "Sum is:" << std::accumulate(result_vec.begin(), result_vec.end(), 0.0) << std::endl;
  std::cout << "Time was:" << (tf-ti)/double(1000000) << std::endl;
  */
  return 0;
}
