#include <vector>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>
#include <omp.h>
#define _USE_MATH_DEFINES


struct point_vector {
  double x;
  double y;
  double theta;
  double magnitude;
};

bool check_collision(const point_vector &pos_vec, \
  const std::vector<point_vector> &matrix_vec, const double &radius){
  //Check to see if a collision occured.
  //Calc collision rectangle
  double alpha = -((M_PI/2)-pos_vec.theta);
  point_vector point_a {pos_vec.x + radius*cos(alpha), pos_vec.y + radius*sin(alpha), 0, 0};
  point_vector vec_ab {point_a.x, point_a.y, M_PI+alpha, radius*2};
  point_vector vec_ad {point_a.x, point_a.y, pos_vec.theta, pos_vec.magnitude};
  //point_vector rec_c {(rec_b.x + cos(pos_vec.theta)), rec_b.y + sin(pos_vec.theta), 0, 0};
  //Method used does not require the fourth point.


  /*See if any of the points are in the rectangle. Thanks to:
  http://math.stackexchange.com/questions/190111/how-to-check-if-a-point-is-inside-a-rectangle
  Raymond Manzoni (http://math.stackexchange.com/users/21783/raymond-manzoni),
  How to check if a point is inside a rectangle?, URL (version: 2012-09-03):
  http://math.stackexchange.com/q/190373) */

  // http://bisqwit.iki.fi/story/howto/openmp/#IntroductionToOpenmpInC
  bool in_rec = false;
  {
  #pragma omp parallel for
  for (auto pt = matrix_vec.begin(); pt < matrix_vec.end(); pt++) {
    // Calculate the A->P vector
    double ap_x = pt->x - point_a.x;
    double ap_y = pt->y - point_a.y;
    point_vector vec_ap {ap_x, ap_y, atan2(ap_y, ap_x), sqrt(pow(ap_x,2)+pow(ap_y,2))};
    double ap_ab_prod = vec_ap.magnitude * vec_ab.magnitude * cos(vec_ap.theta-vec_ab.theta);
    double ab_squared = vec_ab.magnitude * vec_ab.magnitude;
    double ap_ad_prod = vec_ap.magnitude * vec_ad.magnitude * cos(vec_ap.theta - vec_ad.theta);
    double ad_squared = vec_ad.magnitude * vec_ad.magnitude;

    if ((ap_ab_prod < ab_squared) && (ap_ad_prod < ad_squared)){
      //std::cout << "Found one:" << pt->x << "," << pt->y << "\n";
      //std::cout << "ap_x:" << ap_x << ", ap_y:" << ap_y << "\n";
      //std::cout << "point_ax:" << point_a.x << ", point_ay:" << point_a.y << std::endl;
      #pragma omp critical
      {
        in_rec = true;
      }
      #pragma omp cancel for
    }
  }
}//end OpenMP
  return (in_rec);
}


int main(){
  point_vector vec_p {0,0,0,6};
  double rad = 5;
  std::vector<point_vector> mat_vec {};
  srand(std::clock());

  bool cancel_status = omp_get_cancellation();
  if (not cancel_status){
    std::cout << R"(OpenMP cancellation is not enabled, and will slow matrix
      collision calculations considerably. To enable, enter the following in
      your command line: export OMP_CANCELLATION=true)" << std::endl;
  }

  for (int n = 0; n < 100000000; n++){
    mat_vec.push_back({rand()%100, rand()%100, 0, 0});
  }
  double t0 = std::clock();
  bool result = check_collision(vec_p, mat_vec, rad);
  double t1 = std::clock();
  std::cout << "Collions: " << result << "\n";
  std::cout << "Time: " << t1-t0 << "\n";
  std::cout << std::endl;


  return 0;
}
