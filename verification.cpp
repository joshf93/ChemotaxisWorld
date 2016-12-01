#include <vector>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>
#define _USE_MATH_DEFINES


struct point_vector {
  double x;
  double y;
  double theta;
  double magnitude;
};


//Approximation atan2, see http://www.iro.umontreal.ca/~mignotte/IFT2425/Documents/EfficientApproximationArctgFunction.pdf

inline double fast_atan2(const double &ratio, const double &y, const double &x) {
  double result;
  //If the ratio is outside 1, this approx can't be used.
  if (ratio > 1){return atan2(y, x);}
  //Otherwise, use eq 7 from the paper.
  int s_x = (x >= 0) ? 1 : -1;
  int s_y = (y >= 0) ? 1 : -1;

  if (s_x == 1){
    if(fabs(ratio) < (M_PI/4)){
    result = ratio * (1.0584 - s_y*0.273*ratio);
    }
  }





}

//Optimized version
point_vector check_collision(const point_vector &pos_vec, \
  const std::vector<point_vector> &matrix_vec, const double &radius){
  //Check to see if a collision occured.
  //Calc collision rectangle
  point_vector in_rec {0,0,0,0};
  double dist = pos_vec.magnitude;

  for (auto pt = matrix_vec.begin(); pt < matrix_vec.end(); pt++) {
    //Filter out points that are too far away. This is an extra step, but is
    //worthwhile for the following reasons: the full suite of calcs has
    //some slow calculations (atan2 in particular). The second is that the
    //operations needed to exclude these points are very cheap; two additions, two
    //abs, and three compares/logicals, each of which are extremely quick (1-3cycles).
    //Finally, the branch predictor should handle it well, since 95% of points will
    //take the early exit opportunity.
    //Starting from c and moving to p.
    double cp_x;
    double cp_y;
    double cp_dist;
    double r_dist;
    double move_dist;
    double phi;
    //std::cout << "Point:" << pt->x << ", " << pt->y << std::endl;

    cp_x = pt->x - pos_vec.x;
    cp_y = pt->y - pos_vec.y;

    //std::cout << "cp_x:" << cp_x << "," << " cp_y:" << cp_y << std::endl;
    if(fabs(cp_x) > pos_vec.magnitude + radius || std::abs(cp_y) > pos_vec.magnitude + radius) {
      continue; //No possible way to hit it. Go to the next point.
    }
    //Calculate dist from point to movement line.
    cp_dist = sqrt((cp_x*cp_x)+(cp_y*cp_y));
    phi = atan2(cp_y, cp_x) - pos_vec.theta;
    move_dist = cos(phi) * cp_dist;
    r_dist = sin(phi) * cp_dist;

    //std::cout << "cp_dist:" << cp_dist << std::endl;
    //std::cout << "phi:" << phi << " move_dist:" << move_dist << " r_dist:" << r_dist << std::endl;

    if ((fabs(r_dist) < radius) && (cp_dist < (pos_vec.magnitude + radius)) && (fabs(phi) < M_PI/2)) {
      //std::cout << "Passed if, collision detected." << std::endl;
      //Figure out how far we need to move.
      dist = (move_dist < dist) ? move_dist : dist;
      //std::cout << "Collision X,Y: "<< cp_x << "," << cp_y << std::endl;
      //std::cout << "Dist is:" << dist << std::endl;
    }
    //std::cout << std::endl;
  }
  in_rec.magnitude = dist;
  return (in_rec);
}



//Square version.
/*
point_vector check_collision(const point_vector &pos_vec, \
  const std::vector<point_vector> &matrix_vec, const double &radius){
  //Check to see if a collision occured.
  //Calc collision rectangle
  double alpha = -((M_PI/2)-pos_vec.theta);
  double dist = pos_vec.magnitude;
  point_vector point_a {pos_vec.x + radius*cos(alpha), pos_vec.y + radius*sin(alpha), 0, 0};
  point_vector vec_ab {point_a.x, point_a.y, M_PI+alpha, radius*2};
  point_vector vec_ad {point_a.x, point_a.y, pos_vec.theta, pos_vec.magnitude + radius};
  //point_vector rec_c {(rec_b.x + cos(pos_vec.theta)), rec_b.y + sin(pos_vec.theta), 0, 0};
  //Method used does not require the fourth point.


  //See if any of the points are in the rectangle. Thanks to:
  //http://math.stackexchange.com/questions/190111/how-to-check-if-a-point-is-inside-a-rectangle
  //Raymond Manzoni (http://math.stackexchange.com/users/21783/raymond-manzoni),
  //How to check if a point is inside a rectangle?, URL (version: 2012-09-03):
  //http://math.stackexchange.com/q/190373)

  // http://bisqwit.iki.fi/story/howto/openmp/#IntroductionToOpenmpInC
  point_vector in_rec {0,0,0,0};
  {
  //#pragma omp parallel for
  for (auto pt = matrix_vec.begin(); pt < matrix_vec.end(); pt++) {
    //Filter out points that are too far away. This is an extra step, but is
    //worthwhile for the following reasons: the full suite of vector calcs has
    //some *very* slow calculations (atan2 in particular), and it has 4 of them
    //each of which are at least 10x slower than simple ops. The second is that the
    //operations needed to exclude these points are very cheap; two additions, two
    //abs, and three compares/logicals, each of which are extremely quick (1-3cycles).
    //Finally, the branch predictor should handle it well, since 95% of points will
    //take the early exit opportunity. You cannot just use the ap_n below, as that
    //references the A point, which is *not* the same as the cell center.
    //
    double am_x = pt->x - pos_vec.x;
    double am_y = pt->y - pos_vec.y;
    if(std::abs(am_x) > pos_vec.magnitude + radius || std::abs(am_y) > pos_vec.magnitude + radius){
      continue;
    }
    // Calculate the A->P vector
    double ap_x = pt->x - point_a.x;
    double ap_y = pt->y - point_a.y;
    point_vector vec_ap {ap_x, ap_y, atan2(ap_y, ap_x), sqrt(pow(ap_x,2)+pow(ap_y,2))};
    double ap_ab_prod = vec_ap.magnitude * vec_ab.magnitude * cos(vec_ap.theta-vec_ab.theta);
    double ab_squared = vec_ab.magnitude * vec_ab.magnitude;
    double ap_ad_prod = vec_ap.magnitude * vec_ad.magnitude * cos(vec_ap.theta-vec_ad.theta);
    double ad_squared = vec_ad.magnitude * vec_ad.magnitude;

    if ((ap_ab_prod < ab_squared) && (ap_ad_prod < ad_squared)){

      //If true, we need to calculate the distance.
      double mp_x = pt->x - pos_vec.x;
      double mp_y = pt->y - pos_vec.y;
      double mp_dist = sqrt(pow(mp_x,2)+pow(mp_y,2));
      double mp_angle = atan2(mp_y, mp_x);
      double gamma = (pos_vec.theta > mp_angle) ? pos_vec.theta - mp_angle : mp_angle - pos_vec.theta;
      double move_dist = mp_dist * cos(gamma);
      dist = (move_dist < dist) ? move_dist : dist;
      //std::cout << "Collision X,Y: "<< mp_x << "," << mp_y << std::endl;
    }
  }
}
  in_rec.magnitude = dist;
  return (in_rec);
}
*/


int main(){
  point_vector vec_p {0,0,0,10};
  double rad = 1.0;
  int sum = 0;
  double disty = 0;
  point_vector result;
  std::vector<point_vector> mat_vec; //{{2,0,0,0}, {2,2,0,0}, {0,2,0,0}, {-2,0,0,0}, {-2,-2,0,0}, {2,-2,0,0}, {-2,2,0,0}, {0,-2,0,0}};
  srand(std::clock());
  //double user_theta = 0;

  //std::cout << "Enter theta to test:" << std::endl;
  //std::cin >> user_theta;
  vec_p.theta = rand()%10;

  for (int n = 0; n < 1000; ++n){
    mat_vec.push_back({rand()%128, rand()%128, 0,0});
  }

  double t0 = std::clock();
  for (int j = 0; j < 1000000; ++j){
    result = check_collision(vec_p, mat_vec, rad);
    sum += result.magnitude;
    disty = result.magnitude;
  }

  double t1 = std::clock();
  std::cout << "Collisions: " << sum << "\n";
  std::cout << "Time: " << (t1-t0)/1000000 << "\n";
  std::cout << "Move distance: " << disty << "\n";
  std::cout << "Collision points, X: " << result.x << " Y:" << result.y << "\n";
  std::cout << std::endl;


  return 0;
}
