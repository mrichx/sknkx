#include "sknkx_opt/sknkxUtil.hpp"

using namespace sknkx;

void SknkxTimer::tic() {
  this->tictoc_stack.push(omp_get_wtime());
}

double SknkxTimer::toc() {
  double tElapsed = (double)omp_get_wtime() - this->tictoc_stack.top();
  this->tictoc_stack.pop();
  return tElapsed;
}

// void SknkxInterp::init(std::vector<double> xi, std::vector<double> yi) {
//   // check that xi and yi have same number of elements
//   if (xi.size() != yi.size()) {
//     return;
//   } else {
//     // clear existing data
//     this->x_.clear();
//     this->y_.clear();
//     // copy data
//     for (int i=0; i<xi.size(); i++) {
//       x_.push_back(xi[i]);
//       y_.push_back(yi[i]);
//     }
//   }
// }

double SknkxInterp::getLinearInterp(std::vector<double> &xi, std::vector<double> &yi,
                                    double x) {
  // check that x_ and y_ are set
  if (xi.size() == 0) {
    return 0;
  } else {
    // assume that data is sorted (?)
    int i=0;
    while ((xi[i]<x) & (i<xi.size())) {
      i++;
    }
    // perform linear interpolation
    double data;
    if (i==0) {
      data = yi[i];
    } else if (i==xi.size()) {
      data = yi[i];
    } else {
      double dx = xi[i+1] - xi[i];
      double dy = yi[i+1] - yi[i];
      data = yi[i] + dy / dx * (x-xi[i]);
    }
    return data;
  }
}
