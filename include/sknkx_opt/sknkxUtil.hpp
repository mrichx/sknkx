#ifndef __SKNKXUTIL_HPP__
#define __SKNKXUTIL_HPP__

#include <stack>
#include <vector>
#include "omp.h"


namespace sknkx {

  class SknkxTimer{
  public:
    void tic();
    double toc();
  private:
    std::stack<double> tictoc_stack;
  };

  class SknkxInterp{
  public:
    //void init(std::vector<double> xi, std::vector<double> yi);
    double getLinearInterp(std::vector<double> &xi, std::vector<double> &yi,
                           double x);
  private:
    //std::vector<double> x_;
    //std::vector<double> y_;
  };
}

#endif
