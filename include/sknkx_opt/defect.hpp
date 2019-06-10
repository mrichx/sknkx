#ifndef __SKNKXDEFECT_HPP__
#define __SKNKXDEFECT_HPP__

#include "sknkxCore.hpp"

namespace sknkx {

  class Defect : public ConstraintFunc {
    public:
      void evaluate();
      //void init();
  };
}

#endif
