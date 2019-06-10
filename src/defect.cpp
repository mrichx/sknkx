#include "sknkx_opt/defect.hpp"
//#include <cassert>

using namespace sknkx;

/*------------------------------------------------------------------------------
 evaluate defect
 MRich, 2016-12-10
------------------------------------------------------------------------------*/
void Defect::evaluate() {
  // calculate defect value
  this->constraints[0]->value =
    +this->params[1]->value
    -this->params[0]->value
    -this->params[2]->value/2. * (this->additionalInputs[1]->value
                                 +this->additionalInputs[0]->value);

  // calculate jacobian
  this->constraints[0]->jacobian[0] = -1.;
  this->constraints[0]->jacobian[1] = +1.;
  this->constraints[0]->jacobian[2] = -.5 * (this->additionalInputs[1]->value
                                            +this->additionalInputs[0]->value);
  this->constraints[0]->jacobian[3] = -this->params[2]->value/2.;
  this->constraints[0]->jacobian[4] = -this->params[2]->value/2.;

  // calculate hessian
  this->constraints[0]->hessian.at(2,3) = -.5;
  this->constraints[0]->hessian.at(3,2) = -.5;
  this->constraints[0]->hessian.at(2,4) = -.5;
  this->constraints[0]->hessian.at(4,2) = -.5;

}

/*------------------------------------------------------------------------------
 init defect, check that params, inputs and outputs have correct size
 MRich, 2016-12-10
------------------------------------------------------------------------------*/
// void Defect::init() {
//   // check that the defect has exactly three parameters and two additional
//   // inputs and one output
//   assert(this->params.size()==3);
//   assert(this->additionalInputs.size()==2);
//   assert(this->constraints.size()==1);
//   // set correct size of jacobian and hessian
//   this->constraints[0]->jacobian.resize(this->params.size());
//   // set size of hessian
//   //this->constraints[0]->hessian.resize(this->params.size()*this->params.size());
// }
