#include "sknkx_opt/sknkxCore.hpp"
#include <iostream>

using namespace sknkx;

/*------------------------------------------------------------------------------
 add parameter to model function
 MRich, 2016-12-10
------------------------------------------------------------------------------*/
void ModelFunc::addParam(parameter_t *param) {
  // add the parameter at the end of the parameter container
  this->params.push_back( param );
}

/*------------------------------------------------------------------------------
 add output to model function
 MRich, 2016-12-10
------------------------------------------------------------------------------*/
void ModelFunc::addOutput(output_t *output) {
  // add the output to the end of the outputs
  this->outputs.push_back(output);
}

/*------------------------------------------------------------------------------
 set indexes of states, controls and parameters and add pointers to the
 optimization parameters to a vector (-> should be the vector of optimization
 variables)
 MRich, 2016-12-10
------------------------------------------------------------------------------*/
int ModelFunc::setIndexes( int startIndex, vector<parameter_t*> *p ) {

  // iterate over parameters
  for (vector<parameter_t*>::size_type i=0; i != this->params.size(); i++) {
    // set state index if it has none
    if (this->params[i]->index == -1) {
        this->params[i]->index = startIndex++;
        // add parameter to parameter vector p
        if (p->size() <= this->params[i]->index) {
          p->resize( this->params[i]->index+1 );
          p->operator[](this->params[i]->index) = this->params[i];
        } else {
          p->operator[](this->params[i]->index) = this->params[i];
        }
    }
  }
  // return index for next parameter
  return startIndex;
}

/*------------------------------------------------------------------------------
 reset indexes of associated states, controls and parameters
 MRich, 2016-12-10
------------------------------------------------------------------------------*/
void ModelFunc::resetIndexes() {
  // iterate over params
  for (vector<parameter_t*>::size_type i=0; i != this->params.size(); i++) {
    this->params[i]->index = -1;
  }
}

/*------------------------------------------------------------------------------
 set correct size of jacobian and hessian of model outputs
 MRich, 2016-12-10
------------------------------------------------------------------------------*/
void ModelFunc::init() {
  // iterate over model outputs
  for (vector<output_t*>::size_type i=0; i!=this->outputs.size(); i++) {
    // set size of jacobian
    this->outputs[i]->jacobian.resize(this->params.size());
    // set size of hessian
    this->outputs[i]->hessian.resize(this->params.size()*this->params.size());
  }
}

/*------------------------------------------------------------------------------
 show model info
 MRich, 2016-12-10
------------------------------------------------------------------------------*/
void ModelFunc::printInfo() {
    // print number of states, controls and parameters
    std::cout << "number of params:   " << this->params.size() << std::endl;
    std::cout << "param indexes: ";
    for (int i=0; i!=this->params.size(); i++) {
        std::cout << this->params[i]->index << " ";
    }
    std::cout << std::endl;
}
