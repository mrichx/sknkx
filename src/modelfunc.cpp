#include "sknkx_opt/sknkxCore.hpp"
#include <iostream>

using namespace std;
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
    if (this->params.at(i)->index == -1) {
      this->params.at(i)->index = startIndex++;
        // add parameter to parameter vector p
      if (p->size() <= this->params.at(i)->index) {
        p->resize( this->params.at(i)->index+1 );
        p->operator[](this->params.at(i)->index) = this->params.at(i);
        } else {
        p->operator[](this->params.at(i)->index) = this->params.at(i);
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
    this->params.at(i)->index = -1;
  }
}

/*------------------------------------------------------------------------------
 set correct size of jacobian and hessian of model outputs
 MRich, 2016-12-10
------------------------------------------------------------------------------*/
void ModelFunc::init() {
  // reset indexes
  //this->resetIndexes();
  // iterate over model outputs
  for (vector<output_t*>::size_type i=0; i!=this->outputs.size(); i++) {
    // set inputs
    //this->outputs[i]->params.resize(this->params.size());
    this->outputs.at(i)->params.clear();
    for (vector<parameter_t*>::size_type j=0; j!=this->params.size(); j++) {
      //this->outputs[i]->params[j] = this->params[j];
      this->outputs.at(i)->params.push_back(this->params.at(j));
    }
    // set size of jacobian
    // this->outputs.at(i)->jacobian.clear();
    // this->outputs.at(i)->jacobian.reserve(this->params.size());
    this->outputs.at(i)->jacobian.resize(this->params.size());
    // set size of hessian
    //this->outputs[i]->hessian.resize(this->params.size(), vector<double>(this->params.size(), 0.));
    this->outputs.at(i)->hessian.resize(this->params.size(), this->params.size());
  }
}

/*------------------------------------------------------------------------------
 virtual function for model evaluation???
 MRich, 2016-12-17
------------------------------------------------------------------------------*/
void ModelFunc::evaluate() {}

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

/*------------------------------------------------------------------------------
 class OdeFunc
 constructor
 MRich, 2016-12-30
------------------------------------------------------------------------------*/
OdeFunc::OdeFunc() {}

/*------------------------------------------------------------------------------
 class OdeFunc
 init
 MRich, 2016-12-30
------------------------------------------------------------------------------*/
void OdeFunc::init() {
  // add all states & controls to the model func parameter container
  for (vector<parameter_t>::size_type i=0; i!=this->states_.size(); i++) {
    // add state to parameters
    this->addParam(&this->states_[i]);
    // generate state derivative
    // output_t xdot;
    // xdot.name = this->states_[i].name + "_dot";
    // this->statesDot_.push_back(xdot);
    this->addOutput(&this->statesDot_[i]);
  }
  for (vector<parameter_t>::size_type i=0; i!=this->controls_.size(); i++) {
    this->addParam(&this->controls_[i]);
  }
  // add additional outputs to model func
  for (vector<output_t>::size_type i=0; i!=this->outputs_.size(); i++) {
    this->addOutput(&this->outputs_[i]);
  }
  ModelFunc::init();
}

/*------------------------------------------------------------------------------
 class OdeFunc
 add parameter
 MRich, 2016-12-30
------------------------------------------------------------------------------*/
void OdeFunc::addParameter( parameter_t &param ) {
  this->params_.push_back(param);
  this->addParam(&this->params_.back());
}

/*------------------------------------------------------------------------------
 class OdeFunc
 get pointer to states
 MRich, 2016-12-30
------------------------------------------------------------------------------*/
void OdeFunc::getStates( vector<parameter_t*> *states ) {
  states->clear();
  for (vector<parameter_t>::size_type i=0;i!=this->states_.size(); i++) {
    states->push_back( &this->states_[i] );
  }

}

/*------------------------------------------------------------------------------
 class OdeFunc
 get pointer to controls
 MRich, 2016-12-30
------------------------------------------------------------------------------*/
void OdeFunc::getControls( vector<parameter_t*> *controls ) {
  controls->clear();
  for (vector<parameter_t>::size_type i=0; i!=this->controls_.size(); i++) {
    controls->push_back( &this->controls_[i] );
  }
}

/*------------------------------------------------------------------------------
 class OdeFunc
 get pointer to outputs
 MRich, 2016-12-30
------------------------------------------------------------------------------*/
void OdeFunc::getOutputs( vector<output_t*> *outputs ){
  outputs->clear();
  for (vector<output_t>::size_type i=0; i!=this->outputs_.size(); i++) {
    outputs->push_back( &this->outputs_[i] );
  }
}

/*------------------------------------------------------------------------------
 class OdeFunc
 get pointer to state derivatives
 MRich, 2016-12-30
------------------------------------------------------------------------------*/
void OdeFunc::getStatesDot( vector<output_t*> *statesDot ) {
  statesDot->clear();
  for (vector<parameter_t>::size_type i=0; i!=this->statesDot_.size(); i++) {
    statesDot->push_back( &this->statesDot_[i] );
  }
}
