#include "sknkxPhase.hpp"

namespace sknkx {
/*------------------------------------------------------------------------------
 Default constructor / destructor
 MRich, 2016-12-17
------------------------------------------------------------------------------*/
  // SknkxPhase::SknkxPhase()  {}
  // SknkxPhase::~SknkxPhase() {}

/*------------------------------------------------------------------------------
 add state to phase
 MRich, 2016-12-29
------------------------------------------------------------------------------*/
  void SknkxPhase::addState(parameter_t *state){
    // check that initialization of phase is not over yet
    if (!this->isInit) {
      this->state_proto_.push_back(state);
    } else {
      cout << "phase already attached to nlp, can't add anythin anymore!" << endl;
    }
  }

/*------------------------------------------------------------------------------
 add control to phase
 MRich, 2016-12-29
------------------------------------------------------------------------------*/
  void SknkxPhase::addControl(parameter_t *control){
    // check that initialization of phase is not over yet
    if (!this->isInit) {
      this->control_proto_.push_back(control);
    } else {
      cout << "phase already attached to nlp, can't add anythin anymore!" << endl;
    }
  }

/*------------------------------------------------------------------------------
 add parameter to phase
 MRich, 2016-12-29
  ------------------------------------------------------------------------------*/
  void SknkxPhase::addParam(parameter_t *param){
    // check that initialization of phase is not over yet
    if (!this->isInit) {
      this->params_.push_back(param);
    } else {
      cout << "phase already attached to nlp, can't add anythin anymore!" << endl;
    }
  }

/*------------------------------------------------------------------------------
 add model func to phase
 MRich, 2016-12-29
------------------------------------------------------------------------------*/
  void SknkxPhase::addModelFunc(ModelFunc *modelFun){
    // check that initialization of phase is not over yet
    if (!this->isInit) {
      // this->model_proto_.push_back(modelFun);
    } else {
      cout << "phase already attached to nlp, can't add anythin anymore!" << endl;
    }
  }

/*------------------------------------------------------------------------------
 add constraint func to phase
 MRich, 2016-12-29
------------------------------------------------------------------------------*/
  void SknkxPhase::addConstrFunc(ConstraintFunc *constrFun) {
    // check that initialization of phase is not over yet
    if (!this->isInit) {
      // this->constr_proto_->push_back(constrFun);
    } else {
      cout << "phase already attached to nlp, can't add anythin anymore!" << endl;
    }
  }

/*------------------------------------------------------------------------------
 set number of collocation points
 MRich, 2016-12-29
------------------------------------------------------------------------------*/
  void SknkxPhase::setCollocationPts(int ncoll) {
    this->ncoll_ = ncoll;
  }

/*------------------------------------------------------------------------------
 add phase to nlp
 MRich, 2016-12-29
------------------------------------------------------------------------------*/
  void SknkxPhase::addPhaseNlp(SknkxNlp &nlp) {
    // generate containers for states, controls, outputs, and constraints
    this->states_.resize(this->ncoll_, vector<parameter_t>(this->state_proto_.size()));
    this->controls_.resize(this->ncoll_,vector<parameter_t>(this->control_proto_.size()));
    this->outputs_.resize(this->ncoll_,vector<output_t>(this->output_proto_.size()));
    this->constraints_.resize(this->ncoll_,vector<constraint_t>(this->constraint_proto_.size()));
    // generate container for model func and constraint func of phase
    this->modelFuncs_.resize(this->ncoll_);
    this->constrFuncs_.resize(this->ncoll_);
    // set values
    for (int i=0;i!=this->ncoll_;i++) {
      for (vector<parameter_t*>::size_type j=0; j!=this->state_proto_.size(); j++) {
        this->states_[i][j] = parameter_t( *this->state_proto_[j] );
        this->states_[i][j].name += to_string(i);
        // this->states_[i][j].name       = this->state_proto_[j]->name + to_string(i);
        // this->states_[i][j].value      = this->state_proto_[j]->value;
        // this->states_[i][j].lowerBound = this->state_proto_[j]->lowerBound;
        // this->states_[i][j].upperBound = this->state_proto_[j]->upperBound;
      }
      for (vector<parameter_t*>::size_type j=0; j!=this->control_proto_.size(); j++) {
        this->controls_[i][j] = parameter_t( *this->control_proto_[j] );
        this->controls_[i][j].name += to_string(i);
        // this->controls_[i][j].name       = this->control_proto_[j]->name + to_string(i);
        // this->controls_[i][j].value      = this->control_proto_[j]->value;
        // this->controls_[i][j].lowerBound = this->control_proto_[j]->lowerBound;
        // this->controls_[i][j].upperBound = this->control_proto_[j]->upperBound;
      }
      for (vector<output_t*>::size_type j=0; j!=this->output_proto_.size(); j++) {
        this->outputs_[i][j] = output_t( *this->output_proto_[j] );
        this->outputs_[i][j].name += to_string(i);
      }
      for (vector<constraint_t*>::size_type j=0; j!=this->constraint_proto_.size(); j++) {
        // this->constraints_[i][j] = constraint_t( *this->constraint_proto_[j] );
        // this->constraints_[i][j].name += to_string(i);
      }
      // this->modelFuncs_[i] =
    }

    //
  }
}
