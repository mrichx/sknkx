#include "sknkxCore.hpp"
#include "sknkxNlp.hpp"
#include <vector>

using namespace std;

namespace sknkx {

  class SknkxPhase {
  public:
    void addState(parameter_t *state);
    void addControl(parameter_t *control);
    void addParam(parameter_t *param);
    void addOutput(output_t *output);
    void addConstraint(constraint_t *constraint);
    void addModelFunc(ModelFunc *modelFunc);
    void addConstrFunc(ConstraintFunc *constrFunc);
    void setCollocationPts(int ncoll);
    void addPhaseNlp( SknkxNlp &nlp );
  protected:
    // prototypes for states and control variable
    vector<parameter_t*>   state_proto_;
    vector<parameter_t*>   control_proto_;
    vector<output_t*>      output_proto_;
    vector<constraint_t*>  constraint_proto_;
    ModelFunc*      model_proto_;
    ConstraintFunc* constr_proto_;
    // containers to store states and controls at each collocation point
    // and all additional parameters
    vector<vector<parameter_t>> states_;
    vector<vector<parameter_t>> controls_;
    vector<vector<output_t>> outputs_;
    vector<vector<constraint_t>> constraints_;
    vector<parameter_t*> params_;
    // containers to store model and constraint funcs at each collocation point
    vector<ModelFunc> modelFuncs_;
    vector<ConstraintFunc> constrFuncs_;
    // number of collocation points
    int ncoll_ = 0;
    // flag if phase is added to nlp (can't add anything after that)
    bool isInit = false;
  };

}
