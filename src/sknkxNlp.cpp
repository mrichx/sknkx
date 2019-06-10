#include "sknkx_opt/sknkxNlp.hpp"

using namespace std;

namespace sknkx {
  /*------------------------------------------------------------------------------
   Default constructor / destructor
   MRich, 2016-12-17
  ------------------------------------------------------------------------------*/
  SknkxNlp::SknkxNlp() {}
  SknkxNlp::~SknkxNlp() {}

  /*------------------------------------------------------------------------------
   Finalize solution
   MRich, 2016-12-17
  ------------------------------------------------------------------------------*/
  bool SknkxNlp::get_var_con_metadata(Index n,
                                      StringMetaDataMapType& var_string_md,
                                      IntegerMetaDataMapType& var_integer_md,
                                      NumericMetaDataMapType& var_numeric_md,
                                      Index m,
                                      StringMetaDataMapType& con_string_md,
                                      IntegerMetaDataMapType& con_integer_md,
                                      NumericMetaDataMapType& con_numeric_md)
  {
    // check if a bilevel constraint func exists
    if (blConstrFunc==NULL) {
      return true;
    }

    vector<Index> blconstr(m,0);
    vector<int> indexes;
    this->blConstrFunc->returnIndexesG( indexes );
    for (int i=0; i!=indexes.size(); i++) {
      blconstr[indexes[i]] = (i+1);
    }
    con_integer_md["blconstr"] = blconstr;

    return true;
  }

  void SknkxNlp::finalize_metadata(Index n,
                 const StringMetaDataMapType& var_string_md,
                 const IntegerMetaDataMapType& var_integer_md,
                 const NumericMetaDataMapType& var_numeric_md,
                 Index m,
                 const StringMetaDataMapType& con_string_md,
                 const IntegerMetaDataMapType& con_integer_md,
                 const NumericMetaDataMapType& con_numeric_md)
  {}

  /*------------------------------------------------------------------------------
   Finalize solution
   MRich, 2016-12-17
  ------------------------------------------------------------------------------*/
  void SknkxNlp::finalize_solution(SolverReturn status,
                                   Index n, const Number* x,
                                   const Number* z_L, const Number* z_U,
                                   Index m, const Number* g, const Number* lambda,
                                   Number obj_value,
                                   const IpoptData* ip_data,
                                   IpoptCalculatedQuantities* ip_cq)
  {
    // save return status
    this->solveStatus = status;
    // save data calculated by IPOPT
    this->IpData_ = ip_data;
    this->IpCq_   = ip_cq;
  }

  /*------------------------------------------------------------------------------
   Return number of optimization parameters, constraints and number of nonzeros
   in jacobian and hessian
   MRich, 2016-12-14
  ------------------------------------------------------------------------------*/
  bool SknkxNlp::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                              Index& nnz_h_lag, IndexStyleEnum& index_style) {
    // set number of variables in problem
    n = this->params.size();

    // set number of constraints in problem
    m = this->constraints.size();

    // number of nonzeros in jacobian
    nnz_jac_g = this->irow.size();

    // nummber of nonzeros in hessian
    nnz_h_lag = this->hrow.size();

    // use zero-bases indexing (c-style)
    index_style = C_STYLE;

    return true;
  }

  /*------------------------------------------------------------------------------
   Return lower and upper bounds on optimization parameters and constraints
   MRich, 2016-12-14
  ------------------------------------------------------------------------------*/
  bool SknkxNlp::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                 Index m, Number* g_l, Number* g_u) {

    // set lower and upper bounds on optimization variables
    for (vector<parameter_t*>::size_type i=0; i!=this->params.size(); i++) {
      x_l[i] = this->params[i]->lowerBound;
      x_u[i] = this->params[i]->upperBound;
    }

    // set lower and upper bounds on constraints
    for (vector<constraint_t*>::size_type i=0; i!=this->constraints.size(); i++){
      g_l[i] = this->constraints[i]->lowerBound;
      g_u[i] = this->constraints[i]->upperBound;
    }

    return true;
  }

  /*------------------------------------------------------------------------------
   Return initial guess for optimization parameters
   MRich, 2016-12-14
  ------------------------------------------------------------------------------*/
  bool SknkxNlp::get_starting_point(Index n, bool init_x, Number* x,
                                    bool init_z, Number* z_L, Number* z_U,
                                    Index m, bool init_lambda,
                                    Number* lambda){
    // set starting point of nlp
    for (vector<parameter_t>::size_type i=0; i!=this->params.size(); i++) {
      // use current values of parameters as starting point
      x[this->params[i]->index] = this->params[i]->value;
    }

    return true;
  }

  /*------------------------------------------------------------------------------
   Return values of cost function
   MRich, 2016-12-14
  ------------------------------------------------------------------------------*/
  bool SknkxNlp::eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
    // if a new x-vector is provided calculate values/jacobian of cost and
    // constraints

    if (new_x) {
      // set values of optimization parameters
      this->setValues(x);
      // evaluate problem
      this->evaluate();
    }

    // reset objective value (might be set from previous call)
    obj_value = 0.0;

    // set cost function value
    double tmp_value = 0.0;
    for (vector<ConstraintFunc*>::size_type i=0; i!=this->costFunc.size();++i)
    {
      this->costFunc[i]->returnValue(&tmp_value);
      obj_value += tmp_value;
    }

    return true;
  }

  /*------------------------------------------------------------------------------
   Return gradient of cost function
   MRich, 2016-12-14
  ------------------------------------------------------------------------------*/
  bool SknkxNlp::eval_grad_f(Index n, const Number* x, bool new_x, Number* gradf) {
    //
    if (new_x) {
      // set values of optimization parameters
      this->setValues(x);
      // evaluate problem
      this->evaluate();
    }

    // reset grad_f as it might not be all zeros
    memset( gradf, 0.0, n*sizeof(Number) );

    // set gradient of cost function
    for (vector<ConstraintFunc*>::size_type i=0; i!=this->costFunc.size(); ++i)
    {
      this->costFunc[i]->returnJacobian(this->grad_f[i].data());
      for (vector<int>::size_type j=0; j!=this->fcol[i].size(); ++j) {
        gradf[this->fcol[i][j]] += this->grad_f[i][j];
      }
    }

    return true;
  }

  /*------------------------------------------------------------------------------
   Return values of constraints
   MRich, 2016-12-14
  ------------------------------------------------------------------------------*/
  bool SknkxNlp::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
    // evaluate problem if new x
    if (new_x) {
      // set values of optimization parameters
      this->setValues(x);
      // evaluate problem
      this->evaluate();
    }

    // set constraint values
    for (vector<ConstraintFunc*>::size_type i=0; i!=this->constrFuncs.size(); i++){
      this->constrFuncs[i]->returnValue(g);
    }

    return  true;
  }

  /*------------------------------------------------------------------------------
   Return values or structure of jacobian of constraints
   "values" == NULL : return structure
   "values" != NULL : return values
   MRich, 2016-12-14
  ------------------------------------------------------------------------------*/
  bool SknkxNlp::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index *iRow, Index *jCol,
                           Number* values)
  {
    // evaluate problem if new x
    if (new_x) {
      // set values of optimization parameters
      this->setValues(x);
      // evaluate problem
      this->evaluate();
    }

    // set structure of constraint jacobian if "values" == NULL
    if (values==NULL) {
      for (vector<int>::size_type i=0; i!=this->irow.size(); i++) {
        iRow[i] = this->irow[i];
        jCol[i] = this->jcol[i];
      }
    } else {
      // start timer
      this->timer.tic();
      // reset values as we use += to fill them
      memset(values,0.0,nele_jac*sizeof(Number));
      // set values of constraint jacobian if "values" != NULL
      #ifdef _OPENMP
      #pragma omp parallel for shared(values)
      #endif
      for (int i=0; i<this->constrFuncs.size(); i++) {
        this->constrFuncs[i]->returnJacobian(values);
      }
      // get time for execution
      this->tjac += this->timer.toc();
    }

    return true;
  }

  /*------------------------------------------------------------------------------
   Evaluate hessian of lagrangian
   MRich, 2016-12-17
  ------------------------------------------------------------------------------*/
  bool SknkxNlp::eval_h(Index n, const Number* x, bool new_x,
                        Number obj_factor, Index m, const Number* lambda,
                        bool new_lambda, Index nele_hess, Index* iRow,
                        Index* jCol, Number* values) {
    // if a new x-vector is provided, recalculate cost & constraints
    if (new_x) {
      // set values of optimization parameters
      this->setValues(x);
      // evaluate problem
      this->evaluate();
    }

    // get value of hessian if value is not NULL, else return structure of hessian
    if (values == NULL) {
      // return structure of hessian matrix
      for (vector<int>::size_type i=0; i<this->hrow.size(); i++) {
        iRow[i] = this->hrow[i];
        jCol[i] = this->hcol[i];
      }
    } else {
      // start timer
      this->timer.tic();
      // reset hessian values
      memset(values,0.0,nele_hess*sizeof(Number));
      // return hessian of cost function
      for (int i=0; i<this->costFunc.size();++i)
      {
        this->costFunc[i]->returnHessian(values,&obj_factor);
      }
      // set values of hessian
      #ifdef _OPENMP
      #pragma omp parallel for shared(values,lambda)
      #endif
      //for (vector<ConstraintFunc*>::size_type i=0; i!=this->constrFuncs.size();i++){
      for (int i=0; i<this->constrFuncs.size();i++){
        //std::cout << i << std::endl;
        //#pragma omp task
        this->constrFuncs[i]->returnHessian(values,lambda);
      }
      // get time for execution
      this->thess += this->timer.toc();
    }

    return true;
  }

/*------------------------------------------------------------------------------
 Get scaling parameters for ipopt
 MRich, 2016-12-30
------------------------------------------------------------------------------*/
  bool SknkxNlp::get_scaling_parameters(Number& obj_scaling,
                                    bool& use_x_scaling, Index n,
                                    Number* x_scaling,
                                    bool& use_g_scaling, Index m,
                                    Number* g_scaling)
  {
    // acticate scaling for parameters and constraints
    use_x_scaling = true;
    use_g_scaling = true;
    // TODO: THIS IS STUPID HERE!!!
    // set scaling factor for objective function
    vector<double> costScaling = this->costFunc[0]->getScaling();
    for (vector<double>::size_type i=0; i!=costScaling.size(); i++) {
      obj_scaling = costScaling[i];
    }

    // iterate over parameters and return scaling parameters
    for (vector<parameter_t*>::size_type i=0; i!=this->params.size(); i++) {
      x_scaling[i] = this->params[i]->scaling;
    }
    // set scaling factors for constraints
    for (vector<constraint_t*>::size_type i=0; i!=this->constraints.size(); i++) {
      g_scaling[i] = this->constraints[i]->scaling;
    }
    return true;
  }


/*------------------------------------------------------------------------------
 Set values of optimization parameters
 MRich, 2016-12-14
------------------------------------------------------------------------------*/
  void SknkxNlp::setValues(const Number* x) {
    // iterate over optimization parameters and set value
    for (vector<parameter_t*>::size_type i=0; i!=this->params.size(); i++) {
      // assign value to parameter
      this->params[i]->value = x[this->params[i]->index];
    }
  }

  /*------------------------------------------------------------------------------
   Initialize problem
    1.) set indexes of optimization parameters
    2.) set indexes of constraints
    3.) set indexes of constraint jacobian
   MRich, 2016-12-14
  ------------------------------------------------------------------------------*/
  void SknkxNlp::init() {
    this->timer.tic();
    // init indexes for optimization parameters and constraints
    this->nextx=0;
    this->nextg=0;
    this->nextj=0;

    // iterate over model funcs
    for (vector<ModelFunc*>::size_type i=0; i!=this->modelFuncs.size(); i++) {
      // init model func
      this->modelFuncs[i]->init();
      // set indexes of optimization parameters
      this->nextx = this->modelFuncs[i]->setIndexes( this->nextx,
                                                    &(this->params) );
    }

    // set indexes for cost functions
    for (vector<ConstraintFunc*>::size_type i=0; i!=this->costFunc.size(); ++i)
    {
      // add space to store gradient indexes of current cost function
      this->frow.push_back(vector<int>());
      this->fcol.push_back(vector<int>());
      // init cost function
      this->costFunc[i]->init();
      this->nextx = this->costFunc[i]->setIndexesP(this->nextx, &(this->params));
      this->costFunc[i]->setIndexesF(0, NULL);
      this->costFunc[i]->setIndexesJ(0, &(this->frow[i]), &(this->fcol[i]));
      this->hessianStructure.resize(this->params.size(),vector<vector<int>>());
      this->costFunc[i]->setIndexesH(&(this->hrow), &(this->hcol), (this->hessianStructure));
      this->grad_f.push_back( vector<double>(this->frow[i].size(),0.0) );
    }

    // iterate over constraint functions
    for (vector<ConstraintFunc*>::size_type i=0; i!=this->constrFuncs.size(); i++) {
      // init constraint func
      this->constrFuncs[i]->init();
      // set indexes of optimization parameters
      this->nextx = this->constrFuncs[i]->setIndexesP( this->nextx,
                                                       &(this->params) );
      // set indexes of constraints
      this->nextg = this->constrFuncs[i]->setIndexesF( this->nextg,
                                                       &(this->constraints) );
      // set indexes of jacobian
      this->nextj = this->constrFuncs[i]->setIndexesJ( this->nextj,
                                                       &(this->irow),
                                                       &(this->jcol) );
      // set indexes of hessian
      this->hessianStructure.resize(this->params.size(),vector<vector<int>>());
      this->constrFuncs[i]->setIndexesH(&(this->hrow), &(this->hcol),
                                         (this->hessianStructure));
    }
    this->tinit = this->timer.toc();
  }

  /*------------------------------------------------------------------------------
   Evaluate model funcs and constraint funcs of problem
   MRich, 2016-12-17
  ------------------------------------------------------------------------------*/
  void SknkxNlp::evaluate() {
    this->timer.tic();
    #ifdef _OPENMP
    //cout << "parallel calculation" << endl;
    #pragma omp parallel
    #endif
    {
      // evaluate all model functions
      #ifdef _OPENMP
      #pragma omp for
      #endif
      for (int i=0; i<this->modelFuncs.size(); i++) {
        this->modelFuncs[i]->evaluate();
      }
      // evaluate all constraint functions
      #ifdef _OPENMP
      #pragma omp for
      #endif
      for(int i=0; i<this->constrFuncs.size(); i++) {
        this->constrFuncs[i]->evaluate();
      }
    }
    // evaluate cost function
    for (size_t i=0; i<this->costFunc.size(); ++i)
    {
      this->costFunc[i]->evaluate();
    }
    // get time for execution
    this->teval += this->timer.toc();
  }
  /*------------------------------------------------------------------------------
   Add model function to problem
   MRich, 2016-12-17
  ------------------------------------------------------------------------------*/
  void SknkxNlp::addModelFunc(ModelFunc* modelFunc) {
    // add model func to the model function container
    this->modelFuncs.push_back(modelFunc);
  }

  /*------------------------------------------------------------------------------
   Add constraint function to problem
   MRich, 2016-12-17
  ------------------------------------------------------------------------------*/
  void SknkxNlp::addConstraintFunc(ConstraintFunc* constrFunc) {
    // add constraint func to constraint function container
    this->constrFuncs.push_back(constrFunc);
  }

  /*------------------------------------------------------------------------------
   Set bilevel constraint function to problem
   MRich, 2016-12-29
   ------------------------------------------------------------------------------*/
   void SknkxNlp::setBlConstraintFunc(ConstraintFunc* constrFunc) {
    // add constraint func to constraint function container
    this->constrFuncs.push_back(constrFunc);
    // and set it as bilevel constraint function of problem
    this->blConstrFunc = constrFunc;
  }

  /*------------------------------------------------------------------------------
   Add cost function to problem
   MRich, 2016-12-17
  ------------------------------------------------------------------------------*/
  void SknkxNlp::addCostFunc(ConstraintFunc* costFunc) {
    // set cost function
    this->costFunc.push_back(costFunc);
  }

  /*------------------------------------------------------------------------------
   plot function
   MRich, 2017-02-11
  ------------------------------------------------------------------------------*/
  void SknkxNlp::plot(int window) {

  }

  /*------------------------------------------------------------------------------
   Print collected timing
   MRich, 2016-12-17
  ------------------------------------------------------------------------------*/
  void SknkxNlp::printTimer() {
    std::cout << "t init: " << this->tinit << std::endl;
    std::cout << "t eval: " << this->teval << std::endl;
    std::cout << "t jac:  " << this->tjac  << std::endl;
    std::cout << "t hess: " << this->thess << std::endl;
  }

  /*------------------------------------------------------------------------------
   Check derivatives of cost & constraint functions
   MRich, 2016-12-20
  ------------------------------------------------------------------------------*/
  // void SknkxNlp::checkDerivatives() {
  //   // init local variables
  //   double *x, val, *g, *constrJacAnl, *constrJac;
  //
  //   // random number generator
  //   std::default_random_engine generator;
  //   std::uniform_real_distribution<double> distribution(0.0,1.0);
  //   // set up x vector
  //   for (vector<parameter_t*>::size_type i=0; i!=this->params.size(); i++) {
  //     // set parameter value to random number (uniform distribution) that is
  //     // between lower and upper bound
  //     this->params[i]->value = this->params[i]->lowerBound +
  //                               distribution(generator) *
  //                               (this->params[i]->upperBound -
  //                                this->params[i]->lowerBound);
  //   }
  //
  //   // evaluate cost & constraint functions for the unperturbed parameters
  //   this->eval_g(int n, const int *x, bool new_x, int m, int *g)
  //
  //   // iterate over
  //
  // }
}
