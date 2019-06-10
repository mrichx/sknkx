#ifndef __SKNKXNLP_HPP__
#define __SKNKXNLP_HPP__

#include "sknkxCore.hpp"
#include "sknkxUtil.hpp"
#include "IpTNLP.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpIpoptData.hpp"
#include "IpDenseVector.hpp"
#include <array>
#include <iostream>
#include <random>
// #include "omp.h"

using namespace Ipopt;

namespace sknkx{

  class SknkxNlp : public TNLP {
  public:
    // default constructor
    SknkxNlp();
    // default deconstructor
    virtual ~SknkxNlp();

    // functions necessary for ipopt, get problem size
    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                              Index& nnz_h_lag, IndexStyleEnum& index_style);

    /** Method to return the bounds for my problem */
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                                 Index m, Number* g_l, Number* g_u);

    /** Method to return the starting point for the algorithm */
    virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                    bool init_z, Number* z_L, Number* z_U,
                                    Index m, bool init_lambda,
                                    Number* lambda);

    /** Method to return the objective value */
    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

    /** Method to return the gradient of the objective */
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

    /** Method to return the constraint residuals */
    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

    /** Method to return:
     *   1) The structure of the jacobian (if "values" is NULL)
     *   2) The values of the jacobian (if "values" is not NULL)
     */
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                            Index m, Index nele_jac, Index* iRow, Index *jCol,
                            Number* values);

    /** Method to return:
     *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
     *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
     */
    virtual bool eval_h(Index n, const Number* x, bool new_x,
                        Number obj_factor, Index m, const Number* lambda,
                        bool new_lambda, Index nele_hess, Index* iRow,
                        Index* jCol, Number* values);

    virtual bool get_scaling_parameters(Number& obj_scaling,
                                        bool& use_x_scaling, Index n,
                                        Number* x_scaling,
                                        bool& use_g_scaling, Index m,
                                        Number* g_scaling);

    virtual void finalize_solution(SolverReturn status,
                                   Index n, const Number* x,
                                   const Number* z_L, const Number* z_U,
                                   Index m, const Number* g, const Number* lambda,
                                   Number obj_value,
                                   const IpoptData* ip_data,
                                   IpoptCalculatedQuantities* ip_cq);

    virtual bool get_var_con_metadata(Index n,
                                      StringMetaDataMapType& var_string_md,
                                      IntegerMetaDataMapType& var_integer_md,
                                      NumericMetaDataMapType& var_numeric_md,
                                      Index m,
                                      StringMetaDataMapType& con_string_md,
                                      IntegerMetaDataMapType& con_integer_md,
                                      NumericMetaDataMapType& con_numeric_md);

    virtual void finalize_metadata(Index n,
                                   const StringMetaDataMapType& var_string_md,
                                   const IntegerMetaDataMapType& var_integer_md,
                                   const NumericMetaDataMapType& var_numeric_md,
                                   Index m,
                                   const StringMetaDataMapType& con_string_md,
                                   const IntegerMetaDataMapType& con_integer_md,
                                   const NumericMetaDataMapType& con_numeric_md);
    // set values of optimization parameters
    void setValues(const Number *x);

    // initialize problem (set indexes, ...)
    virtual void init();

    // evaluate cost&constraints and jacobian of problem
    virtual void evaluate();

    // add model function to problem
    void addModelFunc(ModelFunc *modelFunc);

    // add constraint function to problem
    void addConstraintFunc(ConstraintFunc *constrFunc);
    void setBlConstraintFunc(ConstraintFunc *constrFunc);

    // add cost function to problem
    void addCostFunc(ConstraintFunc *costFunc);

    // plot solution
    virtual void plot(int window);

    // get timing data
    void printTimer();

    // calculated data from ipopt
    const IpoptData *IpData_;
    IpoptCalculatedQuantities *IpCq_;

    // container for optimization parameters
    std::vector<parameter_t*> params;
    // container for constraints
    std::vector<constraint_t*> constraints;

    // solver status
    int solveStatus = -999;

  protected:
    /**
     * \brief next index for optimization parameters
     */
    int nextx;

    /**
     * \brief next index for constraints
     */
    int nextg;

    /**
     * \brief next index into linear jacobian     std::vector
     */
    int nextj;

    // containers for models and constraints
    std::vector<ModelFunc*> modelFuncs;
    std::vector<ConstraintFunc*> constrFuncs;

    /**
     * \brief container for cost functions pointers of nlp
     */
    std::vector<ConstraintFunc*> costFunc;

    /**
     * \brief pointer to bi-level constraint function,
     *        used to pass upper level parameters to nlp
     */
    ConstraintFunc *blConstrFunc=0;

    //     std::vectors for jacobian structure (row and column indexes)
    std::vector<int> irow, jcol;
    std::vector<     std::vector<int> > frow, fcol;
    std::vector<int> hrow, hcol;
    std::vector<     std::vector<     std::vector<int> > > hessianStructure;
    std::vector<     std::vector<double> > grad_f;
    // default copy and assignment operator
    SknkxNlp(const SknkxNlp&);
    SknkxNlp& operator=(const SknkxNlp&);
    // wall-timer for different functions
    SknkxTimer timer;
    double teval=0;
    double tjac =0;
    double thess=0;
    double tinit=0;
  };

}
#endif
