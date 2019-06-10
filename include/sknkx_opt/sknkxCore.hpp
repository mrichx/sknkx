#ifndef __SKNKXCORE_HPP__
#define __SKNKXCORE_HPP__

#include <string>
#include <vector>
#include <string.h>
#include <iostream>

// using namespace::std;

template <class T>
class matrix{
private:
  std::vector<T> data_;
  size_t nrow_=0, ncol_=0;
public:
  matrix(){ data_.resize(0,T());}
  const T& at(size_t x, size_t y) const { return data_.at(x + y * nrow_); };
  T& at(size_t x, size_t y) { return data_.at(x + y * nrow_); };
  const T& at(size_t l) const {return data_.at(l); };
  T& at(size_t l) {return data_.at(l); }
  const T& operator[](size_t l) const {return data_.at(l); };
  T& operator[](size_t l) {return data_.at(l); };
  const T* data() const {return data_.data(); };
  T* data() {return data_.data(); }
  void resize(size_t nrow, size_t ncol) {
    data_.resize(nrow*ncol);
    nrow_ = nrow;
    ncol_ = ncol;
  };
};

namespace sknkx {

  struct parameter_t {
    std::string name;
    int    index;
    double scaling;
    double value;
    double lowerBound;
    double upperBound;
    parameter_t() : index(-1), scaling(1.0), value(0.0), lowerBound(0.0), upperBound(0.0) {}
    parameter_t(std::string nm, double val, double sc, double lb, double ub){
      name  = nm;
      index = -1;
      scaling = sc;
      value = val;
      lowerBound = lb;
      upperBound = ub;
    };
    parameter_t(const parameter_t& obj){
      name       = obj.name;
      index      = -1;
      value      = obj.value;
      scaling    = obj.scaling;
      lowerBound = obj.lowerBound;
      upperBound = obj.upperBound;
    };
    // parameter_t& operator=(const parameter_t& obj)
    // {
    //   name       = obj.name;
    //   index      = obj.index;
    //   value      = obj.value;
    //   scaling    = obj.scaling;
    //   lowerBound = obj.lowerBound;
    //   upperBound = obj.upperBound;
    //   return *this;
    // }
  };

  struct output_t {
    std::string name;
    double value;
    std::vector<double> jacobian;
    //vector<vector<double>> hessian;
    matrix<double> hessian;
    std::vector<parameter_t*> params;
    output_t() : name(""), value(0.0), jacobian(0) {};
    output_t(const output_t& obj) {
      name  = obj.name;
      value = obj.value;
    };
  };

  struct constraint_t {
    std::string name;
    int index;
    double value;
    double scaling;
    double lowerBound;
    double upperBound;
    std::vector<double> jacobian;
    //vector<vector<double>> hessian;
    matrix<double> hessian;
    std::vector<int> indexesJacobian;
    std::vector<int> indexesHessian;
    std::vector<parameter_t*> params;
    std::vector<output_t*>    additionalInputs;
    constraint_t() : index(-1), value(0.0), scaling(1.0), lowerBound(0.0), upperBound(0.0) {};
    constraint_t(double sc, double lb, double ub) : index(-1), value(0.0), scaling(sc), lowerBound(lb), upperBound(ub) {};
  };

  class ModelFunc {
  public:
    virtual void evaluate();
    virtual void init();
    void addParam( parameter_t *param);
    void addDerivative(output_t *derivative);
    void addOutput(output_t *output);
    int setIndexes(int startIndex, std::vector<parameter_t*> *p);
    void resetIndexes();
    void printInfo();
  protected:
    std::vector<parameter_t*> params;
    //vector<output_t*> stateDerivatives;
    std::vector<output_t*> outputs;
    matrix<double> test;
  };

  class OdeFunc : public ModelFunc {
  public:
    OdeFunc();
    virtual void evaluate()=0;
    virtual void init()=0;
    void addParameter( parameter_t &param );
    void getStates( std::vector<parameter_t*> *states );
    void getControls( std::vector<parameter_t*> *controls );
    void getStatesDot( std::vector<output_t*> *statesDot );
    void getOutputs( std::vector<output_t*> *outputs );
  public:
    std::vector<parameter_t> states_;
    std::vector<parameter_t> controls_;
    std::vector<parameter_t> params_;
    std::vector<output_t>    statesDot_;
    std::vector<output_t>    outputs_;
  };

  class ConstraintFunc {
  public:
    virtual void evaluate();
    virtual void init();
    void returnValue( double *f );
    void returnJacobian( double *g );
    void returnHessian( double *h, const double *lambda );
    void addParam( parameter_t *param );
    void addInput( output_t *input );
    void addConstraint(constraint_t *constraint);
    int setIndexesF (int startIndex, std::vector<constraint_t*> *f);
    int setIndexesP (int startIndex, std::vector<parameter_t*> *p);
    int setIndexesJ (int startIndex, std::vector<int> *irow, std::vector<int> *jcol);
    void setIndexesH (std::vector<int> *hrow, std::vector<int> *hcol,
                      std::vector< std::vector< std::vector<int> > > &hs);
    void resetIndexes();
    void returnIndexesG( std::vector<int> &index );
    std::vector<double> getScaling();
    void printInfo();
  protected:
    std::vector<parameter_t*> params;
    std::vector<output_t*> additionalInputs;
    std::vector<constraint_t*> constraints;
  };

}

#endif
