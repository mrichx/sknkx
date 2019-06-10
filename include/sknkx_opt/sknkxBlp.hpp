#ifndef __SKNKXBLP_HPP__
#define __SKNKXBLP_HPP__

#include "sknkxCore.hpp"
#include "sknkxNlp.hpp"
#include "IpIpoptApplication.hpp"
#include "IpIpoptAlg.hpp"
#include "IpSearchDirCalculator.hpp"
#include "IpPDSearchDirCalc.hpp"
#include "IpPDSystemSolver.hpp"
#include "IpVector.hpp"
#include "IpSmartPtr.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpIpoptData.hpp"
#include "IpTNLPAdapter.hpp"
#include "IpOrigIpoptNLP.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

// using namespace Ipopt;
// using namespace std;

namespace sknkx
{

class SknkxBlp : public SknkxNlp {
public:
  // default constructor
  SknkxBlp();
  // default destructor
  virtual ~SknkxBlp();
  // add lower level problem
  //void addNlp(sknkxNlp *nlp);
  // get sensitivity of nlp wrt params
  //void getSensitivity();

protected:
  // ipopt application (one for each lower level problem)
//    vector<SmartPtr<IpoptApplication>> ipopt;
  // list of lower level problems included in bilevel problem
//    vector<sknkxNlp*> nlp;

};

class BlModelFunc : public ModelFunc {
public:
  //void addParam( parameter_t *param );
  void addNlp(Ipopt::SmartPtr<SknkxNlp> nlp, Ipopt::SmartPtr<IpoptApplication> &app);
  void init();
  // set ipopt options
  // void setIpoptNumericOptions(string option, double value);
  // void setIpoptIntegerOptions(string option, int value);
  // void setIpoptStringOptions(string option, string value);
public: // (for debugging)
  void getSensitivity();
  std::vector<std::vector<double> > dxdp;
  std::vector<std::vector<double> > dldp;
  std::vector<int> indexesConstr;
  SmartPtr<IpoptApplication> ipopt;
  SmartPtr<SknkxNlp> nlp;
};

class BlConstraintFunc : public ConstraintFunc {
public:
  //void evaluate();
  void addBlParam( parameter_t *blparam );
protected:
  std::vector<parameter_t*> blparams;
};

} /* namespace sknkx */

#endif
