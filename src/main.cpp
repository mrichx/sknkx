//#include "../include/SkunkOpt.hpp"
#include "../include/defect.hpp"
#include "../include/sknkxNlp.hpp"
#include "gnuplot-iostream.h"
#include <iostream>
#include "IpIpoptApplication.hpp"
#include "math.h"

//using namespace::std;
using namespace::SkunkWorkx;
//using namespace::std;

class BrachCost : public ConstraintFunc {
public:
  void evaluate() {
    // minimize time (i.e. time step)
    this->constraints[0]->value = this->params[0]->value * (this->ncoll-1.);
    // jacobian
    this->constraints[0]->jacobian[0] = 1. * (this->ncoll-1.);
  }
  double ncoll;
};

class BrachOde : public ModelFunc {
    void evaluate() {
      double g = 9.81;
      // x-position ode
      this->outputs[0]->value = sqrt(2*g)*
                                this->params[1]->value *
                                cos(this->params[2]->value);
      this->outputs[1]->value = sqrt(g/2)*
                                sin(this->params[2]->value);

      // jacobian
      this->outputs[0]->jacobian[1] = sqrt(2*g)*cos(this->params[2]->value);
      this->outputs[0]->jacobian[2] =-sqrt(2*g)*
                                      this->params[1]->value *
                                      sin(this->params[2]->value);
      this->outputs[1]->jacobian[2] = sqrt(g/2)*cos(this->params[2]->value);

      // hessian
      this->outputs[0]->hessian[1][2] =-sqrt(2*g)*sin(this->params[2]->value);
      this->outputs[0]->hessian[2][1] =-sqrt(2*g)*sin(this->params[2]->value);
      this->outputs[0]->hessian[2][2] =-sqrt(2*g)*this->params[1]->value*
                                        cos(this->params[2]->value);
      this->outputs[1]->hessian[2][2] =-sqrt(g/2)*sin(this->params[2]->value);

    }

};

int main() {

  // number of collocation points
  int ncoll = 1001;
  Gnuplot gp;

  // container for position, speed and acceleration
  parameter_t *xpos, *ypos, *vel, *theta;
  xpos  = new parameter_t[ncoll];
  ypos  = new parameter_t[ncoll];
  vel   = new parameter_t[ncoll];
  theta = new parameter_t[ncoll];
  // container for state derivatives
  output_t *xdot, *ydot, *vdot;
  xdot = new output_t[ncoll];
  ydot = new output_t[ncoll];
  vdot = new output_t[ncoll];

  // model functions
  BrachOde *model;
  model = new BrachOde[ncoll];
  // constraint functions
  constraint_t *xeqout, *yeqout, *veqout;
  xeqout = new constraint_t[ncoll-1];
  yeqout = new constraint_t[ncoll-1];
  veqout = new constraint_t[ncoll-1];
  Defect *xeq, *yeq, *veq;
  xeq = new Defect[ncoll-1];
  yeq = new Defect[ncoll-1];
  veq = new Defect[ncoll-1];

  // nlp
  SmartPtr<SknkxNlp> brachistochrone = new SknkxNlp();
  // ipopt
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
  app->RethrowNonIpoptException(true);
  app->Options()->SetNumericValue("tol", 1e-12);
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("hessian_approximation", "exact"); // limited-memory
  //app->Options()->SetStringValue("hessian_approximation", "limited-memory"); //
  //app->Options()->SetStringValue("output_file", "ipopt.out");
  app->Options()->SetStringValue("linear_solver", "ma57");
  app->Options()->SetStringValue("nlp_scaling_method", "none");
  app->Options()->SetStringValue("fixed_variable_treatment","make_constraint");
  //app->Options()->SetStringValue("derivative_test", "second-order");
  app->Options()->SetIntegerValue("max_iter",100);
  app->Options()->SetIntegerValue("print_level",5);

  // create time parameter
  parameter_t timeStep;
  timeStep.lowerBound = 1E-6;
  timeStep.upperBound = 1E-0;
  timeStep.value      = 0.5 / ((double)ncoll-1.);

  // add cost function to problem
  BrachCost costFunc;
  constraint_t costValue;
  costFunc.addParam(&timeStep);
  costFunc.addConstraint(&costValue);
  costFunc.ncoll = (double) ncoll;
  brachistochrone->addCostFunc(&costFunc);

  // create model functions
  for (int i=0;i<ncoll;i++) {
    // create new xpos state
    // set name
    xpos[i].name = "xpos" + to_string(i);
    // set lower and upper bound
    xpos[i].lowerBound = 0.;
    xpos[i].upperBound = 1.;
    xpos[i].value      = 1. * (double) i / (ncoll-1.);

    // create new ypos state
    ypos[i].name = "ypos" + to_string(i);
    ypos[i].lowerBound = 0.;
    ypos[i].upperBound =+1.;
    ypos[i].value      = 1. - 1. * (double) i / (ncoll-1.);

    // create new theta control
    theta[i].name = "acc" + to_string(i);
    theta[i].lowerBound = -3.1416/2;
    theta[i].upperBound = +3.1416/2;
    theta[i].value      = 0.1 + 0.1 * ((double) i / (ncoll-1));

    // create outputs for position and velocity derivative
    xdot[i].name = "xdot_" + to_string(i);
    ydot[i].name = "ydot_" + to_string(i);

    // create new model function and add states, controls and outputs
    model[i].addParam( &xpos[i] );
    model[i].addParam( &ypos[i] );
    model[i].addParam( &theta[i] );
    model[i].addOutput( &xdot[i] );
    model[i].addOutput( &ydot[i] );

    // add model function to nlp
    brachistochrone->addModelFunc( &model[i] );
  }

  // set initial boundary conditions
  xpos[0].lowerBound = 0.;
  xpos[0].upperBound = 0.;
  ypos[0].lowerBound = 0.;
  ypos[0].upperBound = 1.;

  // set final boundary conditions
  xpos[ncoll-1].lowerBound = 1.0;
  xpos[ncoll-1].upperBound = 1.0;
  ypos[ncoll-1].lowerBound = 0.0;
  ypos[ncoll-1].upperBound = 0.0;

  // create state integration defects
  for (int i=0; i<ncoll-1; i++) {
    // create new defect
    // add output, states and state derivatives to defect
    xeq[i].addParam( &xpos[i] );
    xeq[i].addParam( &xpos[i+1] );
    xeq[i].addParam( &timeStep );
    xeq[i].addInput( &xdot[i] );
    xeq[i].addInput( &xdot[i+1] );
    xeq[i].addConstraint( &xeqout[i] );
    xeqout[i].name = "xDefect" + to_string(i);
    // add defect to nlp
    brachistochrone->addConstraintFunc( &xeq[i] );

    yeq[i].addParam( &ypos[i] );
    yeq[i].addParam( &ypos[i+1] );
    yeq[i].addParam( &timeStep );
    yeq[i].addInput( &ydot[i] );
    yeq[i].addInput( &ydot[i+1] );
    yeq[i].addConstraint( &yeqout[i] );
    yeqout[i].name = "yDefect" + to_string(i);
    // add defect to nlp
    brachistochrone->addConstraintFunc( &yeq[i] );

  }

  cout << "init problem" << endl;
  // init problem (set indexes, etc.)
  brachistochrone->init();

  //
  cout << "timeStep: " << timeStep.value << endl;
  // solve problem
  ApplicationReturnStatus status;
  status = app->Initialize();
  status = app->OptimizeTNLP(brachistochrone);


  // plot solution
  vector<double> time;
  vector<double> xstate, ystate, tctrl;
  for (int i=0;i<ncoll; i++) {
    time.push_back( double(i)*timeStep.value );
    xstate.push_back( xpos[i].value ); //
    ystate.push_back( ypos[i].value );
    tctrl.push_back( theta[i].value );
  }
  gp << "plot '-' with lines title 'xpos',"
     <<       "'-' with lines title 'ypos',"
     <<       "'-' with lines title '{theta}' \n";
  gp.send1d(boost::make_tuple(time, xstate));
  gp.send1d(boost::make_tuple(time, ystate));
  gp.send1d(boost::make_tuple(time, tctrl));

  gp << "plot '-' with lines title 'trajectory'\n";
  gp.send1d(boost::make_tuple(ystate,xstate));

  // delete allocated arrays
  delete[] xpos;
  delete[] ypos;
  delete[] theta;
  delete[] model;
  delete[] xeqout;
  delete[] yeqout;
  delete[] xeq;
  delete[] yeq;

  return 0;
}
