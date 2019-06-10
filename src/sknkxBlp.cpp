#include "sknkx_opt/sknkxBlp.hpp"

namespace sknkx {

  /*------------------------------------------------------------------------------
   Class BlModelFunc
   add nlp
   MRich, 2016-12-28
  ------------------------------------------------------------------------------*/
  void BlModelFunc::addNlp(SmartPtr<SknkxNlp> nlp, SmartPtr<IpoptApplication> &app) {
    // add nlp to model func
    this->nlp = nlp;
    // clone app
    this->ipopt = app->clone();
  }

  /*------------------------------------------------------------------------------
   Class BlModelFunc::addParam
   add parameter to upper level model function and pass index of constraint
   used within the nlp to enforce the value of the parameter
   MRich, 2016-12-28
  ------------------------------------------------------------------------------*/
  // void BlModelFunc::addParam(parameter_t *param) {
  //   // add parameter just as it is done in ModelFunc
  //   ModelFunc::addParam(param);
  //   // save index of constraint
  //   //this->indexesConstr.push_back(index);
  // }

  /*------------------------------------------------------------------------------
   Class BlModelFunc::init
   init model func
   MRich, 2016-12-29
  ------------------------------------------------------------------------------*/
  void BlModelFunc::init() {
    // set size of dxdp and dldp
    // get number of parameters in nlp
    size_t npNlp = nlp->params.size();
    size_t ncNlp = nlp->constraints.size();
    // resize
    dxdp.resize(params.size(), std::vector<double>(npNlp, 0.));
    dldp.resize(params.size(), std::vector<double>(ncNlp, 0.));
    // call standard init function
    ModelFunc::init();
  }

  /*------------------------------------------------------------------------------
   Class BlModelFunc::getSensitivity
   compute sensitivity of optimal nlp solution wrt parameters
   MRich, 2016-12-28
  ------------------------------------------------------------------------------*/
  void BlModelFunc::getSensitivity() {

    // check if ipopt solved the problem successfully
    int nIter = 0;
    // init status to Internal_Error, just to make sure it is initialized
    ApplicationReturnStatus status=Internal_Error;
    while (!(status == Solve_Succeeded ||
             status == Solved_To_Acceptable_Level)
           and nIter < 3)
    {
      std::cout << "Solving lower level problem! Iteration " << nIter++ << std::endl;
      // solve nlp
      if (this->nlp->solveStatus == -999) {
        // solve nlp (this should maybe go to a separate function)
        status = this->ipopt->OptimizeTNLP(this->nlp);
      } else {
        status = this->ipopt->ReOptimizeTNLP(this->nlp);
      }
      // if problem is infeasible, write constraint values to file
      if (status == Infeasible_Problem_Detected)
      {
        // open file for writing
        std::ofstream constr_file;
        constr_file.open("constraints.info", std::ios::out);
        for (auto& constraint : this->nlp->constraints)
        {
          constr_file << std::setw(15) << std::setprecision(10) << std::left << constraint->name << ": ";
          constr_file << std::setw(15) << std::setprecision(10) << std::left << constraint->value;
          constr_file << std::endl;
        }
      }
    }

    // we should now ensure that the nlp was solved succesfully before we
    // calculate the sensitivity matrix... but who cares...

    // get tnlp adapter for resorting to original problem formulation
    TNLPAdapter* tnlp_adapter = NULL;
    OrigIpoptNLP* orignlp;
    orignlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(this->nlp->IpCq_->GetIpoptNLP()));
    tnlp_adapter = dynamic_cast<TNLPAdapter*>(GetRawPtr(orignlp->nlp()));

    // get linear solver
    SmartPtr<IpoptAlgorithm> alg = this->ipopt->AlgorithmObject();
    SmartPtr<PDSearchDirCalculator> pd_search;
    pd_search = dynamic_cast<PDSearchDirCalculator*>(GetRawPtr(alg->SearchDirCalc()));
    SmartPtr<PDSystemSolver> pd_solver = pd_search->PDSolver();

    // init right hand side
    SmartPtr<IteratesVector> rhs = this->nlp->IpData_->curr()->MakeNewContainer();
    //SmartPtr<IteratesVector> rhs = this->nlp->IpData_->trial()->MakeNewIteratesVector();

    // SmartPtr<const DenseVectorSpace> x_owner_space_ =
    //   dynamic_cast<const DenseVectorSpace*>(GetRawPtr(this->nlp->IpData_->curr()->x()->OwnerSpace())) ;
    // const std::vector<Index> idx_ipopt =
    //     x_owner_space_->GetIntegerMetaData("sens_state_1");

    SmartPtr<const DenseVectorSpace> c_space =
        dynamic_cast<const DenseVectorSpace*>(GetRawPtr(this->nlp->IpData_->curr()->y_c()->OwnerSpace())) ;
    std::vector<Index> idx_c = c_space->GetIntegerMetaData("blconstr");
    SmartPtr<DenseVector> rhs_c = new DenseVector(GetRawPtr(ConstPtr(c_space)));
    Number* rhs_c_val = rhs_c->Values();
    //
    // iterate over the parameters of the model func and compute the
    // sensitivity of the nlp wrt to the parameters
    for (size_t i=0; i!=this->params.size(); i++) {
      // reset values of rhs c
      memset(rhs_c_val,0.0,rhs_c->Dim()*sizeof(Number));
      // set values of right hand side
      for (int j=0; j!=idx_c.size(); j++) {
        if (idx_c[j]==(i+1)) {
          rhs_c_val[j] = -1.0;
        }
      }
      SmartPtr<Vector> rhs_x = this->nlp->IpCq_->curr_grad_lag_x()->MakeNewCopy();
      SmartPtr<Vector> rhs_s = this->nlp->IpCq_->curr_grad_lag_s()->MakeNewCopy();
      SmartPtr<Vector> rhs_yd= this->nlp->IpCq_->curr_d_minus_s()->MakeNewCopy();
      SmartPtr<Vector> rhs_zl= this->nlp->IpCq_->curr_compl_x_L()->MakeNewCopy();
      SmartPtr<Vector> rhs_zu= this->nlp->IpCq_->curr_compl_x_U()->MakeNewCopy();
      SmartPtr<Vector> rhs_vl= this->nlp->IpCq_->curr_compl_s_L()->MakeNewCopy();
      SmartPtr<Vector> rhs_vu= this->nlp->IpCq_->curr_compl_s_U()->MakeNewCopy();

      rhs_x->Set(0.0);
      rhs_s->Set(0.0);
      rhs_yd->Set(0.0);
      rhs_zl->Set(0.0);
      rhs_zu->Set(0.0);
      rhs_vl->Set(0.0);
      rhs_vu->Set(0.0);

      // set kkt residuals
      rhs->Set_x_NonConst(*rhs_x);
      rhs->Set_s_NonConst(*rhs_s);
      //rhs->Set_y_c_NonConst(*this->nlp->IpCq_->curr_c()->MakeNewCopy());
      rhs->Set_y_c_NonConst(*rhs_c->MakeNewCopy());
      rhs->Set_y_d_NonConst(*rhs_yd);
      rhs->Set_z_L_NonConst(*rhs_zl);
      rhs->Set_z_U_NonConst(*rhs_zu);
      rhs->Set_v_L_NonConst(*rhs_vl);
      rhs->Set_v_U_NonConst(*rhs_vu);
      //cout << "sum of rhs: " << rhs->Sum() << endl;
      // get space to store sensitivities
      SmartPtr<IteratesVector> delta =
         this->nlp->IpData_->curr()->MakeNewIteratesVector(true);

      // compute sensitivity
      bool improve_solution = false;
      bool allow_inexact = false;
      bool retval;
      retval = pd_solver->Solve(-1.0, 0.0, *rhs, *delta, allow_inexact,
                                 improve_solution);

      // get sensitivities of optimal solution
      tnlp_adapter->ResortX(*delta->x(), this->dxdp[i].data());
      tnlp_adapter->ResortG(*delta->y_c(), *delta->y_d(), this->dldp[i].data());

      // unscale sensitivities
      // cout << "unscaling sensitivities ... ";
      for (size_t j=0; j!=this->dxdp[i].size(); ++j)
        this->dxdp[i][j] /= this->nlp->params[j]->scaling;
      // cout << " done!" << endl;
    }

  }

/*------------------------------------------------------------------------------
 Class BlConstraintFunc
 addBlParam
 MRich, 2016-12-28
------------------------------------------------------------------------------*/
void BlConstraintFunc::addBlParam( parameter_t *blparam ) {
  // add to blparams
  this->blparams.push_back(blparam);
}

}
