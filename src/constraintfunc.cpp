#include "sknkx_opt/sknkxCore.hpp"

using namespace std;
using namespace sknkx;

/*------------------------------------------------------------------------------
 add parameter to constraint function
 MRich, 2016-12-10
------------------------------------------------------------------------------*/
void ConstraintFunc::addParam(parameter_t *param) {
  // add the parameter at the end of the parameter container
  this->params.push_back( param );
}

/*------------------------------------------------------------------------------
 add additional input to constraint function
 MRich, 2016-12-10
------------------------------------------------------------------------------*/
void ConstraintFunc::addInput(output_t *input) {
  // add the parameter at the end of the parameter container
  this->additionalInputs.push_back( input );
}

/*------------------------------------------------------------------------------
 add constraint to constraint function
 MRich, 2016-12-10
------------------------------------------------------------------------------*/
void ConstraintFunc::addConstraint (constraint_t *constraint) {
  // add the parameter at the end of the parameter container
  this->constraints.push_back( constraint );
}

/*------------------------------------------------------------------------------
 set indexes of states, controls and parameters and add pointers to the
 optimization parameters to a vector (-> should be the vector of optimization
 variables)
 MRich, 2016-12-10
------------------------------------------------------------------------------*/
int ConstraintFunc::setIndexesP( int startIndex, vector<parameter_t*> *p ) {
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

  // return index for next variable
  return startIndex;
}

/*------------------------------------------------------------------------------
 reset indexes of associated states, controls and parameters
 MRich, 2016-12-10
------------------------------------------------------------------------------*/
void ConstraintFunc::resetIndexes() {
  // iterate over params
  for (vector<parameter_t*>::size_type i=0; i != this->params.size(); i++) {
    this->params[i]->index = -1;
  }
}

/*------------------------------------------------------------------------------
 set indexes of associated constraints into constraint vector
 MRich, 2016-12-10
------------------------------------------------------------------------------*/
int ConstraintFunc::setIndexesF(int startIndex, vector<constraint_t*> *f) {
  // iterate over constraints and check if the constraint has already an index
  for (vector<constraint_t*>::size_type i=0; i!= this->constraints.size(); i++) {
    // check if constraint already has an index
    if (this->constraints[i]->index == -1) {
      // set index
      this->constraints[i]->index = startIndex++;
      // add constraint to constraint vector f
      if (f!=NULL) {
        if (f->size() <= this->constraints[i]->index) {
          f->resize( this->constraints[i]->index+1 );
          f->operator[](this->constraints[i]->index) = this->constraints[i];
        } else {
          f->operator[](this->constraints[i]->index) = this->constraints[i];
        }
      } else {
      }
    }
  }
  // return index for next constraint_t
  return startIndex;
}

/*------------------------------------------------------------------------------
 set indexes of jacobian for constraints
 MRich, 2016-12-10
------------------------------------------------------------------------------*/
int ConstraintFunc::setIndexesJ(int startIndex, vector<int> *irow, vector<int> *jcol) {
  int index;
  // iterate over constraints
  for (vector<constraint_t*>::size_type i=0; i!= this->constraints.size(); i++) {
    // empty jacobian index vector of current constraint
    this->constraints[i]->indexesJacobian.clear();
    // set index to zero
    index = 0;
    // iterate over parameters
    for (vector<parameter_t*>::size_type j=0; j!= this->params.size(); j++) {
      // add entry for jacobian
      if (irow!=NULL) {
        irow->push_back(this->constraints[i]->index);
        jcol->push_back(this->params[j]->index);
      }
      auto it = this->constraints[i]->indexesJacobian.begin();
      this->constraints[i]->indexesJacobian.insert(it+index++, startIndex++);
    }
    // iterate over additional inputs (outputs of model func)
    for (vector<output_t*>::size_type j=0; j!=this->additionalInputs.size(); j++) {
      // iterate over all params the output depends on
      for (vector<parameter_t*>::size_type k=0;
           k!=this->additionalInputs[j]->params.size(); k++) {
        //std::cout << "k = " << k << std::endl;
        // set flag that signalizes if the current param k is also directly
        // associated to the constraint
        int isAssociated = 0;
        // now we need to check if the current param k is also directly
        // associated to the constraint
        // for (vector<parameter_t*>::size_type l=0; l!=this->params.size(); l++) {
        //   // compare indexes of param k and param l
        //   if (this->additionalInputs[j]->params[k]->index == this->params[l]->index) {
        //     // set flag to true
        //     isAssociated = 1;
        //     // set correct index into jacobian index vector of current constraints
        //     auto it = this->constraints[i]->indexesJacobian.begin();
        //     this->constraints[i]->indexesJacobian.insert(it+index++,
        //       this->constraints[i]->indexesJacobian[l]);
        //     // leave for loop
        //     break;
        //   }
        // }
        // check if the index has not been set (isAssociated=0)
        if (isAssociated==0) {
          // add new entry into irow, jcol
          if (irow!=NULL) {
            irow->push_back(this->constraints[i]->index);
            jcol->push_back(this->additionalInputs[j]->params[k]->index);
          }
          auto it = this->constraints[i]->indexesJacobian.begin();
          this->constraints[i]->indexesJacobian.insert(it+index++, startIndex++);
        }
      }
    }
  }

  return startIndex;
}

/*------------------------------------------------------------------------------
 set indexes of hessian for constraints
 MRich, 2016-12-25
------------------------------------------------------------------------------*/
void ConstraintFunc::setIndexesH(vector<int> *hrow, vector<int> *hcol, vector<vector<vector<int>>> &hs) {
  // linear index
  int index = 0;
  int nParams = this->params.size();
  // iterate over constraints
  //#pragma omp parallel for
  for (int i=0; i<this->constraints.size(); i++) {
  //for (vector<constraint_t*>::size_type i=0; i!=this->constraints.size(); i++) {
    // iterate over all parameters
    for (vector<parameter_t*>::size_type j=0; j!=this->params.size(); j++) {
      // iterate over all parameters
      for (vector<parameter_t*>::size_type k=0; k<=j; k++) {
        // set hessian index
        hrow->push_back(this->params[j]->index);
        hcol->push_back(this->params[k]->index);
        // set hessian index
        this->constraints[i]->indexesHessian.push_back(hrow->size()-1);
      }
    }
    // iterate over additional inputs
    for (vector<output_t*>::size_type j=0; j!=this->additionalInputs.size(); j++) {
      // iterate over parameters of additional input j
      for (vector<parameter_t>::size_type k=0; k!=this->additionalInputs[j]->params.size(); k++){
        // iterate over all direct parameters of constraint func
        for (vector<parameter_t*>::size_type l=0; l!=this->params.size(); l++) {
          // get row and column index
          int rowIndex = min(this->additionalInputs[j]->params[k]->index,
                                this->params[l]->index);
          int colIndex = max(this->additionalInputs[j]->params[k]->index,
                                this->params[l]->index);
          // set hessian index
          hrow->push_back(rowIndex);
          hcol->push_back(colIndex);
          // set hessian index
          this->constraints[i]->indexesHessian.push_back(hrow->size()-1);
        }
        // iterate over all additional inputs
        //for (vector<output_t*>::size_type l=0; l!=this->additionalInputs.size(); l++) {
        for (vector<output_t*>::size_type l=0; l!=this->additionalInputs.size(); l++) {
          // iterate over all parameters of additional input
          for (vector<parameter_t*>::size_type m=0; m!=this->additionalInputs[l]->params.size(); m++) {
            // skip of index of parameter m is greater than index of parameter k
            if (this->additionalInputs[l]->params[m]->index >
                this->additionalInputs[j]->params[k]->index) {
              continue;
            }  else if (j==l) {
              // set hessian index
              hrow->push_back(this->additionalInputs[l]->params[m]->index);
              hcol->push_back(this->additionalInputs[j]->params[k]->index);
              // set hessian index
              this->constraints[i]->indexesHessian.push_back(hrow->size()-1);
            } else {
              // set hessian index
              hrow->push_back(this->additionalInputs[l]->params[m]->index);
              hcol->push_back(this->additionalInputs[j]->params[k]->index);
              // set hessian index
              this->constraints[i]->indexesHessian.push_back(hrow->size()-1);
            }
          }
        }
      }
    }
  }
  // // vector of parameters that are used in constraint function
  // vector<int> pindex;
  // int rowIndex, colIndex;
  // bool foundCombination = false;
  // // iterate over parameters and add indexes to pindex
  // for (vector<parameter_t*>::size_type i=0; i!=this->params.size(); i++) {
  //   pindex.push_back(this->params[i]->index);
  // }
  // // iterate over additional inputs and add all index of parameters the input
  // // depends on to pindex
  // for (vector<output_t*>::size_type i=0; i!=this->additionalInputs.size(); i++) {
  //   // iterate over parameters of additional input
  //   for (vector<parameter_t*>::size_type j=0; j!=this->additionalInputs[i]->params.size(); j++) {
  //     pindex.push_back(this->additionalInputs[i]->params[j]->index);
  //   }
  // }
  //
  // // iterate over constraints
  // for (vector<constraint_t*>::size_type i=0; i!=this->constraints.size(); i++) {
  //   // reset hessian index vector of current constraints
  //   this->constraints[i]->indexesHessian.clear();
  //   // iterate over all parameters of constraint functions -> pindex
  //   for (vector<int>::size_type j=0; j!=pindex.size(); j++) {
  //     for (vector<int>::size_type k=0; k<=j; k++) {
  //       // get smaller index of parameters at k and j as we only return the
  //       // lower triangular hessian to ipopt (hessian is symmetric)
  //       rowIndex = min(pindex[j],pindex[k]);
  //       colIndex = max(pindex[j],pindex[k]);
  //       // check if combination of row and column indexes is already present
  //       foundCombination = false;
  //       // for (vector<int>::size_type l=0;l!=hs[rowIndex].size();l++) {
  //       //   // check column index
  //       //   if (hs[rowIndex][l][0]==colIndex) {
  //       //     // set index of hessian
  //       //     this->constraints[i]->indexesHessian.push_back(hs[rowIndex][l][1]);
  //       //     foundCombination = true;
  //       //     break;
  //       //   }
  //       // }
  //       // for (vector<int>::size_type l=0; l!=hrow->size(); l++) {
  //       //   if ((rowIndex==hrow->data()[l]) && (colIndex==hcol->data()[l])) {
  //       //     // set index of hessian
  //       //     this->constraints[i]->indexesHessian.push_back(l);
  //       //     foundCombination = true;
  //       //     break;
  //       //   }
  //       // }
  //       // if (row,col) combination is not present so far, we need to add it
  //       if (!foundCombination) {
  //         // add indexes of row and column to hrow,hcol
  //         hrow->push_back(rowIndex);
  //         hcol->push_back(colIndex);
  //         // update hessian structure
  //         vector<int> v = {colIndex,(int)hrow->size()-1};
  //         hs[rowIndex].push_back(v);
  //         // set hessian index
  //         this->constraints[i]->indexesHessian.push_back(hrow->size()-1);
  //       }
  //     }
  //   }
  // }
}

/*------------------------------------------------------------------------------
 insert value of constraints into global constraint vector g
 MRich, 2016-12-10
------------------------------------------------------------------------------*/
void ConstraintFunc::returnValue( double *g ) {
  // iterate over constraints
  for (vector<constraint_t*>::size_type i=0; i!=this->constraints.size(); i++) {
    // add value of current constraint into f
    g[this->constraints[i]->index] = this->constraints[i]->value;
  }
}

/*------------------------------------------------------------------------------
 return indexes of constraints into g
 MRich, 2016-12-29
------------------------------------------------------------------------------*/
void ConstraintFunc::returnIndexesG( vector<int> &indexes ) {
  // reset indexes
  indexes.clear();
  // iterate over constraints
  for (vector<constraint_t*>::size_type i=0; i!=this->constraints.size(); i++) {
    indexes.push_back(this->constraints[i]->index);
  }
}

/*------------------------------------------------------------------------------
 insert jacobian of constraints into global jacobian vector g
 MRich, 2016-12-10
------------------------------------------------------------------------------*/
void ConstraintFunc::returnJacobian( double *g ) {
  // iterate over associated constraints
  for (vector<constraint_t*>::size_type i=0; i!=this->constraints.size(); i++) {
    //cout << "constraint #" << this->constraits[i]->index << endl;
    // iterate over parameters
    for (vector<parameter_t*>::size_type j=0; j!=this->params.size(); j++) {
      // add jacobian of current constraint to g
      //cout << "param #" << j << endl;
      g[this->constraints[i]->indexesJacobian[j]] =
        this->constraints[i]->jacobian[j];
    }
    // offset into indexesJacobian (first params, then additional inputs)
    int nParams = this->params.size();
    int idxJac  = this->params.size();
    // iterate over additional inputs
    for (vector<output_t*>::size_type j=0; j!=this->additionalInputs.size(); j++) {
      // iterate over parameters the input depends on and get the jacobian
      // by applying the chain rule (df/dx = df/dy * dy/dx)
      //cout << "additional input " << j << endl;
      for (vector<parameter_t*>::size_type k=0; k!=this->additionalInputs[j]->params.size(); k++){
        //cout << "g: " << this->constraints[i]->indexesJacobian[idxJac] << " =";
        g[this->constraints[i]->indexesJacobian[idxJac++]] =
          this->constraints[i]->jacobian[j+nParams] *
          this->additionalInputs[j]->jacobian[k];
        //cout << g[this->constraints[i]->indexesJacobian[j+idxJac]] << endl;
        //idxJac++;
      }
    }
  }
}

/*------------------------------------------------------------------------------
 insert hessian of constraints into global jacobian vector h
 MRich, 2016-12-25
------------------------------------------------------------------------------*/
void ConstraintFunc::returnHessian( double *h, const double *lambda ) {
  // linear index
  int nParams = this->params.size();
  // iterate over constraints
  //#pragma omp parallel for
  for (int i=0; i<this->constraints.size(); i++) {
    int index = 0;
  //for (vector<constraint_t*>::size_type i=0; i!=this->constraints.size(); i++) {
    // iterate over all parameters
    for (vector<parameter_t*>::size_type j=0; j!=this->params.size(); j++) {
      // iterate over all parameters
      for (vector<parameter_t*>::size_type k=0; k<=j; k++) {
        // set hessian
        h[this->constraints[i]->indexesHessian[index++]] =
            this->constraints[i]->hessian.at(j,k)*
            lambda[this->constraints[i]->index];
      }
    }
    // iterate over additional inputs
    for (vector<output_t*>::size_type j=0; j!=this->additionalInputs.size(); j++) {
      // iterate over parameters of additional input j
      for (vector<parameter_t>::size_type k=0; k!=this->additionalInputs[j]->params.size(); k++){
        // iterate over all direct parameters of constraint func
        for (vector<parameter_t*>::size_type l=0; l!=this->params.size(); l++) {
          h[this->constraints[i]->indexesHessian[index++]] =
            (this->constraints[i]->hessian.at(nParams+j,l)*
             this->additionalInputs[j]->jacobian[k])*
            lambda[this->constraints[i]->index];
        }
        // iterate over all additional inputs
        //for (vector<output_t*>::size_type l=0; l!=this->additionalInputs.size(); l++) {
        for (vector<output_t*>::size_type l=0; l!=this->additionalInputs.size(); l++) {
          // iterate over all parameters of additional input
          for (vector<parameter_t*>::size_type m=0; m!=this->additionalInputs[l]->params.size(); m++) {
            // skip of index of parameter m is greater than index of parameter k
            if (this->additionalInputs[l]->params[m]->index >
                this->additionalInputs[j]->params[k]->index) {
              continue;
            }  else if (j==l) {
              // use second derivative of output wrt parameters if j==l
              h[this->constraints[i]->indexesHessian[index++]] =
                (this->additionalInputs[j]->jacobian[k]*
                 this->constraints[i]->hessian.at(nParams+j,nParams+l)*
                 this->additionalInputs[l]->jacobian[m]
                +this->constraints[i]->jacobian[nParams+j]*
                 this->additionalInputs[j]->hessian.at(k,m))*
                lambda[this->constraints[i]->index];
            } else {
              //std::cout << j << " " << l << " " << k << " " << m << std::endl;
              h[this->constraints [i]->indexesHessian[index++]] =
                (this->additionalInputs[j]->jacobian[k]*
                 this->constraints[i]->hessian.at(nParams+j,nParams+l)*
                 this->additionalInputs[l]->jacobian[m])*
                lambda[this->constraints[i]->index];
            }
          }
        }
      }
    }
  }
}

/*------------------------------------------------------------------------------
 init correct sizes of jacobian and hessian of constraints
 MRich, 2016-12-10
------------------------------------------------------------------------------*/
void ConstraintFunc::init() {
    // reset indexes of constraints and optimization parameters
    //this->resetIndexes();
    // associate parameters of constraint func to each constraint
    for (vector<constraint_t*>::size_type i=0; i!=this->constraints.size(); i++) {
      // iterate over all parameters and add parameter to constraint
      for (vector<parameter_t*>::size_type j=0; j!=this->params.size(); j++) {
        this->constraints[i]->params.push_back( this->params[j] );
      }
      // iterate over additional inputs of constraint func and add additional
      // input to constraint
      for (vector<output_t*>::size_type j=0; j!=this->additionalInputs.size(); j++) {
        this->constraints[i]->additionalInputs.push_back( this->additionalInputs[j] );
      }
    }
    // get number of params and additional inputs
    int nInputs = this->params.size() + this->additionalInputs.size();
    // iterate over all constraints and set the size of the jacobian and
    // hessian
    for (vector<constraint_t*>::size_type i=0; i!=this->constraints.size(); i++) {
      // set size of jacobian
      this->constraints[i]->jacobian.resize(nInputs);
      // set size of hessian
      //this->constraints[i]->hessian.resize(nInputs, vector<double>(nInputs, 0.));
      this->constraints[i]->hessian.resize(nInputs, nInputs);
    }
}

/*------------------------------------------------------------------------------
 virtual function for constraint evaluation???
 MRich, 2016-12-17
------------------------------------------------------------------------------*/
void ConstraintFunc::evaluate() {}

/*------------------------------------------------------------------------------
 Return scaling factors for constraints
 MRich, 2016-12-30
------------------------------------------------------------------------------*/
vector<double> ConstraintFunc::getScaling() {
  vector<double> scaling;// = new vector<double>;
  for (vector<constraint_t*>::size_type i=0; i!=this->constraints.size(); i++) {
    scaling.push_back(this->constraints[i]->scaling);
  }
  return scaling;
}

/*------------------------------------------------------------------------------
 print information about constraint function
 MRich, 2016-12-10
------------------------------------------------------------------------------*/
void ConstraintFunc::printInfo() {
  std::cout << "number of parameters:  " << this->params.size() << std::endl;
  std::cout << "number of add inputs:  " << this->additionalInputs.size() << std::endl;
  std::cout << "number of constraints: " << this->constraints.size() << std::endl;
}
