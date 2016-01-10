#include "prodlim.h"

#include <iterator>
#include <algorithm>
#include <vector>
#include <utility>
#include <iostream>
#include <math.h>

using namespace std;
using namespace Rcpp;


SEXP R_fastkm(SEXP X, SEXP D, SEXP V, SEXP left, SEXP Q)
{
  List result;
  
  try {
    NumericVector xX(X);
    NumericVector xV(V);
    IntegerVector xD(D);
    
    bool ll = as<bool>(left);
    
    vector<double> sX(xX.begin(), xX.end());
    vector<double> sV(xV.begin(), xV.end());
    vector<int> sD(xD.begin(), xD.end());
    
    ProductLimit pl(sX, sD, sV, false);
    
    if(Q == R_NilValue) {
      result = List::create(Named("time") = pl.getTime(),
			    Named("surv") = pl.getSurv(),
			    Named("variance") = pl.getVarHazard(),
			    Named("n.atrisk") = pl.getNAtRisk());
    } else {
      NumericVector xQ(Q);
      vector<double> sQ(xQ.begin(), xQ.end());
      ProdlimEval x = pl.eval(sQ, ll);
      result = List::create(Named("time") = sQ,
			    Named("surv") = x.surv,
			    Named("variance") = x.var,
			    Named("n.atrisk") = x.natrisk);
    }
  } catch(std::exception &ex) {	
    forward_exception_to_r(ex);
  } 

  return result;
}

// compare floating point numbers x and y allowing for numerical tolerance eps
// int _equal(double x, double y, double eps)
// {
//   double xy = fabs(x - y);
//   double xn = fabs(x);

//   if(xn > eps) xy = xy / xn;  // relative difference
//   if(xy > eps) return 0;
//   else return 1;
// }


// X, status, V required to be sorted!
ProductLimit::ProductLimit(const vector<double>& X, const vector<int>& status, const vector<double>& V, bool reverse)
{
  int n = X.size();

  vector<pair<double, int> > tmp(n);
  for(int i = 0; i < n; ++i) {
    tmp[i] = pair<double, int>(X[i], status[i]);
  }
  // count unique elements
  int m = 1;
  for(int i = 0; i < n - 1; ++i) {
    if(tmp[i].first != tmp[i+1].first) ++m;
  }

  // pre-allocate result vectors
  mTime.resize(m); 
  mSurv.resize(m);
  mHazard.resize(m);
  mVarhazard.resize(m);
  mNatrisk.resize(m);

  // number-at-risk is needed to handle left-truncation
  vector<int> atrisk(n);
  for(int j = 0; j < n; ++j) {
    int cur_n = 0;
    for(int k = 0; k < n; ++k) {
      if((V[k] <= X[j]) && (X[j] <= X[k])) cur_n++;
    }
    atrisk[j] = cur_n;
  }

  // calculate product limit estimator
  surv(tmp, atrisk, reverse, 0);
}

// survival probability and variance (of -log(surv)) (Greenwood) estimate at each time in Y
ProdlimEval ProductLimit::eval(const vector<double>& Y, bool left_limit)
{
  int n = Y.size();
  vector<double> res(n);
  vector<double> var(n);
  vector<int> nar(n);
  vector<double>::iterator low;
  int m = mTime.size() - 1;

  if (left_limit) {  // left-continuous piecewise constant interpolation
    for (int i = 0; i < n; ++i) {
      if (Y[i] <= mTime[0]) {
	res[i] = 1;
	var[i] = 0;
	nar[i] = mNatrisk[0];
      } else if(Y[i] > mTime[m]) {
	res[i] = mSurv[m];
	var[i] = mVarhazard[m];
	nar[i] = mNatrisk[m];
      } else {
	low = lower_bound(mTime.begin(), mTime.end(), Y[i]);
	int a = int(low - mTime.begin());
	if (*low == Y[i]) {
	  res[i] = mSurv[a-1];
	  var[i] = mVarhazard[a-1];
	  nar[i] = mNatrisk[a-1];
	} else {
	  res[i] = mSurv[a];
	  var[i] = mVarhazard[a];
	  nar[i] = mNatrisk[a];
	}
      }
    }
  } else {
    for (int i = 0; i < n; ++i) {  // right-continuous piecewise constant interpolation
      if (Y[i] < mTime[0]) {
	res[i] = 1;
	var[i] = 0;
	nar[i] = mNatrisk[0];
      } else if(Y[i] >= mTime[m]) {
	res[i] = mSurv[m];
	var[i] = mVarhazard[m];
	nar[i] = mNatrisk[m];
      } else {
	low = lower_bound(mTime.begin(), mTime.end(), Y[i]);
	int b = int(low - mTime.begin());
	if(Y[i] < mTime[b]) b--;
	res[i] = mSurv[b];
	var[i] = mVarhazard[b];
	nar[i] = mNatrisk[b];
      }
    }
  }

  return ProdlimEval(res, var, nar, left_limit);
}

// survival probability and variance (of -log(surv)) (Greenwood) estimate at each time in Y (only works when Y is sorted in ascending order!)
pair<vector<double>, vector<double> > ProductLimit::eval_sorted(const vector<double>& Y, bool left_limit)
{
  int n = Y.size();
  vector<double> res(n);
  vector<double> var(n);
  int m = mTime.size() - 1;
  vector<double>::iterator prev = mTime.begin();

  if (left_limit) {  // left-continuous piecewise constant interpolation
    for (int i = 0; i < n; ++i) {
      if (Y[i] <= mTime[0]) {
	res[i] = 1;
	var[i] = 0;
      } else if(Y[i] > mTime[m]) {
	res[i] = mSurv[m];
	var[i] = mVarhazard[m];
      } else {
	prev = lower_bound(prev, mTime.end(), Y[i]); // start search from previous position since Y is sorted
	int a = int(prev - mTime.begin());
	if (*prev == Y[i]) {
	  res[i] = mSurv[a-1];
	  var[i] = mVarhazard[a-1];
	} else {
	  res[i] = mSurv[a];
	  var[i] = mVarhazard[a];
	}
      }
    }
  } else {
    
    for (int i = 0; i < n; ++i) {  // right-continuous piecewise constant interpolation
      if (Y[i] < mTime[0]) {
	res[i] = 1;
	var[i] = 0;
      } else if(Y[i] >= mTime[m]) {
	res[i] = mSurv[m];
	var[i] = mVarhazard[m];
      } else {
	prev = lower_bound(prev, mTime.end(), Y[i]); // start search from previous position since Y is sorted
	int b = int(prev - mTime.begin());
	if(Y[i] < mTime[b]) b--;
	res[i] = mSurv[b];
	var[i] = mVarhazard[b];
      }
    }
  }
  return make_pair(res, var);
}

// adapted from prodlim_surv from prodlim package
void ProductLimit::surv(const vector<pair<double, int> >& data, const vector<int>& atrisk_all, bool reverse, double trunc)
{
  //int atrisk = data.size(); // this does not work for left-truncated data!
  int s = 0;
  int event = data[0].second;
  int loss = 1 - data[0].second;
  int n = data.size();
  double surv_temp = 1;
  double hazard_temp = 0;
  double varhazard_temp = 0;

  double y0 = 0;
  double status = 0;
  
  for (int i = 0; i < n; i++) {
    y0 = data[i].first;
    int atrisk = atrisk_all[i];
    
    if(i < n-1) { // look at next element -> tie handling
      status = data[i+1].second;
      
      if (y0 == data[i+1].first) { // handle ties
	event += status;
	loss  += (1 - status);
	continue;
      }    
    }
    
    mTime[s] = y0;
    
    if (reverse) {
      if (loss > 0 && atrisk > trunc) {
	hazard_temp = loss / (double) (atrisk - event);
	surv_temp *= (1 - hazard_temp);
	varhazard_temp += (double) loss / ((double) (atrisk - event) * (double) (atrisk - event - loss));
      }
    } else {	
      if (event > 0 && atrisk > trunc) {
	hazard_temp = event / (double) atrisk;
	surv_temp *= (1 - hazard_temp);
	varhazard_temp += (double) event / ((double) atrisk * (double) (atrisk - event));
      }
    }
    
    mHazard[s] = hazard_temp;
    mVarhazard[s] = varhazard_temp;
    mSurv[s] = surv_temp;
    mNatrisk[s] = atrisk;
    
    //atrisk -= (event + loss);
    if(i < n-1) {
      s++;
      event = status;
      loss = (1 - status);
    }
  }
}
