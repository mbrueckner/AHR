#ifndef _PRODLIM_H
#define _PRODLIM_H

#include <utility>
#include <vector>
#include <RcppArmadillo.h>

RcppExport SEXP R_fastkm(SEXP X, SEXP D, SEXP V, SEXP left, SEXP Q);

class ProdlimEval {
public:
  std::vector<double> surv;
  std::vector<double> var;
  std::vector<int> natrisk;
  bool left_limit;

  ProdlimEval(std::vector<double> s, std::vector<double> v, std::vector<int> n, bool l) : surv(s), var(v), natrisk(n), left_limit(l) {};
};


class ProductLimit {
public:
  ProductLimit(const std::vector<double>& X, const std::vector<int>& status, const std::vector<double>& V, bool reverse);

  ProdlimEval eval(const std::vector<double>& Y, bool left_limit);
  std::pair<std::vector<double>, std::vector<double> > eval_sorted(const std::vector<double>& Y, bool left_limit);

  std::vector<double> getSurv() { return mSurv; }
  std::vector<double> getTime() { return mTime; }
  std::vector<double> getHazard() { return mHazard; }
  std::vector<double> getVarHazard() { return mVarhazard; }
  std::vector<int> getNAtRisk() { return mNatrisk; }

private:
  void surv(const std::vector<std::pair<double, int> >& data, const std::vector<int>& atrisk_all, bool reverse, double trunc);

  static bool cmp(std::pair<double, int> a, std::pair<double, int> b) { return (a.first < b.first); }

  std::vector<double> mTime;
  std::vector<double> mSurv;
  std::vector<double> mHazard;
  std::vector<double> mVarhazard;
  std::vector<int> mNatrisk;
};


#endif
