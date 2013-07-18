#include "../interface/GaussianMixture.h"

#include <iostream>
#include <cmath>

#define sqr(x) ((x) * (x))

using namespace std;

/*****************************************************************************/
double GaussianMixture::logProb(const pair<double,double> & p,
                                double mu, double P)
{
  return - 0.5 * sqr(p.first - mu) / p.second + log(P);
}

/*****************************************************************************/
void GaussianMixture::estimateAverages
  (const vector<pair<double,double> > & points,
   const TMatrixD & p,
   TVectorD & mu, TVectorD & var, TVectorD & P)
{
  int K = mu.GetNrows();

  // Recalculate averages
  for(int k = 0; k < K; k++)
  {
    double sum[3] = {0, 0, 0};

    for(unsigned int n = 0; n < points.size(); n++)
    {
      sum[2] += p(n,k) * points[n].first / points[n].second;
      sum[1] += p(n,k)                   / points[n].second;
      sum[0] += p(n,k);
    }

    mu(k) = sum[2]/sum[1];
   var(k) =      1/sum[1];
     P(k) = sum[0]/points.size();
  }
}

/*****************************************************************************/
double GaussianMixture::poisson(double mu, double n)
{
  return pow(mu,n)*exp(-mu - lgamma(n+1));
}

/*****************************************************************************/
double GaussianMixture::estimateChiSquare
  (const vector<pair<double,double> > & points,
   const TVectorD & mu, const TVectorD & P,
   TMatrixD & p, bool estimateResponsibility)
{
  int K = mu.GetNrows();

  double chi2 = 0;
  for(unsigned int n = 0; n < points.size(); n++)
  {
    vector<double> logP(K);

    bool isFirst = true;
    double maxLogP = 0.;
    for(int k = 0; k < K; k++)
    if(P(k) > 0)
    {
      logP[k] = logProb(points[n],mu(k),P(k));

      if(logP[k] > maxLogP || isFirst)
      { maxLogP = logP[k]; isFirst = false; }
    }

    double sum = 0.;
    for(int k = 0; k < K; k++)
    if(P(k) > 0)
      sum += exp(logP[k] - maxLogP);

    if(estimateResponsibility)
      for(int k = 0; k < K; k++)
        if(P(k) > 0) p(n,k) = exp(logP[k] - maxLogP - log(sum));
                else p(n,k) = 0;
    else
      chi2 += -2 * (maxLogP + log(sum));
  }

  return chi2;
}

/*****************************************************************************/
double GaussianMixture::getChiSquare
  (const vector<pair<double,double> > & points,
   const TVectorD & mu, const TVectorD & P)
{
  TMatrixD p;

  return estimateChiSquare(points,mu,P, p,false);
}

/*****************************************************************************/
void GaussianMixture::estimateResponsibility
  (const vector<pair<double,double> > & points,
   const TVectorD & mu, const TVectorD & P,
   TMatrixD & p)
{
  estimateChiSquare(points,mu,P, p,true);
}

/*****************************************************************************/
double GaussianMixture::run
  (int K,
   const vector<pair<double, double> > & points,
   TVectorD mu, TVectorD P,
   VertexCollection & vertices)
{
  // Clear vertices
  vertices.clear();

  int N = points.size();

  TMatrixD p(N,K);
  TVectorD var(K);

  double chi2 = 0.;
  double old_chi2;

  int iter = 0;

  do
  {
    estimateResponsibility(points,mu,P,      p);
    estimateAverages      (points,p, mu,var, P);

    old_chi2 = chi2;
    chi2 = getChiSquare(points,mu,P);
  }
  while(fabs(chi2 - old_chi2) > 1e-3 && ++iter < 100);

  // Fill
  for(int k = 0; k < K; k++)
  {
    Vertex vertex;

    vertex.first = pair<double,double>(mu(k), var(k));

    for(int i = 0; i < N; i++)
    {
      double pmax = 0.;
      int kmax = -1;

      for(int j = 0; j < K; j++)
        if(p(i,j) > pmax)
        { pmax = p(i,j); kmax = j; }

      if(k == kmax)
        vertex.second.push_back(i);
    }

    vertices.push_back(vertex);
  }

  return chi2;
}

