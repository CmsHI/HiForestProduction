#ifndef _FastPairwiseNearestNeighbor_
#define _FastPairwiseNearestNeighbor_

#include <utility>
#include <vector>

#include "TVectorD.h"

class Cluster;

class FastPairwiseNearestNeighbor
{
 public:
  FastPairwiseNearestNeighbor(double dMax);

  void run(const std::vector<std::pair<double,double> > & points,
                 std::vector<std::pair<TVectorD, TVectorD> > & clusters,
                 unsigned int & nOptimal,
                 std::vector<std::vector<int> > & lists,
                 unsigned int maxClusters); 

 private:
  void findClosest(std::vector<Cluster> & clusters,
                   std::vector<Cluster>::iterator c1);
  double dMax;
};

#endif
