#ifndef _OptimalTree_
#define _OptimalTree_

#include <utility>
#include <vector>

#include "TVectorD.h"
#include "TMatrixD.h"

typedef std::pair<std::pair<double,double>, std::vector<int> > Vertex;
typedef std::vector<Vertex> VertexCollection;

class OptimalTree
{
 public:
  OptimalTree() {}

  double run(int K, const std::vector<std::pair<double, double> > & points,
             TVectorD mu, TVectorD P,
             std::vector<std::vector<int> > lists,
             VertexCollection & vertices);

 private:
};

#endif
