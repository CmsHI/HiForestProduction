#include "../interface/OptimalTree.h"

#include <cmath>
#include <iostream>

#define sqr(x) ((x) * (x))

using namespace std;

/*****************************************************************************/
double OptimalTree::run
  (int K,
   const vector<pair<double, double> > & points,
   TVectorD mu, TVectorD P,
   vector<vector<int> > lists,
   VertexCollection & vertices)
{
  // Clear vertices
  vertices.clear();

  // Fill
  for(int k = 0; k < K; k++)
  {
    Vertex vertex;

    double sig2 = 0;
    for(vector<int>::iterator il = lists[k].begin();
                              il!= lists[k].end(); il++)
      sig2 += 1 / points[*il].second;

    vertex.first = pair<double,double>(mu(k), 1/sig2);

    for(vector<int>::iterator il = lists[k].begin();
                              il!= lists[k].end(); il++)
        vertex.second.push_back(*il);

    vertices.push_back(vertex);
  }

  return 0.;
}

