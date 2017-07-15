/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef GRAPHUTILS_H_DEFINED
#define GRAPHUTILS_H_DEFINED


#include <vector>

namespace GMorpher
{


int getConnectedComponents( const std::vector<int>& sorted_adj,
                            const std::vector<int>& sorted_adj_bounds,
                            std::vector<int>&       vertex_component );


// ------------------------------------------
// removes any edge (i,j) where i < j and makes sure that if there is one we have (j,i) stored already
void orientAdjVec( const std::vector<int>&  adj,
                   const std::vector<int>&  adj_boundaries,
                   std::vector<int>&        oadj,
                   std::vector<int>&        oadj_boundaries);

// ------------------------------------------
// find the undirected graph from the oriented graph
void unOrientAdjVec( const std::vector<int>& oadj,
                     const std::vector<int>& oadj_bounds,
                     std::vector<int>&       sorted_adj,
                     std::vector<int>&       sorted_adj_bounds );


// ------------------------------------------
// find the index of the oriented edge for each unoriented edge
// edgeIds will be the same size as adj
void find_OrientedEdgeIds( const std::vector<int>& adj,
                           const std::vector<int>& adj_bounds,
                           const std::vector<int>& oadj,
                           const std::vector<int>& oadj_bounds,
                           std::vector<int>&       edgeIds );




} // end namespace GMorpher



#endif
