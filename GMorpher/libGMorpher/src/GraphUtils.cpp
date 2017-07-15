/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#include <GraphUtils.h>
#include <queue>
#include <cassert>
#include <algorithm>
#include <iostream>

namespace GMorpher
{

int getConnectedComponents( const std::vector<int>& sorted_adj,
                            const std::vector<int>& sorted_adj_bounds,
                            std::vector<int>&       vertex_component )
{
	int numVertices = sorted_adj_bounds.size() - 1;

	vertex_component.resize(numVertices);
	std::fill(vertex_component.begin(), vertex_component.end(), -1);

	int ccId = 0;
	int numUnassigned = numVertices;
	std::queue<int> processQueue;
	while( numUnassigned > 0)
	{
		int seedIndex = 0;
		while( vertex_component[seedIndex] != -1 ) seedIndex++;

		vertex_component[seedIndex] = ccId;
		processQueue.push(seedIndex);
		numUnassigned--;

		while( !processQueue.empty() )
		{
			int visitedId = processQueue.front();
			processQueue.pop();

			const std::vector<int>::const_iterator nBegin = sorted_adj.begin() 
			                                     + sorted_adj_bounds[visitedId];
			const std::vector<int>::const_iterator   nEnd = sorted_adj.begin() 
			                                   + sorted_adj_bounds[visitedId+1];

			for( std::vector<int>::const_iterator n_itr = nBegin; n_itr != nEnd;
			                                                            ++n_itr)
			{
				if( vertex_component[*n_itr] == -1 ) {
					vertex_component[*n_itr] = ccId;
					processQueue.push(*n_itr);
					numUnassigned--;
				}
				// TODO (or not) assert on the fact that it is either -1 or ccId
			}
		}

		ccId++;
	}

	return ccId;
}










// #####################################################
// orients the edges of the graph
// #####################################################
void orientAdjVec( const std::vector<int>&  adj,
                   const std::vector<int>&  adj_boundaries,
                   std::vector<int>&        oadj,
                   std::vector<int>&        oadj_boundaries)
{
	int numVertices = adj_boundaries.size() -1;
	assert( int(adj_boundaries.size()) == numVertices +1);

	oadj.clear();
	oadj.reserve(adj.size());
	oadj_boundaries.resize(numVertices+1);
	
	for(int vi=0; vi<numVertices; ++vi)
	{
		const std::vector<int>::const_iterator nBegin = adj.begin() 
		                                                   + adj_boundaries[vi];
		const std::vector<int>::const_iterator nEnd   = adj.begin() 
		                                                 + adj_boundaries[vi+1];		
		oadj_boundaries[vi] = oadj.size();
		//oadj_boundaries.push_back(oadj.size());
		for(std::vector<int>::const_iterator n_itr = nBegin; n_itr != nEnd; 
		                                                                ++n_itr)
		{
			int ni = *n_itr;
			assert( ni != vi);
			if( ni > vi ) oadj.push_back(ni);
			else // else we make sure that this edge has already been accounted 
			     // for ( the adjacency has to be bidirectionnal )
			{
				std::vector<int>::iterator fbegin = oadj.begin() 
				                                          + oadj_boundaries[ni];
				std::vector<int>::iterator fend   = oadj.begin() 
				                                        + oadj_boundaries[ni+1];
				std::vector<int>::iterator find_itr = std::find(fbegin, fend, vi);
				assert(find_itr != fend);
			}
		}
	}

	oadj_boundaries[numVertices] = oadj.size();
	//oadj_boundaries.push_back(oadj.size());

}







void unOrientAdjVec( const std::vector<int>& oadj,
                     const std::vector<int>& oadj_bounds,
                     std::vector<int>&       sorted_adj,
                     std::vector<int>&       sorted_adj_bounds )
{
	int numNodes = oadj_bounds.size() - 1;
	int numEdges = oadj.size() * 2;
	// ---------
	// allocate
	sorted_adj_bounds = std::vector<int>( numNodes + 1, 0);
	sorted_adj.resize( numEdges );

	// ----------
	// create the bounds vector
		// count the number of guys
	for(int ni = 0; ni <numNodes; ++ni ) {
		const std::vector<int>::const_iterator nBegin = oadj.begin() 
		                                                      + oadj_bounds[ni];
		const std::vector<int>::const_iterator   nEnd = oadj.begin() 
		                                                    + oadj_bounds[ni+1];
		sorted_adj_bounds[ni+1] += nEnd - nBegin;
		for(std::vector<int>::const_iterator nj_itr = nBegin; nj_itr != nEnd; 
		                                                            ++nj_itr ) {
			sorted_adj_bounds[*nj_itr+1]++;
		}
	}
		// accumulate
	for(int ni = 0; ni <numNodes; ++ni ) {
		sorted_adj_bounds[ni+1] += sorted_adj_bounds[ni];
	}


	std::vector<int> offset( sorted_adj_bounds );
	// ----------
	// fill the edges
	for(int ni = 0; ni <numNodes; ++ni ) {
		const std::vector<int>::const_iterator nBegin = oadj.begin() 
		                                                      + oadj_bounds[ni];
		const std::vector<int>::const_iterator   nEnd = oadj.begin() 
		                                                    + oadj_bounds[ni+1];
		for(std::vector<int>::const_iterator nj_itr = nBegin; nj_itr != nEnd;
		                                                            ++nj_itr ) {
			int nj = *nj_itr;
			sorted_adj[ offset[ni]++ ] = nj;
			sorted_adj[ offset[nj]++ ] = ni;
		}
	}

	// ----------
	// sort the edges
	for(int ni = 0; ni <numNodes; ++ni ) {
		std::vector<int>::iterator nBegin = sorted_adj.begin() 
		                                                + sorted_adj_bounds[ni];
		std::vector<int>::iterator   nEnd = sorted_adj.begin() 
		                                              + sorted_adj_bounds[ni+1];
		std::sort( nBegin, nEnd );
	}
}






void find_OrientedEdgeIds( const std::vector<int>& adj,
                           const std::vector<int>& adj_bounds,
                           const std::vector<int>& oadj,
                           const std::vector<int>& oadj_bounds,
                           std::vector<int>&       edgeIds )
{
	int numPatches = adj_bounds.size() -1;
	assert( int( oadj_bounds.size() - 1) == numPatches );
	edgeIds.resize( adj.size() );

	for(int pi=0;pi<numPatches;++pi)
	{
		std::vector<int>::iterator              e_itr = 
		                                       edgeIds.begin() + adj_bounds[pi];
		std::vector<int>::const_iterator        n_itr = adj.begin() 
		                                                       + adj_bounds[pi];
		const std::vector<int>::const_iterator  n_end = adj.begin() 
		                                                     + adj_bounds[pi+1];
		while( n_itr != n_end ) 
		{
			int pj = *n_itr;
			if( pj > pi ) { // j > i : we look for the edge [i,j]
				const std::vector<int>::const_iterator f_itr = std::find( 
				                                 oadj.begin() + oadj_bounds[pi], 
				                         oadj.begin() + oadj_bounds[pi+1], pj );
				assert ( f_itr !=  oadj.begin() + oadj_bounds[pi+1] );
				*e_itr = f_itr - oadj.begin();
			}
			else { // i > j : we look for the edge [j,i]
				const std::vector<int>::const_iterator f_itr = std::find( 
				                                 oadj.begin() + oadj_bounds[pj],
				                         oadj.begin() + oadj_bounds[pj+1], pi );
				assert ( f_itr !=  oadj.begin() + oadj_bounds[pj+1] );
				*e_itr = f_itr - oadj.begin();
			}
			e_itr++;
			n_itr++;
		}
	}
}


} // end namespace GMorpher
