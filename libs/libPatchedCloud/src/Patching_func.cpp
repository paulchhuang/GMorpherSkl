/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#include "Patching.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <queue>

namespace SklRgedPatchedCloud
{

// #####################################################
// builds the patch_vertex vector from the vertex_patch
// #####################################################
void VP_To_PV(int numPatches,
			  const std::vector<int>&         vertex_patch,
			  std::vector<int>&               patch_vertex,
			  std::vector<int>&               patch_vertex_boundaries)
{
	int numVertices = vertex_patch.size();

	patch_vertex.clear();
	patch_vertex.reserve(numVertices);
	patch_vertex_boundaries.resize(numPatches+1);


	// sort PV by patch
	std::vector<std::pair<int, int> > PandV(numVertices);
	for(int vi=0;vi<numVertices; ++vi)
	{
		PandV[vi].first  = vertex_patch[vi];
		PandV[vi].second = vi;
	}
	std::sort(PandV.begin(), PandV.end());

	std::vector<std::pair<int, int> >::const_iterator pair_itr = PandV.begin();
	for(int pi = 0; pi < numPatches; ++pi)
	{
		patch_vertex_boundaries[pi] = patch_vertex.size();
		while ( (pair_itr != PandV.end() ) && ( pair_itr->first == pi ) )
		{
			patch_vertex.push_back(pair_itr->second);
			pair_itr++;
		}
	}
	patch_vertex_boundaries[numPatches] = patch_vertex.size();
	assert( int(patch_vertex.size() ) == numVertices );
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
		const std::vector<int>::const_iterator nBegin = adj.begin() + adj_boundaries[vi];
		const std::vector<int>::const_iterator nEnd   = adj.begin() + adj_boundaries[vi+1];

		oadj_boundaries[vi] = oadj.size();

		for(std::vector<int>::const_iterator n_itr = nBegin; n_itr != nEnd; ++n_itr)
		{
			int ni = *n_itr;
			assert( ni != vi);
			if( ni > vi ) oadj.push_back(ni);
			else // else we make sure that this edge has already been accounted for ( the adjacency has to be bidirectionnal )
			{
				std::vector<int>::iterator fbegin = oadj.begin() + oadj_boundaries[ni];
				std::vector<int>::iterator fend   = oadj.begin() + oadj_boundaries[ni+1];
				std::vector<int>::iterator find_itr = std::find(fbegin, fend, vi);
				assert(find_itr != fend);
			}
		}
	}

	oadj_boundaries[numVertices] = oadj.size();

}



int getConnectedComponents( const std::vector<int>& sorted_patch_adj,
                            const std::vector<int>& sorted_patch_adj_boundaries,
                            std::vector<int>&       patch_component )
{
	int numPatches = sorted_patch_adj_boundaries.size() - 1;

	patch_component.resize(numPatches);
	std::fill(patch_component.begin(), patch_component.end(), -1);

	int ccId = 0;
	int numUnassigned = numPatches;
	std::queue<int> processQueue;
	while( numUnassigned > 0)
	{
		int seedIndex = 0;
		while( patch_component[seedIndex] != -1 ) seedIndex++;

		patch_component[seedIndex] = ccId;
		processQueue.push(seedIndex);
		numUnassigned--;

		while( !processQueue.empty() )
		{
			int visitedId = processQueue.front();
			processQueue.pop();

			const std::vector<int>::const_iterator nBegin = sorted_patch_adj.begin() + sorted_patch_adj_boundaries[visitedId];
			const std::vector<int>::const_iterator nEnd = sorted_patch_adj.begin() + sorted_patch_adj_boundaries[visitedId+1];

			for( std::vector<int>::const_iterator n_itr = nBegin; n_itr != nEnd; ++n_itr)
			{
				if( patch_component[*n_itr] == -1 ) {
					patch_component[*n_itr] = ccId;
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


} // end namespace PatchedCloud
