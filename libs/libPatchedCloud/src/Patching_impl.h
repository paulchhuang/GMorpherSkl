/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef PATCHING_IMPL_H_DEFINED
#define PATCHING_IMPL_H_DEFINED

#include "Patching_func.h"
#include <queue>
#include <algorithm>
#include <limits>
#include <iostream>
#include <cassert>

namespace SklRgedPatchedCloud
{

	// ##################################
	// ##################################
	// SEEDPATCHES
	// ##################################
	// ##################################
	template <class IndexedMeshAdjacency_Wrapper>
	void SeedPatches_redux( const typename IndexedMeshAdjacency_Wrapper::Data& M, const int maxPatchSize,
	                        std::vector<int> &vertex_patch,
	                        std::vector<int> &patch_centers_indices,
	                        std::vector<int> &geoDistance)
	{
		typedef typename IndexedMeshAdjacency_Wrapper::IndexIterator IndexIterator;

		int numVertices         = IndexedMeshAdjacency_Wrapper::numVertices(M);
		int unassignedVertices  = numVertices;

		// resize the output vectors
		patch_centers_indices.clear();
		vertex_patch.resize(numVertices);
		std::fill(vertex_patch.begin(), vertex_patch.end(), -1);

		// this will hold the best distance data
		geoDistance.resize(numVertices);
		std::fill(geoDistance.begin(), geoDistance.end(), std::numeric_limits<int>::max());

		// the vector which will hold the vertices to visit
		std::queue<std::pair<int, int> >    verticesToVisit;
		std::vector<bool>                   vertexVisited(numVertices);
		std::vector<int>                    borderCount(numVertices, 0);



		/* unassignedVertices, borderCount, patch_centers_indices is for the process of whole mesh patching 
		   vertexVisited, verticesToVisit is for the process of single patch growing
		 */

		while(unassignedVertices)
		{
			// the vector must be reseted for this patch
			std::fill(vertexVisited.begin(), vertexVisited.end(), false);
			// the queue of vertices to process must be empty
			assert( verticesToVisit.empty() );


			// 1 - Find a new patch Center
			std::vector<int>::iterator maxBorder_itr = std::max_element( borderCount.begin(), borderCount.end() );
			int new_Center_id = std::distance( borderCount.begin(), maxBorder_itr);
			if( *maxBorder_itr == 0 ){ // if we have no border point left.. can happen if multiple disconnected components
				new_Center_id = 0;	// the start of the first patch, 1907
				while( vertex_patch[new_Center_id] != -1 ) new_Center_id++;	// if the current componet is done with patching, find the first vertex of the next component and start with it
			}
			//~ std::cout<<new_Center_id<<" "<<*maxBorder_itr<<" "<<unassignedVertices<<std::endl;

			// 2 - Assign the center to the new patch
			patch_centers_indices.push_back(new_Center_id);
			int patchId = patch_centers_indices.size() - 1;


			// 3 - Front propagation
			verticesToVisit.push(std::make_pair(new_Center_id, 0));
			vertexVisited[new_Center_id] = true;
			borderCount[new_Center_id] = 0;

			// now starts the process of single patch growing
			while( !verticesToVisit.empty() )
			{
				// pop a vertice from the que
				int visited_index = verticesToVisit.front().first;
				int visited_distance = verticesToVisit.front().second;
				verticesToVisit.pop();

				if(visited_distance <= geoDistance[visited_index] ) //if this new patch is better
				{
					// if this vertex was previously unassigned we can decrement the count of unassigned patches
					if(vertex_patch[visited_index] == -1) unassignedVertices--;

					// lets write ourselves as the best patch
					vertex_patch[visited_index] = patchId;
					geoDistance[visited_index] = visited_distance;

					// if we are still close enough, add our neighbours to be processed
					if(visited_distance < maxPatchSize)
					{// the following two lines extract all neighors of the vertex
						IndexIterator nbegin 	= IndexedMeshAdjacency_Wrapper::neighboursCWBegin(M, visited_index);
						IndexIterator nend		= IndexedMeshAdjacency_Wrapper::neighboursCWEnd(M, visited_index);
						for(IndexIterator ni_itr = nbegin; ni_itr != nend; ++ni_itr)
						{
							if(!vertexVisited[*ni_itr]) //if that vertex is not yet visited in this growing pass, add it into "verticesToVisit"
							{ 
								verticesToVisit.push(std::make_pair(*ni_itr, visited_distance + 1));
								vertexVisited[*ni_itr] = true; //just to make sure we only add it once to the queue
								borderCount[*ni_itr] = 0;
							}
						}
					}
					// if not.. increment the border count of our neighbours    
					else // Paul: if we already reach the border, add the neighbors as borders (increase borderCount)
					{
						IndexIterator nbegin 	= IndexedMeshAdjacency_Wrapper::neighboursCWBegin(M, visited_index);
						IndexIterator nend		= IndexedMeshAdjacency_Wrapper::neighboursCWEnd(M, visited_index);
						for(IndexIterator ni_itr = nbegin; ni_itr != nend; ++ni_itr)
						{
							if((!vertexVisited[*ni_itr]) && (vertex_patch[*ni_itr] == -1))
							{
								vertexVisited[*ni_itr] = true; //just to make sure we only add it once to the queue
								borderCount[*ni_itr] ++;
							}
						}
					}
				}
			}
		}
		std::cout<<"GeoPatcher returned "<<patch_centers_indices.size()<<" patches"<<std::endl;
	}





	// ##################################
	// ##################################
	// SEEDPATCHES
	// ##################################
	// ##################################
	template <class IndexedMeshAdjacency_Wrapper>
	void SeedPatches( const typename IndexedMeshAdjacency_Wrapper::Data& M, const int maxPatchSize,
	                  std::vector<int> &vertex_patch,
	                  std::vector<int> &patch_centers_indices,
	                  std::vector<int> &geoDistance,
	                  const std::vector<int>&  priorityCenters ) // = std::vector<int>())
	{
		typedef typename IndexedMeshAdjacency_Wrapper::IndexIterator IndexIterator;



		int numVertices         = IndexedMeshAdjacency_Wrapper::numVertices(M);
		int unassignedVertices  = numVertices;

		// resize the output vectors
		patch_centers_indices.clear();
		vertex_patch.resize(numVertices);
		std::fill(vertex_patch.begin(), vertex_patch.end(), -1);

		// this will hold the best distance data
		geoDistance.resize(numVertices);
		std::fill(geoDistance.begin(), geoDistance.end(), std::numeric_limits<int>::max());

		// the vector which will hold the vertices to visit
		std::queue<std::pair<int, int> > 	verticesToVisit;
		std::vector<bool>					vertexVisited(numVertices);

		std::vector<int>::const_iterator nextPriority = priorityCenters.begin();

		while(unassignedVertices)
		{
			// the vector must be reseted for this patch
			std::fill(vertexVisited.begin(), vertexVisited.end(), false);
			// the queue of vertices to process must be empty
			assert( verticesToVisit.empty() );

			// 1 - Find a new patch Center
			int new_Center_id = -1;
			while( (nextPriority!= priorityCenters.end()) && (vertex_patch[*nextPriority] != -1)) nextPriority++;
			if(nextPriority != priorityCenters.end())  new_Center_id = *nextPriority;
			else new_Center_id = rand()%numVertices;
			while(vertex_patch[new_Center_id] != -1) new_Center_id = (new_Center_id+1)%numVertices;


			// 2 - Assign the center to the new patch
			patch_centers_indices.push_back(new_Center_id);
			int patchId = patch_centers_indices.size() - 1;


			// 3 - Front propagation
			verticesToVisit.push(std::make_pair(new_Center_id, 0));
			vertexVisited[new_Center_id] = true;

			while( !verticesToVisit.empty() )
			{
				// pop a vertice from the que
				int visited_index = verticesToVisit.front().first;
				int visited_distance = verticesToVisit.front().second;
				verticesToVisit.pop();

				if(visited_distance <= geoDistance[visited_index] ) //if this new patch is better
				{
					// if this vertex was previously unassigned we can decrement the count of unassigned patches
					if(vertex_patch[visited_index] == -1) unassignedVertices--;

					// lets write ourselves as the best patch
					vertex_patch[visited_index] = patchId;
					geoDistance[visited_index] = visited_distance;

					// if we are still close enough, add our neighbours to be processed
					if(visited_distance < maxPatchSize)
					{
						IndexIterator nbegin 	= IndexedMeshAdjacency_Wrapper::neighboursCWBegin(M, visited_index);
						IndexIterator nend		= IndexedMeshAdjacency_Wrapper::neighboursCWEnd(M, visited_index);
						for(IndexIterator ni_itr = nbegin; ni_itr != nend; ++ni_itr)
						{
							if(!vertexVisited[*ni_itr])
							{
								verticesToVisit.push(std::make_pair(*ni_itr, visited_distance + 1));
								vertexVisited[*ni_itr] = true; //just to make sure we only add it once to the queue
							}
						}
					}
				}
			}
		}
		std::cout<<"GeoPatcher returned "<<patch_centers_indices.size()<<" patches"<<std::endl;
	}













	// ##################################
	// ##################################
	// BUILD PATCHADJACENCY
	// ##################################
	// ##################################
	//~ template <class IndexedMeshAdjacency_Wrapper>
	//~ void buildPatchAdjacencyVec(const typename IndexedMeshAdjacency_Wrapper::Data&  M,
	                            //~ const std::vector<int>&                             vertex_patch,
	                            //~ const std::vector<int>&                             patch_vertex,
	                            //~ const std::vector<int>&                             patch_vertex_boundaries, // has size Numpatches+1
	                            //~ std::vector<int>&                                   patch_adj,
	                            //~ std::vector<int>&                                   patch_adj_boundaries)
	//~ {
		//~ typedef typename IndexedMeshAdjacency_Wrapper::IndexIterator       IndexIterator;
		//~ typedef std::vector<int>::iterator                                 IntIterator;
		//~ typedef std::vector<std::vector<int> >::iterator                   IntVecIterator;

		//~ int numPatches = std::max(0, int(patch_vertex_boundaries.size()) - 1);

		//~ patch_adj.clear();
		//~ patch_adj.reserve(4*numPatches);
		//~ patch_adj_boundaries.resize(numPatches+1);
		//~ int currId = 0;

		//~ // build the adjacency
		//~ std::vector<int> tempAdj(numPatches);
		//~ for (int pi = 0; pi < numPatches; ++pi)
		//~ {
			//~ std::fill(tempAdj.begin(), tempAdj.end(), 0);

			//~ int np_begin = patch_vertex_boundaries[pi];
			//~ int np_end = patch_vertex_boundaries[pi+1];

			//~ for(int v = np_begin; v != np_end; ++v)
			//~ {
				//~ int id                 = patch_vertex[v];
				//~ IndexIterator nbegin   = IndexedMeshAdjacency_Wrapper::neighboursCWBegin(M, id);
				//~ IndexIterator nend     = IndexedMeshAdjacency_Wrapper::neighboursCWEnd(M,   id);
				//~ for(IndexIterator ni_itr = nbegin; ni_itr != nend; ++ni_itr)
				//~ {
					//~ tempAdj[vertex_patch[*ni_itr]] = 1;
				//~ }
			//~ }
			//~ patch_adj_boundaries[pi] = currId;
			//~ for(int pi_n=0;pi_n<numPatches;++pi_n)
			//~ {
				//~ if((pi_n != pi) && (tempAdj[pi_n]))
				//~ {
					//~ patch_adj.push_back(pi_n);
					//~ currId++;
				//~ }
			//~ }
			//~ patch_adj_boundaries[pi+1] = currId;
		//~ }
	//~ }


	template <class IndexedMeshAdjacency_Wrapper>
	void buildPatchAdjacencyVec(const typename IndexedMeshAdjacency_Wrapper::Data&  M,
	                            const std::vector<int>&                             vertex_patch,
	                            const std::vector<int>&                             patch_vertex,
	                            const std::vector<int>&                             patch_vertex_boundaries, // has size Numpatches+1
	                            std::vector<int>&                                   patch_adj,
	                            std::vector<int>&                                   patch_adj_boundaries)
	{
		typedef typename IndexedMeshAdjacency_Wrapper::IndexIterator       IndexIterator;
		typedef std::vector<int>::iterator                                 IntIterator;
		typedef std::vector<std::vector<int> >::iterator                   IntVecIterator;

		int numPatches = std::max(0, int(patch_vertex_boundaries.size()) - 1);

		patch_adj.clear();
		patch_adj.reserve(4*numPatches);
		patch_adj_boundaries.resize(numPatches+1);


		// build the adjacency
		std::vector<int> tempAdj(numPatches);
		for (int pi = 0; pi < numPatches; ++pi)
		{
			std::fill(tempAdj.begin(), tempAdj.end(), 0);

			const std::vector<int>::const_iterator vBegin = patch_vertex.begin() + patch_vertex_boundaries[pi];
			const std::vector<int>::const_iterator vEnd   = patch_vertex.begin() + patch_vertex_boundaries[pi+1];

			for(std::vector<int>::const_iterator v_itr = vBegin; v_itr != vEnd; ++ v_itr)
			{
				int id                 = *v_itr;
				IndexIterator nbegin   = IndexedMeshAdjacency_Wrapper::neighboursCWBegin(M, id);
				IndexIterator nend     = IndexedMeshAdjacency_Wrapper::neighboursCWEnd(M,   id);
				for(IndexIterator ni_itr = nbegin; ni_itr != nend; ++ni_itr) {
					tempAdj[vertex_patch[*ni_itr]] = 1;
				}
			}


			patch_adj_boundaries[pi] = patch_adj.size();

			for(int pi_n=0; pi_n<numPatches; ++pi_n) {
				if((pi_n != pi) && (tempAdj[pi_n] != 0)) patch_adj.push_back(pi_n);
			}

		}
		patch_adj_boundaries[numPatches] = patch_adj.size();
	}
















	// ##################################
	// ##################################
	// BUILD SMOOTH PATCHING WEIGHTS
	// ##################################
	// ##################################
	template <class IndexedMeshAdjacency_Wrapper>
	void buildSmoothPatchingWeights ( const typename IndexedMeshAdjacency_Wrapper::Data&  M,
	                                  const double alpha,
	                                  const int maxPatchSize,
	                                  const std::vector<int>& patch_centers_indices,
	                                  const std::vector<int>& patch_vertex,
	                                  const std::vector<int>& patch_vertex_boundaries,
	                                  const std::vector<int>& patch_adj,
	                                  const std::vector<int>& patch_adj_boundaries,
	                                  std::vector<double>&    patch_vertex_smooth_weights,
	                                  std::vector<int>&       patch_vertex_smooth_boundaries )
	{
		typedef typename IndexedMeshAdjacency_Wrapper::IndexIterator       IndexIterator;

		int numPatches = patch_centers_indices.size();
		int numVertices = patch_vertex.size();
		assert( int(patch_vertex_boundaries.size()) == numPatches+1 );
		assert( int(patch_adj_boundaries.size()) == numPatches+1 );

		// 1 - count the number of smooth vertex assignments for each patch (basically its own vertices + those of its neighbours)
		int numSmoothVertices = 0;
		for(int pi = 0; pi < numPatches; ++ pi)
		{
			numSmoothVertices += patch_vertex_boundaries[pi+1] - patch_vertex_boundaries[pi];

			const std::vector<int>::const_iterator nBegin = patch_adj.begin() + patch_adj_boundaries[pi];
			const std::vector<int>::const_iterator nEnd = patch_adj.begin() + patch_adj_boundaries[pi+1];

			for(std::vector<int>::const_iterator ni_itr = nBegin; ni_itr != nEnd; ++ ni_itr)
			{
				int ni = *ni_itr;
				numSmoothVertices += patch_vertex_boundaries[ni+1] - patch_vertex_boundaries[ni];
			}
		}

		// -----------------------------------
		// 2 - allocate the space
		patch_vertex_smooth_weights.resize(numSmoothVertices);
		patch_vertex_smooth_boundaries.resize(numPatches+1);

		// -----------------------------------
		// 3 - compute the weights
		std::vector<double>::iterator smooth_weights_itr = patch_vertex_smooth_weights.begin();
		int distanceLimit = 3*maxPatchSize; // this is an upper boundary on the vertex distance of our neighbours
			// workspace
		std::queue<std::pair<int, int> >    verticesToVisit;
		std::vector<int>                    geoDistance(numVertices);
		std::vector<bool>                   vertexVisited(numVertices);

		for(int pi = 0; pi < numPatches; ++ pi)
		{
			// write the boundary
			patch_vertex_smooth_boundaries[pi] = std::distance(patch_vertex_smooth_weights.begin(), smooth_weights_itr);

			// 3.a - compute geodistances to our patchcenter
			// clean the workspace
			std::fill(geoDistance.begin(), geoDistance.end(), std::numeric_limits<int>::max());
			std::fill(vertexVisited.begin(), vertexVisited.end(), false);
			//verticesToVisit.clear();

			// compute the geodistances to us on a region biiger than us + neighbours
			verticesToVisit.push(std::make_pair(patch_centers_indices[pi], 0));
			vertexVisited[patch_centers_indices[pi]] = true;

			while( !verticesToVisit.empty() )
			{
				// pop a vertice from the queue
				int visited_index = verticesToVisit.front().first;
				int visited_distance = verticesToVisit.front().second;
				verticesToVisit.pop();

				if(visited_distance <= distanceLimit ) //if this new patch is better
				{
					// lets write ourselves as the best patch
					geoDistance[visited_index] = visited_distance;

					// if we are still close enough, add our neighbours to be processed
					if(visited_distance < distanceLimit)
					{
						IndexIterator nbegin 	= IndexedMeshAdjacency_Wrapper::neighboursCWBegin(M, visited_index);
						IndexIterator nend		= IndexedMeshAdjacency_Wrapper::neighboursCWEnd(M, visited_index);
						for(IndexIterator ni_itr = nbegin; ni_itr != nend; ++ni_itr)
						{
							if(!vertexVisited[*ni_itr])
							{
								verticesToVisit.push(std::make_pair(*ni_itr, visited_distance + 1));
								vertexVisited[*ni_itr] = true; //just to make sure we only add it once to the queue
							}
						}
					}
				}
			}


			// 3.b fill the actual weights
			{
				const std::vector<int>::const_iterator vBegin = patch_vertex.begin() + patch_vertex_boundaries[pi];
				const std::vector<int>::const_iterator vEnd = patch_vertex.begin() + patch_vertex_boundaries[pi+1];

				for(std::vector<int>::const_iterator v_itr = vBegin; v_itr != vEnd; ++v_itr) {
					*smooth_weights_itr++ = exp(-double(geoDistance[*v_itr]) / double(maxPatchSize));
				}
			}

			const std::vector<int>::const_iterator nBegin = patch_adj.begin() + patch_adj_boundaries[pi];
			const std::vector<int>::const_iterator nEnd = patch_adj.begin() + patch_adj_boundaries[pi+1];

			for(std::vector<int>::const_iterator ni_itr = nBegin; ni_itr != nEnd; ++ ni_itr)
			{
				int ni = *ni_itr;
				const std::vector<int>::const_iterator vBegin = patch_vertex.begin() + patch_vertex_boundaries[ni];
				const std::vector<int>::const_iterator vEnd = patch_vertex.begin() +patch_vertex_boundaries[ni+1];

				for(std::vector<int>::const_iterator v_itr = vBegin; v_itr != vEnd; ++v_itr) {
					*smooth_weights_itr++ = exp(-(alpha*geoDistance[*v_itr]) / double(maxPatchSize));
				}
			}

		}

		patch_vertex_smooth_boundaries[numPatches] = std::distance(patch_vertex_smooth_weights.begin(), smooth_weights_itr);
		assert( patch_vertex_smooth_boundaries[numPatches] == numSmoothVertices);
	}


} // end namespace PatchedCloud





#endif
