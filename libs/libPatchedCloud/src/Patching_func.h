/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef PATCHING_FUNC_H_DEFINED
#define PATCHING_FUNC_H_DEFINED

#include <vector>

namespace SklRgedPatchedCloud
{
	// ------------------------------------------
	// seeds random patches on the surface of a maximal radius
	template <class IndexedMeshAdjacency_Wrapper>
	void SeedPatches( const typename IndexedMeshAdjacency_Wrapper::Data& M, const int maxPatchSize,
	                  std::vector<int> &vertex_patch,
	                  std::vector<int> &patch_centers_indices,
	                  std::vector<int> &geoDistance,
	                  const std::vector<int>&  priorityCenters = std::vector<int>());


	// ------------------------------------------
	// computes the adjacency of the patches
	template <class IndexedMeshAdjacency_Wrapper>
	void buildPatchAdjacencyVec(const typename IndexedMeshAdjacency_Wrapper::Data&  M,
	                            const std::vector<int>&                             vertex_patch,
	                            const std::vector<int>&                             patch_vertex,
	                            const std::vector<int>&                             patch_vertex_boundaries, // has size Numpatches+1
	                            std::vector<int>&                                   patch_adj,
	                            std::vector<int>&                                   patch_adj_boundaries);

	// ------------------------------------------
	// builds the patch_vertex vector from the vertex_patch
	void VP_To_PV(int numPatches,
	              const std::vector<int>&         vertex_patch,
	              std::vector<int>&               patch_vertex,
	              std::vector<int>&               patch_vertex_boundaries);

	// ------------------------------------------
	// removes any edge (i,j) where i > j and makes sure that if there is one we have (j,i) stored already
	void orientAdjVec( const std::vector<int>&  adj,
	                   const std::vector<int>&  adj_boundaries,
	                   std::vector<int>&        oadj,
	                   std::vector<int>&        oadj_boundaries);


	template <class IndexedMeshAdjacency_Wrapper>
	void buildSmoothPatchingWeights ( const typename IndexedMeshAdjacency_Wrapper::Data&  M,
	                                  const double alpha,
	                                  const int maxPatchSize,
	                                  const std::vector<int>& patchesCenter,
	                                  const std::vector<int>& patch_vertex,
	                                  const std::vector<int>& patch_vertex_boundaries,
	                                  const std::vector<int>& patch_adj,
	                                  const std::vector<int>& patch_adj_boundaries,
	                                  std::vector<double>&    patch_vertex_smooth_weights,
	                                  std::vector<int>&       patch_vertex_smooth_boundaries );





	int getConnectedComponents( const std::vector<int>& sorted_patch_adj,
	                            const std::vector<int>& sorted_patch_adj_boundaries,
	                            std::vector<int>&       patch_component );





} // end namespace PatchedCloud

#endif
