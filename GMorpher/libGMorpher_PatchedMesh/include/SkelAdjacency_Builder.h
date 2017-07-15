/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/

#ifndef SKELADJACENCY_BUILDER_H_DEFINED
#define SKELADJACENCY_BUILDER_H_DEFINED


#include "PatchedCloud.h"
#include <vector>
#include <string>

namespace SklRgedPatchedCloud
{

/**
	The idea behind this is that you will simulate an articulated rigidity by 
	- loading the mesh
	- making a patched mesh
	- building a skel adjacency
	- create a solver with the skel adjacency
	- creating a Energy_Rigidity_PAULY with that skel adjacency 
*/



// 1 - load the vj vector and make it numeric 
// 2 - vote for each patch to find the majority bone and assign it 
// 3 - cluster the patches together by bones and write the adjacency
void buildSkelAdjacency( const std::string&      filename,
                         const PatchedCloud&     PC,
                         std::vector<int>&       skelAdj,
                         std::vector<int>&       skelAdj_bounds );

void mergeSkelAdjacency( const std::vector<int>& skelAdj,
                         const std::vector<int>& skelAdj_bounds,
                         const std::vector<int>& patchAdj,
                         const std::vector<int>& patchAdj_bounds,
                         std::vector<int>&       mergeAdj,
                         std::vector<int>&       mergeAdj_bounds );



} // end namespace PatchedCloud

#endif
