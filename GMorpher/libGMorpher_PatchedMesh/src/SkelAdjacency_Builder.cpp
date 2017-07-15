#include <SkelAdjacency_Builder.h>
#include <stdexcept>
#include <fstream>
#include <set>
#include <map>
#include <algorithm>
#include <cassert>

namespace SklRgedPatchedCloud
{


int load_Patch_to_Joint( const std::string&   filename,
                         const PatchedCloud&  PC,
                         std::vector<int>&    patch_joint )
{
	// --------------------
	// 1 - load the vj vector and make it numeric 
	std::ifstream fin( filename.c_str() );
	if( !fin.is_open() ) 
		            throw std::runtime_error("(E): could not open " + filename);

		// check that we have the correct vertex count
	int numVertices = PC.numVertices();

		// read the file
	std::vector<std::string> vj_string(numVertices);
	std::set<std::string> joint_names;
	for(int vi=0;vi<numVertices;++vi) {
		std::string name;
		fin >> vj_string[vi];
		joint_names.insert( vj_string[vi] );
	}
	int numJoints = joint_names.size();

		// convert the names to indices
	std::map<std::string, int> name_to_id;
	{
		int id = 0;
		for(std::set<std::string>::const_iterator n_itr = joint_names.begin();
			                               n_itr != joint_names.end(); ++n_itr){
			name_to_id[*n_itr] = id++;
		}
	}
	std::vector<int> vj(numVertices);
	for(int vi=0;vi<numVertices;++vi) {
		vj[vi] = name_to_id[vj_string[vi]];
	}


	// --------------------
	// 2 - vote for each patch to find the majority bone and assign it
	int numPatches = PC.numPatches();
	patch_joint.resize(numPatches);
	for(int pi=0;pi<numPatches;++pi)
	{
		std::vector<int>::const_iterator vi_itr = PC.pv().begin() 
		                                                   + PC.pv_bounds()[pi];
		std::vector<int>::const_iterator vi_end = PC.pv().begin() 
		                                                 + PC.pv_bounds()[pi+1];

		std::vector<int> vCount(numJoints, 0);
		while( vi_itr != vi_end ) {
			int jointId = vj[*vi_itr++];
			vCount[jointId]++;
		}

		patch_joint[pi] = std::max_element(vCount.begin(), vCount.end() ) 
		                                                       - vCount.begin();
	}

	return numJoints;
}



void buildSkelAdjacency( const std::string&      filename,
                         const PatchedCloud&     PC,
                         std::vector<int>&       skelAdj,
                         std::vector<int>&       skelAdj_bounds )
{
	// ----------
	// first get patch_joint
	std::vector<int> patch_joint;
	int numJoints  = load_Patch_to_Joint( filename, PC, patch_joint );
	int numPatches = PC.numPatches();
	
	// --------------------
	// create joint_patch and joint_patch_bounds
	std::vector<int> joint_patch;
	//joint_patch.reserve(numPatches);		//Paul
	joint_patch.resize(numPatches);
	std::vector<int> joint_patch_bounds(numJoints+1, 0);
	
		// first we count how many patch per joint
	for(int pi=0;pi<numPatches;++pi) joint_patch_bounds[patch_joint[pi]+1] ++;
	for(int ji=0;ji<numJoints;++ji)  joint_patch_bounds[ji+1] += joint_patch_bounds[ji];
	
		// second we write the joint_patch vector 
	std::vector<int> offset = joint_patch_bounds;
	
	//system("pause");
	for(int pi=0;pi<numPatches;++pi) {
		int jointId = patch_joint[pi];
		joint_patch[offset[jointId]] = pi;
		offset[jointId]++;
	}
	//system("pause");
	// --------------------
	//	cluster the patches together by bones and write the adjacency
	skelAdj.clear();
	skelAdj_bounds.resize(numPatches+1);
	skelAdj_bounds[0] = 0;
	
	for(int pi=0;pi<numPatches;++pi) {
		int jointId = patch_joint[pi];
		std::vector<int>::const_iterator pj_itr = joint_patch.begin() + 
		                                            joint_patch_bounds[jointId];
		std::vector<int>::const_iterator pj_end = joint_patch.begin() + 
		                                          joint_patch_bounds[jointId+1];
		while( pj_itr != pj_end ) {
			if( *pj_itr != pi ) skelAdj.push_back(*pj_itr);
			pj_itr++;
		}
		skelAdj_bounds[pi+1] = skelAdj.size();
	}	
}



void mergeSkelAdjacency( const std::vector<int>& skelAdj,
                         const std::vector<int>& skelAdj_bounds,
                         const std::vector<int>& patchAdj,
                         const std::vector<int>& patchAdj_bounds,
                         std::vector<int>&       mergeAdj,
                         std::vector<int>&       mergeAdj_bounds )
{
	assert( skelAdj_bounds.size() == patchAdj_bounds.size() );
	int numPatches = skelAdj_bounds.size() -1 ;

	// allocate output vectors
	mergeAdj.clear();
	mergeAdj.reserve( skelAdj.size() + patchAdj.size() );
	mergeAdj_bounds.resize( numPatches + 1);

	// fill output vectors
	std::set<int> sorter;
	mergeAdj_bounds[0] = 0;
	for(int pi=0;pi<numPatches;++pi) {
		sorter.clear();
		// insert the skel neighbors
		for( std::vector<int>::const_iterator 
		                        n_itr  = skelAdj.begin() + skelAdj_bounds[pi];
		                        n_itr != skelAdj.begin() + skelAdj_bounds[pi+1];
		                        n_itr++ ) sorter.insert( *n_itr);

		// insert the patch neighbors
		for( std::vector<int>::const_iterator 
		                      n_itr  = patchAdj.begin() + patchAdj_bounds[pi];
		                      n_itr != patchAdj.begin() + patchAdj_bounds[pi+1];
		                      n_itr++ ) sorter.insert( *n_itr);

		// read out 
		for(std::set<int>::const_iterator n_itr = sorter.begin(); 
		            n_itr != sorter.end(); ++n_itr ) mergeAdj.push_back(*n_itr);

		// write the bound 
		mergeAdj_bounds[pi+1] = mergeAdj.size();
	}
}


} // end namespace PatchedCloud
