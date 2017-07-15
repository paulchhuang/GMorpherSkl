/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#include <RiggedMesh.h>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <cassert>

namespace SklRgedPatchedCloud
{


using IndexedMesh3D::Triangle;
using IndexedMesh3D::float3;
using IndexedMesh3D::ColorVec;

using CC3D::make_float3;

RiggedMesh::RiggedMesh( const std::vector<JointNode>&     KBTree,
                        const std::vector<std::string>& KBTreeNames,
                        const std::string&              meshFile,
                        const char*                     vjFile, 
						int								numSub)
{
	// 1 - load mesh
	std::vector<Triangle>                triangles0;
	std::vector<float3>                  coords0;
	std::vector<ColorVec>                colors0;
	loadIndexedMesh3D(meshFile, coords0, colors0, triangles0);

	

	int numVertices = coords0.size();
	int numJoints   = KBTree.size();

	// 2 - load the vjFile and build the vj vector
	std::ifstream fin(vjFile);
	if( !fin.is_open() ) throw std::runtime_error(vjFile);
	m_vj.resize(numVertices);

	if (numSub > 1){
		for( int vi=0;vi<numVertices;++vi) {
			int buff;
			fin >> buff;			
			m_vj[vi] = buff;
		}
	}
	else{
		for( int vi=0;vi<numVertices;++vi) {
			std::string buff;
			fin >> buff;
			std::vector<std::string>::const_iterator f_itr = std::find( KBTreeNames.begin(), KBTreeNames.end(), buff );
			assert( f_itr != KBTreeNames.end() );
			m_vj[vi] = f_itr - KBTreeNames.begin();
		}
	}		
	fin.close();
	// 3 - build the jv_bounds vector
	m_jv_bounds = std::vector<int>(numJoints+1, 0);
	for( int vi=0;vi<numVertices;++vi) m_jv_bounds[m_vj[vi]+1]++;
	for( int ji=0;ji<numJoints;++ji) m_jv_bounds[ji+1] += m_jv_bounds[ji];

	// 4 - builds the jv_vector
	m_jv.resize( numVertices );
	std::vector<int> offset = m_jv_bounds;
	std::vector<int> newid( numVertices );
	for( int vi=0;vi<numVertices;++vi) {
		int ji = m_vj[vi];
		m_jv[offset[ji]] = vi;
		newid[vi] = offset[ji];
		offset[ji]++;
	}

	// 5 -reindexes triangles, coords, colors
	m_triangles0 = triangles0; //Paul: save the original topology, for the further saving usage
	m_triangles = triangles0;
	
	for( std::vector<IndexedMesh3D::Triangle>::iterator t_itr = m_triangles.begin(); t_itr != m_triangles.end(); ++t_itr) {
		t_itr->v0 = newid[t_itr->v0];
		t_itr->v1 = newid[t_itr->v1];
		t_itr->v2 = newid[t_itr->v2];
	}

	m_coords.resize(numVertices);	
	for(int vi=0;vi<numVertices;++vi) {
		m_coords[vi] = coords0[m_jv[vi]];
		//if (m_jv[vi] == 0){ std::cout << "vi: " << vi << std::endl; return; }

	}


	m_colors.resize(numVertices);	
	for(int vi=0;vi<numVertices;++vi) {
		m_colors[vi] = colors0[m_jv[vi]];
	}
	

}




void RiggedMesh::transformCloud( const double*                       TTis,
                                 std::vector<IndexedMesh3D::float3>& coords_out )
{
	coords_out.resize( numVertices() );
	
	for (int ji = 0; ji < numJoints(); ++ji) {
		const double* Ri =  TTis + 12*ji;
		const double* ci =  TTis + 12*ji+9;
		float3 cif = make_float3(ci[0], ci[1], ci[2]);
		float Rif[9];
		std::copy(Ri, Ri+9, Rif);

		std::vector<float3>::const_iterator xin_itr = m_coords.begin() + m_jv_bounds[ji];
		std::vector<float3>::const_iterator xin_end = m_coords.begin() + m_jv_bounds[ji+1];
		std::vector<float3>::iterator      xout_itr = coords_out.begin() + m_jv_bounds[ji];
		
		while( xin_itr != xin_end )  {
			*xout_itr++ = prod3(Rif, *xin_itr++) + cif;
		}
	}
}


} // end namespace Kineben 
