/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef RIGGEDMESH_H_DEFINED
#define RIGGEDMESH_H_DEFINED

#include <IndexedMesh/IndexedMesh3D.h>
#include "Skl.h"

namespace SklRgedPatchedCloud
{
	
	class RiggedMesh 
	{
		public :

		RiggedMesh(const std::vector<JointNode>&     KBTree,
					const std::vector<std::string>& KBTreeNames,
					const std::string&              meshFile,
					const char*                     vjFile,
					int								numSub);

		inline int numJoints() const { return m_jv_bounds.size() -1; }
		inline int numVertices() const { return m_coords.size(); }

		inline const std::vector<int>& vj()        const { return m_vj; }
		inline const std::vector<int>& jv()        const { return m_jv; }
		inline const std::vector<int>& jv_bounds() const { return m_jv_bounds; }

		inline const std::vector<IndexedMesh3D::Triangle>& triangles() const { return m_triangles; }
		inline const std::vector<IndexedMesh3D::Triangle>& triangles0() const { return m_triangles0; }
		inline const std::vector<IndexedMesh3D::float3>&   coords() const { return m_coords; }
		inline const std::vector<IndexedMesh3D::ColorVec>& colors() const { return m_colors; }


		void transformCloud( const double*                       TTis,
							 std::vector<IndexedMesh3D::float3>& coords_out );

		protected :
		std::vector<int> m_vj;
		std::vector<int> m_jv;
		std::vector<int> m_jv_bounds;

		std::vector<IndexedMesh3D::Triangle>  m_triangles;
		std::vector<IndexedMesh3D::Triangle>  m_triangles0;
		std::vector<IndexedMesh3D::float3>    m_coords;
		std::vector<IndexedMesh3D::ColorVec>  m_colors;
	};

} 

#endif
