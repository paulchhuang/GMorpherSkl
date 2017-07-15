/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef INDEXEDMESH3D_H_DEFINED
#define INDEXEDMESH3D_H_DEFINED


#include <CC3D/float3.h>
#include "IndexedMesh/ColorVec.h"
#include <vector>

#ifdef WIN32 
	#ifdef INDEXEDMESH_EXPORTS
	#define INDEXEDMESH_API __declspec( dllexport )
	#else
	#define INDEXEDMESH_API __declspec( dllimport )
	#endif
#else 
	#define INDEXEDMESH_API
#endif

namespace IndexedMesh3D
{

	typedef CC3D::float3 float3;
	using CC3D::make_float3;

	struct Triangle {
		int v0, v1, v2;
	};

	struct Edge {
		int            v;
		unsigned char  hasNext;
		int            vnext;
	};


	// ##################################################################
	// IO
	// ##################################################################

	void INDEXEDMESH_API 
	loadIndexedMesh3D( const std::string&     meshFile,
	                   std::vector<float3>&   coords,
	                   std::vector<ColorVec>& colors,
	                   std::vector<Triangle>& triangles);

	void INDEXEDMESH_API 
	saveIndexedMesh3D( const char*                  filename,
	                   const std::vector<float3>&   coords,
	                   const std::vector<ColorVec>& colors,
	                   const std::vector<Triangle>& triangles);


	// ##################################################################
	// triangles to adjacency
	// ##################################################################
	void INDEXEDMESH_API 
	convertTrianglesToEdgeGraph(const int numVertices,
	                            const std::vector<Triangle>& triangles,
	                            std::vector<Edge>& CWEdges,
	                            std::vector<int>&  CWEdgesBoundaries);

	// ##################################################################
	// ajdacency to triangles
	// ##################################################################
	void INDEXEDMESH_API 
	convertEdgeGraphToTriangles( const std::vector<Edge>& CWEdges,
	                             const std::vector<int>&  CWEdgesBoundaries,
	                             std::vector<Triangle>&   triangles);

	// ##################################################################
	// reindex the mesh
	// It would for sure be easier to do this at the Triangle stage... 
	// but graph algorithms operate on adjacency graphs
	// ##################################################################
	// the new vertex at position i had the indice oldIndices[i] in A0
	void INDEXEDMESH_API 
	reindexAdjacency( const std::vector<Edge>& CWEdges0,
	                  const std::vector<int>&  CWEdgesBoundaries0,
	                  const std::vector<int>&  oldIndices, 
	                  std::vector<Edge>&       CWEdges,
	                  std::vector<int>&        CWEdgesBoundaries);




	// ##################################################################
	// count triangles
	// ##################################################################
	int INDEXEDMESH_API 
	getNumTriangles(const std::vector<Edge>& CWEdges,
	                const std::vector<int>&  CWEdgesBoundaries);

	// ##################################################################
	// Extract the Normals
	// ##################################################################
	void INDEXEDMESH_API 
	extractNormals(const std::vector<Edge>&       CWEdges,
	               const std::vector<int>&        CWEdgesBoundaries,
	               const std::vector<float3>&     coords,
	               std::vector<float3>&           normals);

	void INDEXEDMESH_API 
	averageNormals( const std::vector<Triangle>& triangles,
	                const std::vector<float3>&   coords,
	                std::vector<float3>&         normals);
	
	void INDEXEDMESH_API
	averageNormals( const std::vector<Triangle>& triangles,
	                     const std::vector<float3>&   coords,
	                     std::vector<float3>&         normals, 
						 std::vector<float3>&		  normals_tri);

	void INDEXEDMESH_API		//Paul
	averageNormalsFliped( const std::vector<Triangle>& triangles,
	                const std::vector<float3>&   coords,
	                std::vector<float3>&         normals);

	void INDEXEDMESH_API		//Paul
	averageNormalsFliped( const std::vector<Triangle>& triangles,
	                const std::vector<float3>&   coords,
	                std::vector<float3>&         normals, 
					std::vector<float3>&		  normals_tri);

	// ##################################################################
	// GetMeanEdgeLength
	// ##################################################################
	double INDEXEDMESH_API  
	getMeanEdgeLength(const std::vector<Triangle>& triangles,
	                  const std::vector<float3>& coords);

	double INDEXEDMESH_API 
	getMeanEdgeLength(const std::vector<Edge>&     CWEdges,
	                  const std::vector<int>&      CWEdgesBoundaries,
	                  const std::vector<float3>& coords);



	// ##################################################################
	// GetConnectedComponents
	// ##################################################################
	int INDEXEDMESH_API 
	getConnectedComponents( const std::vector<Edge>& CWEdges,
	                        const std::vector<int>&  CWEdgesBoundaries,
	                        std::vector<int>&       vertex_ccomponent );


	// ##################################################################
	// GetBBox
	// ##################################################################
	void INDEXEDMESH_API 
	getBBox(const std::vector<float3>& coords, 
	        double& xmin, double& ymin, double& zmin, 
	        double& xmax, double& ymax, double& zmax);
}



#endif
