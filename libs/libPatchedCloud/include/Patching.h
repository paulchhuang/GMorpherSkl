/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef PATCHING_H_DEFINED
#define PATCHING_H_DEFINED

#include <IndexedMesh/IndexedMesh3D.h>

namespace SklRgedPatchedCloud
{

	class Patching
	{
		public :

		struct PEdge {
				inline PEdge( const int i, const int j, int ni, int nj, int TiPj, int TjPi, int TiPi, int TjPj):
			pi(i), pj(j), numVertices_i(ni), numVertices_j(nj), TiPjOffset(TiPj), TjPiOffset(TjPi),TiPiOffset(TiPi), TjPjOffset(TjPj){}
			int pi;
			int pj;
			int numVertices_i;
			int numVertices_j;
			int TiPjOffset;
			int TjPiOffset;
			int TiPiOffset;
			int TjPjOffset;
		};


		void init(const int                              MatchPatchSize,
		          const std::vector<IndexedMesh3D::Edge> CWEdges,
		          const std::vector<int>                 CWEdgesBoundaries,
		          const char*                            vp_filename = 0);

		void saveToFile( const char* filename );

		// -------------------------
		// accessors

		/// how many patches are there ?
		inline size_t numPatches()                                         const { return m_numPatches; }
		/// how many vertices are there ?
		inline size_t numVertices()                                        const { return m_vp.size(); }
		/// how many patchEdges are there
		inline size_t numPatchEdges()                                      const { return m_sorted_patch_adj.size(); }
		/// how many smooth vertices are there ?
		inline size_t numSmoothVertices()                                  const { return m_smooth_boundaries[numPatches()]; }

		//inline const std::vector<int>&    patch_centers_indices()          const { return m_patch_centers_indices; }
		inline const std::vector<int>&    vp()                             const { return m_vp; }
		inline const std::vector<int>&    pv()                             const { return m_pv; }
		inline const std::vector<int>&    pv_boundaries()                  const { return m_pv_boundaries; }

		inline const std::vector<int>&    sorted_patch_adj()               const { return m_sorted_patch_adj; }
		inline const std::vector<int>&    sorted_patch_adj_boundaries()    const { return m_sorted_patch_adj_boundaries; }

		inline const std::vector<int>&    smooth_boundaries()              const { return m_smooth_boundaries; }

		inline const std::vector<PEdge>&  sorted_pe()                      const { return m_sorted_pe; }
		inline const std::vector<PEdge>&  oriented_pe()                    const { return m_oriented_pe; }

		// -------------------------
		// connectivity stuff
		inline const std::vector<IndexedMesh3D::Edge>& CWEdges()               const {return m_CWEdges; }
		inline const std::vector<int>&                 CWEdgesBoundaries()     const {return m_CWEdgesBoundaries; }
		inline const std::vector<IndexedMesh3D::Edge>& pv_CWEdges()            const {return m_pv_CWEdges; }
		inline const std::vector<int>&                 pv_CWEdgesBoundaries()  const {return m_pv_CWEdgesBoundaries; }


		/// utility function to get the patch_vertex index of a vector from his original vertex_index
		int get_patch_vertex_index( int vertex_index ) const;
		/// utility function to have colours on patches
		void genColors( std::vector<IndexedMesh3D::ColorVec>& colors ) const;


		protected :
		size_t              m_numPatches;
		// the normal patching stuff
		std::vector<int>    m_vp;
		//std::vector<int>    m_patch_centers_indices;
		std::vector<int>    m_geoDistance;

		// the patch_vertex vector
		std::vector<int>    m_pv;
		std::vector<int>    m_pv_boundaries;

		// the patch adjacency
		std::vector<int>    m_sorted_patch_adj;
		std::vector<int>    m_sorted_patch_adj_boundaries;

		// the smooth boundaries
		std::vector<int>    m_smooth_boundaries;

		// the patch edges
		std::vector<PEdge>  m_sorted_pe;   // all patch edges (i,j), (j,i) and (i,i) sorted
		std::vector<PEdge>  m_oriented_pe; // only patch edges (i,j) where (j>i)


		// the original connectivity
		std::vector<IndexedMesh3D::Edge> m_CWEdges;
		std::vector<int>                 m_CWEdgesBoundaries;
		// the new connectivity
		std::vector<IndexedMesh3D::Edge> m_pv_CWEdges;
		std::vector<int>                 m_pv_CWEdgesBoundaries;

	};
} // end namespace MeshNLMorpher








#endif
