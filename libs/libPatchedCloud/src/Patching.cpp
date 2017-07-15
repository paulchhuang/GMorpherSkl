/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#include <Patching.h>
#include "Patching_func.h"
#include "Patching_impl.h"
#include <fstream>
#include <stdexcept>


namespace SklRgedPatchedCloud
{
	class IndexedMesh_Wrapper
	{
		public :

			// there is surely something better to be done here without copy.. with an adapter iterator
		struct Data {
			Data(const std::vector<IndexedMesh3D::Edge>& e, const std::vector<int>& eb) :
			CWEdgesBoundaries(eb)
			{
				CWEdges.resize(e.size());
				for( std::size_t ei =0; ei< e.size(); ++ei) { CWEdges[ei] = e[ei].v; }
			}
			std::vector<int>                 CWEdges;
			std::vector<int>                 CWEdgesBoundaries;
		};

		typedef std::vector<int>::const_iterator   IndexIterator;

		inline static IndexIterator neighboursCWBegin(const Data& A, int id){
			return A.CWEdges.begin() + A.CWEdgesBoundaries[id];
		}
		inline static IndexIterator neighboursCWEnd(const Data& A, int id){
			return A.CWEdges.begin() + A.CWEdgesBoundaries[id+1];
		}
		inline static int           numNeighbours(const Data& A, int id) {
			return A.CWEdgesBoundaries[id+1] - A.CWEdgesBoundaries[id];
		}
		inline static int           numVertices(const Data& A){
			return A.CWEdgesBoundaries.size() - 1;
		}
	};






	int Patching::get_patch_vertex_index( int vertex_index ) const
	{
		int pi = m_vp[vertex_index];
		std::vector<int>::const_iterator piBegin = m_pv.begin() + m_pv_boundaries[pi];
		std::vector<int>::const_iterator piEnd   = m_pv.begin() + m_pv_boundaries[pi+1];

		std::vector<int>::const_iterator f_itr   = std::find(piBegin, piEnd, vertex_index);
		if( f_itr == piEnd) {
			std::cout<<"wrong index : "<<vertex_index<<std::endl;
			throw ( std::runtime_error("inconsistency in the patch_vertex and vertex_patch vectors") );
		}

		return f_itr - m_pv.begin();
	}








	void Patching::init(const int MatchPatchSize,
	                    const std::vector<IndexedMesh3D::Edge> CWEdges,
	                    const std::vector<int>                 CWEdgesBoundaries,
	                    const char*                            vp_filename)
	{
		IndexedMesh_Wrapper::Data A( CWEdges, CWEdgesBoundaries );

		if( vp_filename == 0 ) {
			std::vector<int> p_centers_indices;
			SeedPatches_redux<IndexedMesh_Wrapper>( A, MatchPatchSize, m_vp, p_centers_indices, m_geoDistance);
			m_numPatches = p_centers_indices.size();
		}
		else {
			std::ifstream fvp( vp_filename);
			//if( ! fvp.is_open() ) throw std::runtime_error(std::string("could not open : ") + std::string(vp_filename) ) ;
			size_t numVertices = CWEdgesBoundaries.size() - 1;
			m_vp.resize(numVertices);
			m_numPatches = 0;
			for( size_t vi=0;vi<numVertices;++vi) {
				fvp >> m_vp[vi];
				if( m_vp[vi] > m_numPatches) m_numPatches = m_vp[vi];
			}
			fvp.close();
		}



		//  ----------------
		// the patch_vertex vector
		VP_To_PV(numPatches(), m_vp, m_pv, m_pv_boundaries);
		//  ----------------
		// build and sort the patch adjacency vector
		buildPatchAdjacencyVec<IndexedMesh_Wrapper>(A, m_vp, m_pv, m_pv_boundaries, m_sorted_patch_adj, m_sorted_patch_adj_boundaries);
		for(size_t pi = 0; pi < numPatches(); ++pi ) {
			std::vector<int>::iterator nBegin = m_sorted_patch_adj.begin() + m_sorted_patch_adj_boundaries[pi];
			std::vector<int>::iterator nEnd   = m_sorted_patch_adj.begin() + m_sorted_patch_adj_boundaries[pi+1];
			std::sort(nBegin, nEnd);
		}


		m_smooth_boundaries.resize(numPatches()+1);
		m_smooth_boundaries[0] = 0;
		for(size_t pi = 0; pi < numPatches(); ++pi ) {
			int numVertices_i   = m_pv_boundaries[pi+1] - m_pv_boundaries[pi];
			int numNeighbours_i = m_sorted_patch_adj_boundaries[pi+1] - m_sorted_patch_adj_boundaries[pi];
			m_smooth_boundaries[pi+1] = m_smooth_boundaries[pi] + numVertices_i*(numNeighbours_i+1); // the +1 accounts for ourselves
		}


		//  ----------------
		// The sorted patch edges
		m_sorted_pe.clear();
		m_sorted_pe.reserve(numPatches() + m_sorted_patch_adj.size() );
		{
			for(size_t pi = 0; pi < numPatches(); ++pi )
			{
				int numVertices_i                 = m_pv_boundaries[pi+1] - m_pv_boundaries[pi];
					// us
				m_sorted_pe.push_back( PEdge(pi, pi, numVertices_i, numVertices_i, 0, 0, m_smooth_boundaries[pi], m_smooth_boundaries[pi]) ) ;
					// go through beighbours
				std::vector<int>::const_iterator n_itr  = m_sorted_patch_adj.begin() + m_sorted_patch_adj_boundaries[pi];
				std::vector<int>::const_iterator nEnd   = m_sorted_patch_adj.begin() + m_sorted_patch_adj_boundaries[pi+1];
				while ( (n_itr != nEnd) ) {
					int pj            = *n_itr;
					int numVertices_j = m_pv_boundaries[pj+1] - m_pv_boundaries[pj];
					m_sorted_pe.push_back( PEdge(pi, *n_itr, numVertices_i, numVertices_j, 0, 0, m_smooth_boundaries[pi], m_smooth_boundaries[pj]) ) ;
					n_itr++;
				}
			}

			// fill the TjPis
			for( size_t eij = 1; eij < m_sorted_pe.size(); ++eij){
				m_sorted_pe[eij].TjPiOffset = m_sorted_pe[eij-1].TjPiOffset + m_sorted_pe[eij-1].numVertices_i;
			}

			// fill the TiPjs and TiPis and TjPjs
			for( size_t eij = 0; eij < m_sorted_pe.size(); ++eij){
				int pi = m_sorted_pe[eij].pi;
				int pj = m_sorted_pe[eij].pj;
				// we are at the TjPi ( eij) edge.. we look for the TiPj edge (eji)

				//search for eji, the one whose TiPj is our TjPi
				int eji = m_sorted_patch_adj_boundaries[pj] + pj;
				while( m_sorted_pe[eji].pj != pi ) eji++;

				m_sorted_pe[eji].TiPjOffset = m_sorted_pe[eij].TjPiOffset;
			}
		}


		//  ----------------
		// The oriented patch edges
		int numOriented = m_sorted_patch_adj.size()/2;
		m_oriented_pe.clear();
		m_oriented_pe.reserve(numOriented);
		for( size_t eij = 0; eij < m_sorted_pe.size(); ++eij){
			int pi = m_sorted_pe[eij].pi;
			int pj = m_sorted_pe[eij].pj;
			if( pj > pi ) m_oriented_pe.push_back(m_sorted_pe[eij]);
		}


		//  ----------------
		// copy the connectivity data
		m_CWEdges = CWEdges;
		m_CWEdgesBoundaries = CWEdgesBoundaries;
		IndexedMesh3D::reindexAdjacency( m_CWEdges,
		                                 m_CWEdgesBoundaries,
		                                 m_pv,
		                                 m_pv_CWEdges,
		                                 m_pv_CWEdgesBoundaries);
	}




	void Patching::saveToFile( const char* filename )
	{
		std::ofstream fout(filename);
		//if( ! fout.is_open())  throw std::runtime_error(std::string("could not open : ") + std::string(filename) ) ;
		//mark by Paul
		for(size_t vi=0;vi<m_vp.size();++vi) fout<<m_vp[vi]<<"\n";

		fout.close(); // not necessary
	}









	void Patching::genColors( std::vector<IndexedMesh3D::ColorVec>& colors ) const
	{
		// resize output vector
		colors.resize( numVertices() );

		// generate random colors
		std::vector<IndexedMesh3D::ColorVec> patchcolors( numPatches() );
		for(size_t pi=0; pi< numPatches(); ++pi ){
			patchcolors[pi].r  = ( (rand() %1000)/1000.0);
			patchcolors[pi].g  = ( (rand() %1000)/1000.0);
			patchcolors[pi].b  = ( (rand() %1000)/1000.0);
			patchcolors[pi].a = 1.0;
		}

		// fill the vertex colors

		for(size_t vi=0; vi<numVertices(); ++vi) {
			int pid    = m_vp[vi];
			colors[vi] = patchcolors[pid];
		}

	}



} // end namespace PatchedCloud
