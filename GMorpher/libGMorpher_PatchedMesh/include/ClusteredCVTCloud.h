/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef CLUSTEREDCVTCLOUD_H_DEFINED
#define CLUSTEREDCVTCLOUD_H_DEFINED
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <CC3D/float3.h>
#include <CC3D/rigidT.h>

namespace CVTClusteredCloud
{

	using CC3D::float3;
	using CC3D::quaternion;
	using CC3D::rigidT;
	using namespace Eigen;	


	class CVTClusteredCloud
	{
		public :
			
			CVTClusteredCloud(const std::string& TEMPLATE_FOLDER, const std::vector<float3>& coords_ref, const double avgEdgeL);
			

			void reinit(const std::string& filename, const std::vector<float3>& coords,
						 bool recomputeWeights = false,
						 double distVariance = 0.);

		

			inline const std::vector<float3>&	X0()					const { return m_X0; }		
			inline const std::vector<rigidT>&	RT0()					const { return m_RT0; }
			inline const std::vector<float3>&	dX0_smooth()			const { return m_dX0_smooth; }		
			inline const std::vector<float>&	w_smooth()				const { return m_w_smooth; }		

	
			inline int							numPatches()			const { return m_pv_bounds.size() -1; }
			inline int							numVertices()			const { return m_X0.size(); }
			inline int							numVertices_smooth()	const { return m_dX0_smooth.size(); }	


			inline const std::vector<int>&		vp()					const { return m_vp; }
			inline const std::vector<int>&		pv()					const { return m_pv; }
			inline const std::vector<int>&		pv_bounds()				const { return m_pv_bounds; }
			inline const std::vector<int>&		sadj()					const { return m_sadj; }
			inline const std::vector<int>&		sadj_bounds()			const { return m_sadj_bounds; }
			inline const std::vector<int>&		smooth_bounds()			const { return m_smooth_bounds; }
		

			void blend( const std::vector<rigidT>& RT, std::vector<float3>& X) const;
			void RTfromCloud(const std::vector<float3>& X, std::vector<rigidT>& RT) const;
			void load_vp(const std::string& filename, const int nVertices);
			void load_sadj(const std::string& filename, const int nPatches);
		

		protected :
			// coord data
			std::vector<float3> m_X0;			
			std::vector<rigidT> m_RT0;
			std::vector<float>  m_w_smooth;
			std::vector<float3> m_dX0_smooth;
			// vertex-patch and patch-vertex
			std::vector<int>    m_vp;
			std::vector<int>    m_pv;
			std::vector<int>    m_pv_bounds;
			// sorted patch adj
			std::vector<int>    m_sadj;
			std::vector<int>    m_sadj_bounds;
			// the smooth boundaries
			std::vector<int>    m_smooth_bounds;

};









// ######################################################
// ######################################################
// BUILDING
// ######################################################
// ######################################################

void build_dXN0( int                                         numPatches,
                 const rigidT*                               RT,
                 const int*                                  pv_bounds,
                 const float3*                               XN0,
                 float3*                                     dXN0 );

void build_dX0_smooth( int                        numPatches,
                       const int*                 sadj,
                       const int*                 sadj_bounds,
                       const int*                 pv_bounds,
                       const rigidT*              RT0,
                       const float3*              X0,
                       float3*                    dX0_smooth );



// ######################################################
// ######################################################
// INTERPOLATING
// ######################################################
// ######################################################

void build_w_smooth( int                     numPatches,
                     const double            distVariance,
                     const int*              sadj_bounds,
                     const int*              pv_bounds,
                     const int*              smooth_bounds,
                     const float3*           dX0_smooth,
                     float*                  w_smooth );

void interpolate_cloud_smooth( int           numPatches,
                               const float*  RT,
                               const int*    sadj,
                               const int*    sadj_bounds,
                               const int*    pv_bounds,
                               const float3* dX0_smooth,
                               const float*  w_smooth,
                               float3*       X );
//void interpolate_Cloudsmooth( int           numPatches,
//                              const int*    sadj_bounds, // know how many neighbours
//                              const int*    pv_bounds,   // know how many vertices
//                              const float3* XN_smooth,
//                              const float*  w_smooth,
//                              float3*       X);


// ######################################################
// ######################################################
// TRANSFORMING
// ######################################################
// ######################################################
void transform_Cloud( int           numPatches,
                      const float*  RT,
                      const int*    pv_bounds,
                      const float3* dXN0,
                      float3*       XN );

void rotate_Cloud( int           numPatches,
                   const float*  RT,
                   const int*    pv_bounds,
                   const float3* dXN0,
                   float3*       dXN );

void transform_Cloudsmooth( int           numPatches,
                            const float*  RT,
                            const int*    sadj,
                            const int*    sadj_bounds,
                            const int*    pv_bounds,
                            const float3* dXN0_smooth,
                            float3*       XN_smooth );
void transform_Cloudsmooth( int           numPatches,
                           const float*  RT,
						   const int*	smooth_bound,
                           const int*    sadj,
                           const int*    sadj_bounds,
                           const int*    pv_bounds,
                           const float3* dXN0_smooth,
						   float3*       XN_smooth );

void rotate_dX_smooth( int                     numPatches,
                       const float*            RTf,
                       const int*              sadj,
                       const int*              sadj_bounds,
                       const int*              pv_bounds,
                       const float3*           dX0_smooth,
                       float3*                 dX_smooth );

void rotate_dX_smooth( int                     numPatches,					   
                       const float*            RTf,
					   const int*			   smooth_bound,
                       const int*              sadj,
                       const int*              sadj_bounds,
                       const int*              pv_bounds,
                       const float3*           dX0_smooth,
                       float3*                 dX_smooth );



} // end namespace CVTClusteredCloud


#endif
