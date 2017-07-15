/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef GRAPHSMOOTH_H_DEFINED
#define GRAPHSMOOTH_H_DEFINED

#pragma warning(disable: 4267) 
#pragma warning(disable: 4244) 
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <CC3D/float3.h>
#include <CC3D/rigidT.h>
#include <Skl.h>
#include <RiggedMesh.h>
#include "Patching.h"

namespace SklRgedPatchedCloud
{

	using CC3D::float3;
	using CC3D::quaternion;
	using CC3D::rigidT;
	using IndexedMesh3D::ColorVec;
	using namespace Eigen;	


	class PatchedCloud
	{
		public :
		PatchedCloud( const Patching&            patches,
					  double                     distVariance,
					  const std::vector<float3>& coords,
					  const std::vector<float3>& normals, 
					  const std::vector<KBNode>&   KBNodes,
					  const std::vector<KBJoint>&  KBJoints, 
					  const RiggedMesh&   mesh);
		PatchedCloud(   const Patching&				 patches,
	                    double						 distVariance,
                        const std::vector<float3>&	 coords,
                        const std::vector<float3>&	 normals);
		PatchedCloud( const Patching&				 patches,
                        double						 distVariance,
                        const std::vector<float3>&	 coords,
                        const std::vector<float3>&	 normals, 
						const std::vector<KBNode>&   KBNodes,
						const std::vector<KBJoint>&  KBJoints,
						const RiggedMesh&   mesh,		
						const std::vector<float3>&	 refVectors0);

		void reinit( const std::vector<float3>& coords,
					 const std::vector<float3>& normals,
					 bool recomputeWeights = false,
					 double distVariance = 0.);

		void sklinit(	const Patching&				patches,
						const std::vector<KBNode>&  KBNodes, 
						const std::vector<KBJoint>&  KBJoints,
						const RiggedMesh&  mesh, 
						bool recomputeWeights =		false);

		inline const std::vector<float3>&	X0()					const { return m_X0; }
		inline const std::vector<float3>&	N0()					const { return m_N0; }
		inline const std::vector<rigidT>&	RT0()					const { return m_RT0; }
		inline const std::vector<float3>&	dX0_smooth()			const { return m_dX0_smooth; }
		inline const std::vector<float3>&	dJX0_smooth()			const { return m_dJX0_smooth; }
		inline const std::vector<float>&	w_smooth()				const { return m_w_smooth; }
		inline const std::vector<float>&	w_JX_smooth()			const { return m_w_JX_smooth; }

		inline int							numJoints()				const { return m_joint_coords0.size(); }		// ROOT_T included
		inline int							numPatches()			const { return m_pv_bounds.size() -1; }
		inline int							numVertices()			const { return m_X0.size(); }
		inline int							numVertices_smooth()	const { return m_dX0_smooth.size(); }	


		inline const std::vector<int>&		vp()					const { return m_vp; }
		inline const std::vector<int>&		pv()					const { return m_pv; }
		inline const std::vector<int>&		pv_bounds()				const { return m_pv_bounds; }
		inline const std::vector<int>&		sadj()					const { return m_sadj; }
		inline const std::vector<int>&		sadj_bounds()			const { return m_sadj_bounds; }
		inline const std::vector<int>&		smooth_bounds()			const { return m_smooth_bounds; }
		inline const std::vector<int>&		patch_to_joint()		const { return m_patch_to_joint; }
		inline const std::vector<int>&		patch_to_joint_child()	const { return m_patch_to_joint_child; }

		inline const std::vector<int>&		patch_comp()			const { return m_patch_comp;}

		inline const std::vector<float3>&	joint_coords0()			const { return m_joint_coords0; }
		inline const std::vector<float3>&	ppdicular_ft0()			const { return m_ppdicular_ft0; }
		inline const std::vector<float3>&	beta0()					const { return m_beta0; }
		inline const std::vector<float>&	w_gamma()				const { return m_w_gamma; }
		inline const std::vector<float>&	w_kappa()				const { return m_w_kappa; }

		inline const std::vector<int>&		jOffsp()				const { return m_jOffsp; }
		inline const std::vector<int>&		jOffsp_bounds()			const { return m_jOffsp_bounds; }
		inline const std::vector<int>&		jp()					const { return m_jp; }
		inline const std::vector<int>&		jp_bounds()				const { return m_jp_bounds;}
		
		inline const std::vector<float3>&	ref_vec()				const { return m_refVectors_LCF;}

		void blend( const std::vector<rigidT>& RT, std::vector<float3>& X) const;
		void RTfromCloud( const std::vector<float3>& X, std::vector<rigidT>& RT) const;
		void genJointBasedColor( std::vector<ColorVec>& colors ) const;
		void JXfromRT_Bone(const std::vector<rigidT>& RT, std::vector<float3>& JX);
		void JXfromRT_Joint(const std::vector<rigidT>& RT, std::vector<float3>& JX);

		protected :
		// coord data
		std::vector<float3> m_X0;
		std::vector<float3> m_N0;
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


		// the skeleton thing		
		std::vector<float3> m_joint_coords0;		// length: 16, 32..., ROOT_T included
		std::vector<float3> m_ppdicular_ft0;		// length: numPatch
		std::vector<float3> m_beta0;
		std::vector<float3> m_dJX0_smooth;
		std::vector<float>	m_w_JX_smooth;
		std::vector<float>  m_w_gamma;				// used for each patch to find a perpendicular foot		
		std::vector<float>  m_w_kappa;				// used for each patch when summing up bone-binding energy

		std::vector<int>	m_patch_to_joint;		// ROOT_T included, but no one attached to it, i.e., no one shows "0"
		std::vector<int>	m_patch_to_joint_child;	// child of the attched joint for evert patch
		
		// skeleton adj
		std::vector<int>	m_jOffsp;				// start with the offspring of ROOT_R, ROOT_T not included
		std::vector<int>	m_jOffsp_bounds;		// length: 16 - 1 (no ROOT_T) + 1 (bounds), n*16 - n + 1
		std::vector<int>	m_jp;
		std::vector<int>	m_jp_bounds;

		//how many connected component per subject
		std::vector<int>	m_patch_comp;		
		int					m_numComponents;

		//LCF test		
		std::vector<float3> m_refVectors_LCF;

		// Joint output thing
		Eigen::MatrixXf		m_B;		
		Eigen::VectorXf		m_Delta;		

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

void build_dXN0Smooth( int           numPatches,
                       const float*  RT,
                       const int*    sadj,
                       const int*    sadj_bounds,
                       const int*    pv_bounds,
                       const float3* XN0,
                       float3*       dXN0_smooth );

void build_dXN0Smooth_from_dXN0( int           numPatches,
                                 const float*  RT,
                                 const int*    sadj,
                                 const int*    sadj_bounds,
                                 const int*    pv_bounds,
                                 const float3* dXN0,
                                 float3*       dXN0_smooth );

void build_dJX0_smooth(	const	std::vector<rigidT>& RT0,
						const	std::vector<float3>& JX0, 
						const	std::vector<int>&	jp, 
						const	std::vector<int>&	jp_bound, 
								std::vector<float3>& dJX0_smooth);

void build_w_JX0_smooth(	int numJoint,
								const	std::vector<float3>& dJX0_smooth, 							
								const	std::vector<int>&	jp_bound, 
										std::vector<float>&	w_JX_smooth);


void fillBMatrix(const std::vector<float>&			gamma, 
				 const std::vector<float>&			kappa,
				 const std::vector<int>&			pj,
				 const std::vector<int>&			pj_child,
				 Eigen::MatrixXf&	B);


void fillDeltaVec(		const std::vector<float>&	kappa,
						const std::vector<rigidT>&	RT, 
						const float3*				beta0,
						Eigen::VectorXf&			Delta);



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

void transform_Cloud_X(int           numPatches,
						const float*  RT,
						const int*    pv_bounds,
						const float3* dX0,
						float3*       X);

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

// ######################################################
// ######################################################
// Rigged Mesh stuff
// ######################################################
// ######################################################
int getChildId(int ji, 
			   const std::vector<KBJoint>& KBTree);

int getChildIdKNN(const float3& ci,  
				  const std::vector<KBJoint>& KBTree, 
				  const std::vector<float3>& j_coords);

inline void removeROOT_T(const std::vector<int>& before, std::vector<int>& after);

void compute_patch_to_joint( const SklRgedPatchedCloud::Patching&	patches,
   				             const RiggedMesh&				mesh,
									std::vector<int>&				patch_to_joint );

void compute_patch_to_joint_child(	int							numPatches,
									const rigidT*				RT0,
									const std::vector<KBJoint>& KBTree, 
									const std::vector<float3>&	j_coords, 
									const int*					p2j, 									
									std::vector<int>&			p2j_child);

void compute_jointOffspring(const	std::vector<KBJoint>& KBTree, 
									std::vector<int>& jOffsp, 
									std::vector<int>& jOffsp_bounds);

void find_idx_in_jOffsp(const	std::vector<int>& p2j_child, 
						const	std::vector<int>& jOffsp, 
								std::vector<int>& idx_in_jOff);

void build_joint_to_patch(	int numJoint, 
								const std::vector<int>&	pj, 
								const std::vector<int>&	pj_child, 
									  std::vector<int>&	jp, 
									  std::vector<int>& jp_bound);

void build_joint_coords0(const	std::vector<KBNode>&   KBNodes, 
								std::vector<float3>&		 joint_coords0);

void build_gamma(	int							numPatches,
					const rigidT*				RT0, 
					const int*					p2j,
					const int*					p2j_child,
					const std::vector<float3>&	j_coords, 
					std::vector<float>&			gamma);


void build_kappa(	const std::vector<float>& gamma, 
					const std::vector<int>& pj,
					bool weightflag,
					std::vector<float>& kappa);

void build_ppdicular_ft0(	const std::vector<float3>&	j_coords, 
							const int*					p2j, 
							const int*					p2j_child, 
							const std::vector<float>&	gamma, 
							std::vector<float3>&		ppdicular_ft0);

void build_beta0(	const rigidT*				RT0,
					const std::vector<float3>&	ppdicular_ft0, 
					std::vector<float3>&		beta0);

void rotate_beta0(	int           numPatches,
					const float*  RT,	
				    const float3* beta0,
			        float3*       beta );

void transform_ppdicular_ft0( int           numPatches,
							  const float*  RT,	
							  const float3* beta0,
							  float3*       ppdicular_ft );

void rotate_dJX0_smooth(	const std::vector<int>& jp,
							const std::vector<int>& jp_bound,
							const float*  RT,	
							const float3* dJX0_smooth,
							float3*       dJX_smooth );

void transform_JX0_smooth(		const std::vector<int>& jp,
								const std::vector<int>& jp_bound,
								const float*  RT,	
								const float3* dJX0_smooth,
								float3*       JX0_tfed );

// ###################################################
// ###################################################
// Some output helpers (for check the correctness)
// ###################################################
// ###################################################
void outputPatchCenters(int numPatches, const char* filename, const rigidT* RT);
void outputJoints(int numJoints, const char* filename, const std::vector<float3>&	j_coords);
void outputGamma(int numPatches, const char* filename, const std::vector<float>& gamma);
void outputPPfts(int numPatches, const char* filename, const std::vector<float3>&  ppdicular_ft0);
void outputBeta(int numPatches, const char* filename, const std::vector<float3>&  beta0);


} // end namespace PatchedCloud

#endif
