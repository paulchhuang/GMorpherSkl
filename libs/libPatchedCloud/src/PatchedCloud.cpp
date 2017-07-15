/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#include <PatchedCloud.h>
#include "Patching_func.h"
#include <cassert>
#include <CC3D/RigidMotion.h>


namespace SklRgedPatchedCloud
{

using CC3D::make_float3;
using CC3D::make_quaternion;
using CC3D::RigidTFlatten;

PatchedCloud::PatchedCloud( const Patching&				 patches,
                            double						 distVariance,
                            const std::vector<float3>&	 coords,
                            const std::vector<float3>&	 normals, 
							const std::vector<KBNode>&   KBNodes,
							const std::vector<KBJoint>&  KBJoints,
							const RiggedMesh&   mesh)
{
	int nVertices        = patches.numVertices();
	int nPatches         = patches.numPatches();
	int nVertices_smooth = patches.numSmoothVertices();

	assert( int(coords.size()) == nVertices );
	assert( int(normals.size()) == nVertices );

	// copy the patching stuff
	m_vp            = patches.vp();
	m_pv            = patches.pv();
	m_pv_bounds     = patches.pv_boundaries();
	m_sadj          = patches.sorted_patch_adj();
	m_sadj_bounds   = patches.sorted_patch_adj_boundaries();
	m_smooth_bounds = patches.smooth_boundaries();

	// resize the vects
	m_X0.resize(nVertices);
	m_N0.resize(nVertices);
	m_RT0.resize(nPatches);
	m_dX0_smooth.resize(nVertices_smooth);
	m_w_smooth.resize(nVertices_smooth);
	

	m_numComponents = SklRgedPatchedCloud::getConnectedComponents(m_sadj, m_sadj_bounds, m_patch_comp);
	// reinit the coords
	reinit( coords,normals, true, distVariance);
	sklinit( patches, KBNodes, KBJoints, mesh, true);
}

PatchedCloud::PatchedCloud( const Patching&				 patches,
                            double						 distVariance,
                            const std::vector<float3>&	 coords,
                            const std::vector<float3>&	 normals, 
							const std::vector<KBNode>&   KBNodes,
							const std::vector<KBJoint>&  KBJoints,
							const RiggedMesh&   mesh,		
							const std::vector<float3>&	 refVectors0)
{
	int nVertices        = patches.numVertices();
	int nPatches         = patches.numPatches();
	int nVertices_smooth = patches.numSmoothVertices();

	assert( int(coords.size()) == nVertices );
	assert( int(normals.size()) == nVertices );

	// copy the patching stuff
	m_vp            = patches.vp();
	m_pv            = patches.pv();
	m_pv_bounds     = patches.pv_boundaries();
	m_sadj          = patches.sorted_patch_adj();
	m_sadj_bounds   = patches.sorted_patch_adj_boundaries();
	m_smooth_bounds = patches.smooth_boundaries();

	m_numComponents = SklRgedPatchedCloud::getConnectedComponents(m_sadj, m_sadj_bounds, m_patch_comp);
	// resize the vects
	m_X0.resize(nVertices);
	m_N0.resize(nVertices);
	m_RT0.resize(nPatches);
	m_dX0_smooth.resize(nVertices_smooth);
	m_w_smooth.resize(nVertices_smooth);
	
	// reinit the coords
	reinit( coords,normals, true, distVariance);
	sklinit( patches, KBNodes, KBJoints, mesh, true);

	
	// LCF_test	
	m_refVectors_LCF.resize(nVertices);

	#pragma omp parallel for 
	for(int vi=0;vi<nVertices;++vi){ 		
		m_refVectors_LCF[vi]	= refVectors0[m_pv[vi]];
	}
}

PatchedCloud::PatchedCloud( const Patching&				 patches,
                            double						 distVariance,
                            const std::vector<float3>&	 coords,
                            const std::vector<float3>&	 normals)
{
	int nVertices        = patches.numVertices();
	int nPatches         = patches.numPatches();
	int nVertices_smooth = patches.numSmoothVertices();

	assert( int(coords.size()) == nVertices );
	assert( int(normals.size()) == nVertices );

	// copy the patching stuff
	m_vp            = patches.vp();
	m_pv            = patches.pv();
	m_pv_bounds     = patches.pv_boundaries();
	m_sadj          = patches.sorted_patch_adj();
	m_sadj_bounds   = patches.sorted_patch_adj_boundaries();
	m_smooth_bounds = patches.smooth_boundaries();

	// resize the vects
	m_X0.resize(nVertices);
	m_N0.resize(nVertices);
	m_RT0.resize(nPatches);
	m_dX0_smooth.resize(nVertices_smooth);
	m_w_smooth.resize(nVertices_smooth);
	
	m_numComponents = SklRgedPatchedCloud::getConnectedComponents(m_sadj, m_sadj_bounds, m_patch_comp);
	// reinit the coords
	reinit( coords,normals, true, distVariance);
						
}

void PatchedCloud::reinit( const std::vector<float3>& coords,
                           const std::vector<float3>& normals,
                           bool recomputeWeights, double distVariance)
{
	int nVertices        = numVertices();
	int nPatches         = numPatches();
	int nVertices_smooth = numVertices_smooth();

	assert( int(coords.size()) == nVertices );
	assert( int(normals.size()) == nVertices );

	// reindex vertices
	#pragma omp parallel for 
	for(int vi=0;vi<nVertices;++vi) {
		m_X0[vi] = coords[m_pv[vi]];
		m_N0[vi] = normals[m_pv[vi]];
	}

	// compute mean
	#pragma omp parallel for 
	for( int pi=0;pi<nPatches;++pi ) {
		int vbegin = m_pv_bounds[pi];
		int vend   = m_pv_bounds[pi+1];
		float3 ci = make_float3(0,0,0);
		for(int vi=vbegin;vi!=vend;++vi) ci += m_X0[vi];
		ci = ci / std::max(1, m_pv_bounds[pi+1]- m_pv_bounds[pi] );
		m_RT0[pi].m_q = make_quaternion(1,0,0,0);
		m_RT0[pi].m_t = ci;			// Paul: notice that the mean of the patch is not the coordinate of the center vertex.
	}

	// compute smooth vertices
	build_dX0_smooth( nPatches, &m_sadj[0], &m_sadj_bounds[0], &m_pv_bounds[0], &m_RT0[0], &m_X0[0], &m_dX0_smooth[0] );

	// recompute weights
	if( recomputeWeights ){
		assert( distVariance != 0. );
		build_w_smooth( nPatches, distVariance, &m_sadj_bounds[0], &m_pv_bounds[0], &m_smooth_bounds[0], &m_dX0_smooth[0], &m_w_smooth[0] );
	}
}


void PatchedCloud::sklinit(	const Patching&				patches,
							const std::vector<KBNode>&  KBNodes,
							const std::vector<KBJoint>&  KBJoints,
							const RiggedMesh&  mesh, 
							bool recomputeWeights){
	
	int nPatches         = numPatches();							
	//how many comps are there
	

	// start to compute gamma
	build_joint_coords0(KBNodes, m_joint_coords0);
	compute_patch_to_joint(patches, mesh, m_patch_to_joint);
	compute_patch_to_joint_child(nPatches, &m_RT0[0], KBJoints, m_joint_coords0, &m_patch_to_joint[0], m_patch_to_joint_child);
	compute_jointOffspring(KBJoints, m_jOffsp, m_jOffsp_bounds);


	build_gamma(nPatches, &m_RT0[0], &m_patch_to_joint[0], &m_patch_to_joint_child[0], m_joint_coords0, m_w_gamma);
	build_kappa(m_w_gamma, m_patch_to_joint, recomputeWeights, m_w_kappa);
	build_ppdicular_ft0(m_joint_coords0, &m_patch_to_joint[0],  &m_patch_to_joint_child[0], m_w_gamma, m_ppdicular_ft0);
	build_beta0(&m_RT0[0],m_ppdicular_ft0,m_beta0);

	// helpers
	/*outputJoints(m_joint_coords0.size(), "Joints_coords", m_joint_coords0);
	outputPatchCenters(nPatches, "Patches center", &m_RT0[0]);
	outputPPfts(nPatches, "PPft", m_ppdicular_ft0);
	outputGamma(nPatches, "Gamma", m_w_gamma);
	outputBeta(nPatches, "Beta", m_beta0);*/
	

	// skeleton blender

		// A. beta coordinate
	m_B.resize(3*nPatches, 3*(m_joint_coords0.size()-m_numComponents));
	fillBMatrix(m_w_gamma, m_w_kappa, m_patch_to_joint, m_patch_to_joint_child, m_B);
	m_Delta.resize(3*nPatches, 1);

		//B. patch-based blend
	build_joint_to_patch((m_joint_coords0.size()-m_numComponents), m_patch_to_joint, m_patch_to_joint_child,  m_jp, m_jp_bounds);
	build_dJX0_smooth(m_RT0, m_joint_coords0, m_jp, m_jp_bounds, m_dJX0_smooth);

	m_w_JX_smooth.resize(m_dJX0_smooth.size());
	build_w_JX0_smooth((m_joint_coords0.size()-m_numComponents),m_dJX0_smooth, m_jp_bounds, m_w_JX_smooth);
}

void PatchedCloud::blend( const std::vector<rigidT>& RT, std::vector<float3>& X) const
{
	int nPatches = m_pv_bounds.size()-1;
	float* RTf = new float[12*nPatches];
	RigidTFlatten(RT, RTf);

	X.resize( m_X0.size() );
	interpolate_cloud_smooth( nPatches,  RTf, &m_sadj[0], &m_sadj_bounds[0], &m_pv_bounds[0], &m_dX0_smooth[0], &m_w_smooth[0],  &X[0] );

	delete[] RTf;
}


void PatchedCloud::RTfromCloud( const std::vector<float3>& X, std::vector<rigidT>& RT) const
{
	int nPatches = m_pv_bounds.size()-1;
	RT.resize( nPatches );

	for(int pi=0;pi<nPatches;++pi)
	{
		std::vector<float3>::const_iterator        x_itr = X.begin()    + m_pv_bounds[pi];
		std::vector<float3>::const_iterator       x0_itr = m_X0.begin() + m_pv_bounds[pi];
		const std::vector<float3>::const_iterator x0_end = m_X0.begin() + m_pv_bounds[pi+1];


		CC3D::CovMatBuilder covMat;
		while( x0_itr != x0_end ) {
			covMat.pushCorrespondance(x0_itr->x, x0_itr->y, x0_itr->z, x_itr->x, x_itr->y, x_itr->z, 1.0); // accum
			x0_itr++;
			x_itr++;
		}
		double cov[9], c1[3], c2[3], q[4];
		covMat.getCovMat( cov, c1, c2 );
		CC3D::RigidMotionHorn Horn;
		int INFO = Horn.getQuaternion( cov, q[0], q[1], q[2], q[3]);
		RT[pi].m_q = make_quaternion( q[0], q[1], q[2], q[3] );
		float R[9];
		CC3D::Quat_to_Rot<float>((const float*)&(RT[pi].m_q), R);
		RT[pi].m_t = make_float3( c2[0], c2[1], c2[2] );
	}
}

void PatchedCloud::genJointBasedColor( std::vector<IndexedMesh3D::ColorVec>& colors ) const{
		// resize output vector
		colors.resize( numVertices() );

		// generate random colors
		std::vector<IndexedMesh3D::ColorVec> patchcolors( numPatches() );
		std::vector<IndexedMesh3D::ColorVec> jointcolors( numJoints() );

		for(int ji = 0; ji < numJoints(); ++ji){
			jointcolors[ji].r  = ( (rand() %1000)/1000.0);
			jointcolors[ji].g  = ( (rand() %1000)/1000.0);
			jointcolors[ji].b  = ( (rand() %1000)/1000.0);	
		}

		for(size_t pi=0; pi< numPatches(); ++pi ){
			if(m_w_gamma[pi] < 0.05){	// patch near the child of the attached joint
				patchcolors[pi].r  = 1;
				patchcolors[pi].g  = 1;
				patchcolors[pi].b  = 1;				
			}
			else if(m_w_gamma[pi] > 0.95){	// patch near attached joint
				patchcolors[pi].r  = 0;
				patchcolors[pi].g  = 0;
				patchcolors[pi].b  = 0;	
			}
			else{
				float perturbation = (rand() %100)/1000.0;
				patchcolors[pi].r  = jointcolors[m_patch_to_joint[pi]].r + perturbation;
				patchcolors[pi].g  = jointcolors[m_patch_to_joint[pi]].g + perturbation;
				patchcolors[pi].b  = jointcolors[m_patch_to_joint[pi]].b + perturbation;				
			}
			patchcolors[pi].a = 1.0;
		}

		// fill the vertex colors
		for(size_t vi=0; vi<numVertices(); ++vi) {
			int pid    = m_vp[vi];
			colors[vi] = patchcolors[pid];
		}
	}

	void PatchedCloud::JXfromRT_Bone(const std::vector<rigidT>& RT, std::vector<float3>& JX){
		
		fillDeltaVec(m_w_kappa, RT, &m_beta0[0], m_Delta);

		int nJoints = JX.size();
		
		Eigen::VectorXf J = m_B.jacobiSvd(ComputeThinU | ComputeThinV).solve(m_Delta);

		assert(nJoints==(J.size()/3));

		for (int ji = 0; ji < nJoints; ++ji){
			JX[ji].x = J(3*ji+0,0);
			JX[ji].y = J(3*ji+1,0);
			JX[ji].z = J(3*ji+2,0);
		}

	}


	void PatchedCloud::JXfromRT_Joint(const std::vector<rigidT>& RT, std::vector<float3>& JX){
		
		int numJoints = JX.size();				// here numJoints: 15, 30, 45, etc

		float* RTf = new float[12*RT.size()];
		RigidTFlatten(RT, RTf);

		for (int ji = 0; ji < numJoints; ++ji){
			int begin = m_jp_bounds[ji];
			int end = m_jp_bounds[ji+1];
			int idx = begin;

			JX[ji] = make_float3(0.0,0.0,0.0);

			std::vector<float3>::iterator		dJXi_smooth_itr		= m_dJX0_smooth.begin() + m_jp_bounds[ji];
			std::vector<float3>::const_iterator dJXi_smooth_end		= m_dJX0_smooth.begin() + m_jp_bounds[ji+1];
			
			while (dJXi_smooth_itr!=dJXi_smooth_end)
			{
				const float* Rj  = RTf + 12*m_jp[idx];
				const float3& tj = *(const float3*)(Rj+9);
				JX[ji] += m_w_JX_smooth[idx]*(prod3(Rj, *dJXi_smooth_itr++) + tj);			
				
				idx++;
			}
			assert(idx==end);
		}

		delete[] RTf;
		//
	}

} // end namespace PatchedCloud


