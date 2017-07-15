/* *************************************************
 * Copyright (2016) : Paul Huang
 * *************************************************/
#include <OptionParser.h>
// std
#include <iostream>
#include <fstream>	
// libPatchedMesh
#include <IndexedMesh/IndexedMesh3D.h>
#include <Patching.h>
#include <RiggedMesh.h>

// GMorpher
#include <EnergyTerm_Rigidity_CEDRIC_ceres.h>
#include <EnergyTerm_BoneBinding_ceres.h>
#include <EnergyTerm_3DConstraint_ceres.h>
#include <Solver_EM_ceres.h>	


#define NUM_JOINTs_SUB 16


using namespace IndexedMesh3D;

inline void convertXNToIndexedMesh(const std::vector<int>&				 jv,
	const std::vector<int>&                 pv,
	const std::vector<float3>&    X,
	std::vector<float3>&                    coords);

int main(int argc, char** argv)
{
	OptionParser opt(argc, argv);
	bool        probabilistic       = opt.getOption<bool>("-probabilistic");
	bool        prediction          = opt.getOption<bool>("-prediction");
	std::string meshRef             = opt.getOption<std::string>("-meshRef");
	std::string meshBaseName        = opt.getOption<std::string>("-meshBaseName");
	std::string outBaseName         = opt.getOption<std::string>("-outBaseName");
	std::string outJointName        = opt.getOption<std::string>("-outJointName");
	std::string SKL_FileName        = opt.getOption<std::string>("-skl");	
	std::string VJ_FileName         = opt.getOption<std::string>("-vj");
	int         firstFrame          = opt.getOption<int>("-F");
	int         lastFrame           = opt.getOption<int>("-L");
	int			patchSize           = opt.getOption<int>("-S");
	float		EMForce_ICP         = opt.getOption<float>("-EMF_ICP");
	float		ESForce             = opt.getOption<float>("-ESF");
	float		normThresh          = opt.getOption<float>("-NThresh");
	float		sigma0              = opt.getOption<float>("-sigma0");
	float		Eoutlier            = opt.getOption<float>("-Eoutlier");

	// ###################################################
	// 1- load the skl file
	std::vector<KBNode>              KBNodes;
	std::vector<std::string>         KBJointNames;
	std::vector<KBJoint>			 KBJoints;
	std::vector<double>              KBParamsVec;
	std::vector<double>              KBStateVec;
	// load the skel

	KBParse(SKL_FileName.c_str(), KBNodes);
	KBGenJoints(KBNodes, KBJointNames, KBJoints, KBParamsVec, KBStateVec);

	std::cout << "\nloading tracking Ref. surface: ";
	SklRgedPatchedCloud::RiggedMesh mesh(KBJoints, KBJointNames, meshRef.c_str(), VJ_FileName.c_str(), KBJoints.size() / NUM_JOINTs_SUB);

	// ###################################################
	// 2- load the start mesh
	std::vector<IndexedMesh3D::Edge>     CWEdges_ref;
	std::vector<int>                     CWEdgesBoundaries_ref;
	std::vector<Triangle>                triangles_ref;
	std::vector<float3>                  coords_ref;
	std::vector<float3>                  normals_ref;
	std::vector<float3>					 normals_tri_ref;
	std::vector<ColorVec>                colors_ref;

	triangles_ref = mesh.triangles();												// the mesh representation: coord_ref, triangles_ref
	coords_ref = mesh.coords();


	convertTrianglesToEdgeGraph(coords_ref.size(), triangles_ref,
		CWEdges_ref, CWEdgesBoundaries_ref);	// the graph representation of the mesh: CWEdges_ref, CWEdgesBoundaries_ref
	IndexedMesh3D::averageNormals(triangles_ref, coords_ref, normals_ref, normals_tri_ref);
	float  meanEdge_ref = getMeanEdgeLength(triangles_ref, coords_ref);

	// ###################################################
	// A - We create the Patchings
	// Patching 
	SklRgedPatchedCloud::Patching patches;
	patches.init(patchSize, CWEdges_ref, CWEdgesBoundaries_ref);		// patch the 
	std::vector<ColorVec> colors;
	patches.genColors(colors);
	std::vector<int> vp_reidx = patches.vp();

	colors_ref.resize(colors.size());
	double distVar = meanEdge_ref*meanEdge_ref*patchSize*patchSize;

#pragma omp parallel for
	for (int vi = 0; vi < mesh.coords().size(); ++vi)
	{
		colors_ref[mesh.jv()[vi]] = colors[vi];
		vp_reidx[mesh.jv()[vi]] = patches.vp()[vi];
	}

	// ###################################################
	// B - We create the solver
	// point cloud

	SklRgedPatchedCloud::PatchedCloud PC(patches, distVar, coords_ref, normals_ref, KBNodes, KBJoints, mesh);
	GMorpher::EPtrCeres EBoneCeres = GMorpher::EPtrCeres(new SklRgedPatchedCloud::EnergyTerm_BoneBindingCeres(PC, ESForce));
	GMorpher::EPtrCeres ERigCeres = GMorpher::EPtrCeres(new SklRgedPatchedCloud::EnergyTerm_Rigidity_CEDRIC_ceres(PC));

	Solver_EMCeres SolverCeres_EM(prediction, PC);// , PC.sadj(), PC.sadj_bounds(), PC.jOffsp(), PC.jOffsp_bounds());


	std::vector<GMorpher::rigidT>	RT = PC.RT0();
	std::vector<float3>		JX((KBJoints.size() / NUM_JOINTs_SUB)*(NUM_JOINTs_SUB - 1));
	int count = 0;
	for (int ji = 0; ji < PC.joint_coords0().size(); ++ji)
	{
		if ((ji%NUM_JOINTs_SUB) != 0)
			JX[count++] = PC.joint_coords0()[ji];
	}

	std::vector<float3>	Xtemp, coords_temp;

	for (size_t t = firstFrame; t <= lastFrame; t++)
	{
		std::vector<Triangle>                triangles_t;
		std::vector<float3>                  coords_t;
		std::vector<ColorVec>                colors_t;
		std::vector<float3>                  normals_t;
		
		loadIndexedMesh3D(buildFilename(meshBaseName, t), coords_t, colors_t, triangles_t);
		averageNormals(triangles_t, coords_t, normals_t);					//triangles_dst is not used


		int numTargetPoints = coords_t.size();

	
		std::list<GMorpher::EPtrCeres>  energiesCeres;


		clock_t start, stop;
		
		float sigma = sigma0*meanEdge_ref;

			
		energiesCeres.push_back(EBoneCeres);
		energiesCeres.push_back(ERigCeres);


		start = clock();
		int niterationsEM = SolverCeres_EM.solve(energiesCeres,											
											probabilistic,
											30, //maxIter_EM,
											1, //maxIter_M
											10, //maxSubDiv,
											meanEdge_ref,
											normThresh,
											Eoutlier,
											EMForce_ICP,
											coords_t,
											normals_t,
											RT,
											JX,
											sigma);

		stop = clock();
		std::cout << "spend: " << (stop - start) << " ms, " << niterationsEM << " EM interations on solving." << std::endl;
		

		
		PC.blend(RT, Xtemp);
		convertXNToIndexedMesh(mesh.jv(), PC.pv(), Xtemp, coords_temp);

		
		std::vector<float3> jcoords_temp(JX.begin(), JX.end());		
		SklRgedPatchedCloud::outputJoints(jcoords_temp.size(), buildFilename(outJointName, t).c_str(), jcoords_temp);

		saveIndexedMesh3D((buildFilename(outBaseName, t) + ".off").c_str(), coords_temp, colors_ref, mesh.triangles0());
	}	

	return 0;
}

void convertXNToIndexedMesh(const std::vector<int>&				 jv,
	const std::vector<int>&                 pv,
	const std::vector<float3>&    X,
	std::vector<float3>&                    coords)
{
	int numVertices = X.size();
	coords.resize(numVertices);

	std::vector<float3>		coord_tmp(numVertices);

#pragma omp parallel for
	for (int vi = 0; vi<numVertices; ++vi)
		coord_tmp[pv[vi]] = X[vi];

#pragma omp parallel for
	for (int vi = 0; vi < numVertices; ++vi)
		coords[jv[vi]] = coord_tmp[vi];

}
