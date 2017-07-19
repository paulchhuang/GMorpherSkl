/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#include <OptionParser.h>
// libIndexedMesh
#include <IndexedMesh/IndexedMesh3D.h>
// libGMorpher 
#include <Solver_ceres.h>
// libPatchedCloud_GMorpher
#include <EnergyTerm_Rigidity_CEDRIC_ceres.h>
// local shit
#include <EnergyTerm_SilFit.h>
// std
#include <list>
#include <stdlib.h> // for the getEnv

using namespace IndexedMesh3D;
using SklRgedPatchedCloud::EnergyTerm_Rigidity_CEDRIC_ceres;
using GMorpher::rigidT;

int main(int argc, char** argv)
{
	OptionParser opt(argc, argv);
	std::string meshBasename    = opt.getOption<std::string>("-meshBasename");
	std::string camBasename     = opt.getOption<std::string>("-camBasename");
	std::string silBasename     = opt.getOption<std::string>("-silBasename");
	std::string camIndex        = opt.getOption<std::string>("-camIndex");
	std::string outBasename     = opt.getOption<std::string>("-outBasename");
	int firstCam                = opt.getOption<int>("-firstCam");
	int lastCam                 = opt.getOption<int>("-lastCam");
	int R                       = opt.getOption<int>("-R"); // is the ref file name
	int F                       = opt.getOption<int>("-F");
	int L                       = opt.getOption<int>("-L");
	int patchSize               = opt.getOption<int>("-S");
	double SForce               = opt.getOption<double>("-SF");



	// ----------------------------------------------
	// ----------------------------------------------
	// A - We load the root mesh to have an idea of the dimensions we will be
	// Loading the mesh
	std::vector<Triangle>                triangles_0;
	std::vector<float3>                  coords_0;
	std::vector<float3>                  normals_0;
	std::vector<ColorVec>                colors_0;
	loadIndexedMesh3D(buildFilename(meshBasename, R), coords_0,
	                                                     colors_0, triangles_0);	

	averageNormals( triangles_0, coords_0, normals_0);
	
	// ----------------------------------------------
	// ----------------------------------------------
	// B - create the patching
		// creating adj list 
	std::vector<IndexedMesh3D::Edge>     CWEdges_0;
	std::vector<int>                     CWEdgesBoundaries_0;
	convertTrianglesToEdgeGraph(coords_0.size(), triangles_0, 
	                                           CWEdges_0,  CWEdgesBoundaries_0);
		// the patching 
	SklRgedPatchedCloud::Patching patches;
	patches.init( patchSize, CWEdges_0, CWEdgesBoundaries_0);
	std::vector<ColorVec> colors;
	patches.genColors(colors);
		// the mean Edge 
	double meanEdge_0 = getMeanEdgeLength(triangles_0, coords_0 );
		// the patchedCloud
	double distVar = meanEdge_0*meanEdge_0*patchSize*patchSize;
	SklRgedPatchedCloud::PatchedCloud PC( patches, distVar, coords_0, normals_0 );


	std::vector<Triangle>                triangles_r(triangles_0.size());
	
	for(int ti = 0;ti < triangles_r.size();++ti)
	{
		triangles_r[ti].v0 = std::find(patches.pv().begin(),patches.pv().end(),triangles_0[ti].v0)-patches.pv().begin();
		triangles_r[ti].v1 = std::find(patches.pv().begin(),patches.pv().end(),triangles_0[ti].v1)-patches.pv().begin();
		triangles_r[ti].v2 = std::find(patches.pv().begin(),patches.pv().end(),triangles_0[ti].v2)-patches.pv().begin();
	}
		


	// ----------------------------------------------
	// ----------------------------------------------
	// C - create the solver
	GMorpher::SolverCeres solver(PC.RT0(), std::vector<float3>());
	

	// ----------------------------------------------
	// ----------------------------------------------
	// D - Create the silhouette and rigidity energies
		// silhouette
	std::list<int> ids; 
	for(int i=firstCam; i<=lastCam; ++i) ids.push_back(i);
	
	boost::shared_ptr<EnergyTerm_SilFit> ESil(new EnergyTerm_SilFit(
	                                   "myWindow", SForce, patches, PC, triangles_r));
	ESil->setCameras(ids, camBasename, buildFilename(silBasename, camIndex, F), 0.001, 1000);	
	
		// rigidity
	boost::shared_ptr<EnergyTerm_Rigidity_CEDRIC_ceres> EReg(new EnergyTerm_Rigidity_CEDRIC_ceres(PC));
		// the list of energies 
	std::list<GMorpher::EPtrCeres> energies;
	energies.push_back( ESil );
	energies.push_back( EReg );

	// ----------------------------------------------
	// ----------------------------------------------
	// E - process the sequence
	std::vector<float3> pv_coords;
	std::vector<rigidT> RT = PC.RT0();
	std::vector<Triangle>                triangles_t;
	std::vector<float3>                  coords_t;
	std::vector<float3>                  normals_t;
	std::vector<ColorVec>                colors_t;
	for(int t=F; t<=L; ++t)
	{
		
		loadIndexedMesh3D(buildFilename(meshBasename, t).c_str(), coords_t, 
		                                                 colors_t, triangles_t);
			
		
		//PC.coords_to_pvcoords( coords_t, pv_coords );		
		pv_coords.resize(coords_t.size());

		#pragma omp parallel for
		for(int vi=0;vi<PC.numVertices();++vi) 
			pv_coords[vi] = coords_t[PC.pv()[vi]];								
		
		
		// set the RT that were obtained previously
		PC.RTfromCloud( pv_coords, RT);
		
		// set the silhouettes
		
		ESil->setImages( buildFilename(silBasename, camIndex, t) );
		
		// solve
		solver.solve( 20, // maxIter //		              10, // maxSubdiv 
		              energies,
					  RT, 
					  std::vector<float3>());
		
		// save the result
		PC.blend( RT, pv_coords);
		
		#pragma omp parallel for
		for(int vi=0;vi<PC.numVertices();++vi) 
			coords_t[PC.pv()[vi]] = pv_coords[vi];						
		

		//std::cout<<"before save"<<std::endl;
		saveIndexedMesh3D(buildFilename(outBasename, t).c_str(), coords_t, 
		                                                   colors, triangles_0);
		//std::cout<<"after save"<<std::endl;
	}

	return 0;
}
