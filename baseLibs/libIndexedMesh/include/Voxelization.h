/* *************************************************
 * Copyright (2014) : Paul Huang
 * *************************************************/
#ifndef VOXELIZATION_H_DEFINED
#define VOXELIZATION_H_DEFINED
// std
#include <IndexedMesh/IndexedMesh3D.h>
#include <tuple>

// opencv
#include <cv.h>

enum LCF { SIGN, EVD, MY };
// faust
//#define CUBE_SIZE_LCF 12 //cube size for the local reference frame vector
//#define CUBE_SIZE_FEATURE 2 //cube size for the cube substraction feature
// others
#define CUBE_SIZE_LCF 7 //cube size for the local reference frame vector
#define CUBE_SIZE_FEATURE 2 //cube size for the cube substraction feature
//#define SignDistVolume

namespace IndexedMesh3D
{
	//using namespace cv;

	const float STATS_HUGE_NUMBER = 100000.0f;

	class INDEXEDMESH_API Voxelizer{
		public:

			typedef std::tuple<int, int, int, float, float3>  Neighbors;

			Voxelizer::Voxelizer(	float   Xmin,	float cubeSizeX, float   Xmax,	
									float   Ymin,	float cubeSizeY, float   Ymax,
									float   Zmin,	float cubeSizeZ, float   Zmax, int LCF_FLAG=0, bool SHOW_FLAG=false);

			Voxelizer::~Voxelizer();
			inline float3 getVoxelCoord(int voxelInd);
			// methods

			void voxelization(	const std::vector<float3>&		coords_s_demean, 
								const std::vector<Triangle>&	triangles_s,  
								const std::vector<float3>&		normals_v, 
								const std::vector<float3>&		normals_tri, 
									  std::vector<int>&			Vtx2Vxl, 
									  std::vector<float3>&		normals_Vol, 
									  std::vector<float3>&		sumCube, 
									  std::vector<float3>&		ref_frame_vectors, 
									  std::vector<float3>&		limits);

			void voxelization(	const std::vector<float3>&		coords_s_demean, 
								const std::vector<Triangle>&	triangles_s,  
								const std::vector<float3>&		normals_v, 
								const std::vector<float3>&		normals_tri, 
									  std::vector<int>&			Vtx2Vxl, 
									  std::vector<float3>&		normals_Vol, 
									  std::vector<float3>&		sumCube, 
									  std::vector<std::vector<float3>>& ref_frame_vectors, 
									  std::vector<float3>&		limits);

#ifdef SignDistVolume	

			void voxelization(  const std::vector<float3>&		coords_s_demean,
								const std::vector<Triangle>&	triangles_s,
								const std::vector<float3>&		normals_vtx,
								const std::vector<float3>&		normals_tri,
								std::vector<int>&				Vtx2Vxl,
								std::vector<float3>&			normals_Vol,
								std::vector<float>&				integralVol,
								std::vector<float3>&			ref_frame_vectors,
								std::vector<float3>& limits);
#endif
			// aux. methods
			int bwlabel(const IplImage* const img, int n, int* labels);
			inline void getCornerPointsCoords(const cv::Point3i& cubi, const float* const cubeSize, 
										float Xmin, float Ymin, float Zmin, 
										std::vector<float3>& corners_coords);
		private:

			float m_cubeSize[3];
			float m_Xmin, m_Xmax;            
			float m_Ymin, m_Ymax;
			float m_Zmin, m_Zmax;

			float m_max_length_LCF;
	
			unsigned int m_L, m_W, m_H;
			unsigned int m_WxL;
			unsigned int m_NUM_Vxls;

			int m_LCF_FLAG;

			bool m_SHOW_FLAG;

			std::vector<std::tuple<unsigned int, unsigned int, unsigned int> > m_ArrayIdx2VolCoords;
			std::vector<std::vector<unsigned int> >	m_cube_neighbors_idx_VolumeSub;
			std::vector<Neighbors>	m_cube_neighbors_idx_LCF;
	};
}

#endif