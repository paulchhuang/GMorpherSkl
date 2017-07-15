/* *************************************************
 * Copyright (2014) : Paul Huang
 * *************************************************/
#pragma warning(disable: 4018) 
// std
#include <fstream>
#include <boost/foreach.hpp>
#include <Voxelization.h>
// opencv
#include <opencv2/highgui/highgui.hpp>

// macro
#define NO_OBJECT       0  
#define MIN(x, y)       (((x) < (y)) ? (x) : (y))  
#define ELEM(img, r, c) (CV_IMAGE_ELEM(img, unsigned char, r, c))  
#define ONETWO(L, r, c, col) (L[(r) * (col) + c])  
#define NUM_LOCALMAX 5

inline int find( int set[], int x ){  
    int r = x;  
    while ( set[r] != r )  
        r = set[r];  
    return r;  
}  

namespace IndexedMesh3D{
	Voxelizer::Voxelizer(	float   Xmin,	float cubeSizeX, float   Xmax,	
							float   Ymin,	float cubeSizeY, float   Ymax,
							float   Zmin,	float cubeSizeZ, float   Zmax, int LCF_FLAG, bool SHOW_FLAG) : 
							m_Xmin(Xmin), m_Ymin(Ymin), m_Zmin(Zmin), m_LCF_FLAG(LCF_FLAG),
							m_Xmax(Xmax), m_Ymax(Ymax), m_Zmax(Zmax), m_SHOW_FLAG(SHOW_FLAG){

		// set the volume.
		m_cubeSize[0] = cubeSizeX;
		m_cubeSize[1] = cubeSizeY;
		m_cubeSize[2] = cubeSizeZ;
		
		m_L			= (int)((m_Xmax-m_Xmin)/m_cubeSize[0]);
		m_W			= (int)((m_Ymax-m_Ymin)/m_cubeSize[1]);
		m_H			= (int)((m_Zmax-m_Zmin)/m_cubeSize[2]);
		m_WxL		= m_W*m_L;
		m_NUM_Vxls	= m_WxL*m_H;	


		std::cout<<"initializing a 3D volume: "<< m_L << "x" << m_W << "x" << m_H << std::endl;
		m_ArrayIdx2VolCoords.resize(m_NUM_Vxls);
		//m_cube_neighbors_idx_VolumeSub.resize(m_NUM_Vxls);
		//
		// pre-compute and store all the idx for neighboring voxels.
#pragma omp parallel for 
		for (int vli = 0; vli < m_NUM_Vxls; ++vli)
		{
			int curr_z = vli/m_WxL;
			int curr_x = (vli % m_WxL) / m_W;
			int curr_y = (vli % m_WxL) % m_W;
			m_ArrayIdx2VolCoords[vli] = std::make_tuple(curr_x, curr_y, curr_z);
		}
		
		m_max_length_LCF = sqrtf(3.0*CUBE_SIZE_LCF*CUBE_SIZE_LCF);


		// the support of different LCF method
		if(LCF_FLAG<=LCF::SIGN){			// signDist: visting voxels are border voxels on the cuboids
			{							
				int z = -CUBE_SIZE_LCF;
				do{
					for(int x=-CUBE_SIZE_LCF; x<=CUBE_SIZE_LCF; x++) {
						for(int y=-CUBE_SIZE_LCF; y<=CUBE_SIZE_LCF; y++) {						
							float3 local_coord_vec			= make_float3(x, y, z);
							float3 local_coord_vec_normed	= make_float3(x, y, z);		normalize(local_coord_vec_normed);

							float normalized_l = norm2(local_coord_vec)/m_max_length_LCF;
							m_cube_neighbors_idx_LCF.push_back(std::make_tuple(x, y, z, normalized_l, local_coord_vec_normed));
						}
					}

					z += 2*CUBE_SIZE_LCF;
				}while(z<=CUBE_SIZE_LCF);		
			}

			{
				int x = -CUBE_SIZE_LCF;
				do{
					for(int z=-CUBE_SIZE_LCF; z<=CUBE_SIZE_LCF; z++) {
						for(int y=-CUBE_SIZE_LCF; y<=CUBE_SIZE_LCF; y++) {						
							float3 local_coord_vec			= make_float3(x, y, z);
							float3 local_coord_vec_normed	= make_float3(x, y, z);		normalize(local_coord_vec_normed);

							float normalized_l = norm2(local_coord_vec)/m_max_length_LCF;
							m_cube_neighbors_idx_LCF.push_back(std::make_tuple(x, y, z, normalized_l, local_coord_vec_normed));
						}
					}

					x += 2*CUBE_SIZE_LCF;
				}while(x<=CUBE_SIZE_LCF);		
			}

			{
				int y = -CUBE_SIZE_LCF;
				do{
					for(int z=-CUBE_SIZE_LCF; z<=CUBE_SIZE_LCF; z++) {
						for(int x=-CUBE_SIZE_LCF; x<=CUBE_SIZE_LCF; x++) {						
							float3 local_coord_vec			= make_float3(x, y, z);
							float3 local_coord_vec_normed	= make_float3(x, y, z);		normalize(local_coord_vec_normed);

							float normalized_l = norm2(local_coord_vec)/m_max_length_LCF;
							m_cube_neighbors_idx_LCF.push_back(std::make_tuple(x, y, z, normalized_l, local_coord_vec_normed));
						}
					}

					y += 2*CUBE_SIZE_LCF;
				}while(y<=CUBE_SIZE_LCF);		
			}
		}
		else{	// MY, EVD: visting voxels are all voxels in the cuboids
			for(int z=-CUBE_SIZE_LCF; z<=CUBE_SIZE_LCF; z++) {
				for(int x=-CUBE_SIZE_LCF; x<=CUBE_SIZE_LCF; x++) {
					for(int y=-CUBE_SIZE_LCF; y<=CUBE_SIZE_LCF; y++) {						
						float3 local_coord_vec			= make_float3(x, y, z);
						float3 local_coord_vec_normed	= make_float3(x, y, z);
						normalize(local_coord_vec_normed);
						float normalized_l = norm2(local_coord_vec)/m_max_length_LCF;
						m_cube_neighbors_idx_LCF.push_back(std::make_tuple(x, y, z, normalized_l, local_coord_vec_normed));
					}
				}
			}
		}

		/*std::cout<< m_cube_neighbors_idx_LCF.size() << std::endl;
		system("pause");*/
		/*std::cout<< numNeighbor_LCF << std::endl;*/
	}

Voxelizer::~Voxelizer(){}

#ifdef SignDistVolume	
void Voxelizer::voxelization(	const std::vector<float3>& coords_s_demean, const std::vector<Triangle>&	triangles_s,
								const std::vector<float3>& normals_vtx, const std::vector<float3>&		normals_tri,
								std::vector<int>& Vtx2Vxl, std::vector<float3>& normals_Vol, std::vector<float>& integralVol, std::vector<float3>& ref_frame_vectors,
								std::vector<float3>& limits){

	size_t numVertex = coords_s_demean.size();

	// ##############################################################
	// 1. create vertice-voxels mapping list
#pragma omp parallel for 
	for (int vi = 0; vi < numVertex; ++vi){
		int idxX = floor((coords_s_demean[vi].x - m_Xmin) / m_cubeSize[0]);		// start counting at 0.
		int idxY = floor((coords_s_demean[vi].y - m_Ymin) / m_cubeSize[1]);
		int idxZ = floor((coords_s_demean[vi].z - m_Zmin) / m_cubeSize[2]);

		assert(idxX < m_L);
		assert(idxY < m_W);
		assert(idxZ < m_H);
		if (idxX < limits[0].x) {
			limits[0].x = idxX;
		}
		if (idxX > limits[1].x) {
			limits[1].x = idxX;
		}
		if (idxY < limits[0].y) {
			limits[0].y = idxY;
		}
		if (idxY > limits[1].y) {
			limits[1].y = idxY;
		}
		if (idxZ < limits[0].z) {
			limits[0].z = idxZ;
		}
		if (idxZ > limits[1].z) {
			limits[1].z = idxZ;
		}

		Vtx2Vxl[vi] = idxZ*m_WxL + idxX*m_W + idxY;
	}


	// ##############################################################
	// 2. create shell by checking triangle intersection
	std::vector<float3> axes(13);
	axes[0] = make_float3(1.0f, 0.0f, 0.0f);		axes[1] = make_float3(0.0f, 1.0f, 0.0f);	axes[2] = make_float3(0.0f, 0.0f, 1.0f);

	std::vector<float3>::const_iterator tnl_itr = normals_tri.begin();
	cv::Mat AABC(2, 3, CV_32SC1);

	std::vector<cv::Mat>	occ_slices(m_H);
#pragma omp parallel for 
	for (int zi = 0; zi < m_H; ++zi){
		occ_slices[zi].create(m_L, m_W, CV_8UC1);
		occ_slices[zi].setTo(255);
	}

	for (std::vector<Triangle>::const_iterator t_itr = triangles_s.begin(); t_itr != triangles_s.end(); ++t_itr){


		axes[3] = *tnl_itr++;

		std::vector<float3>	edges(3);
		edges.push_back(coords_s_demean[t_itr->v1] - coords_s_demean[t_itr->v0]);	// v1_v0		
		edges.push_back(coords_s_demean[t_itr->v2] - coords_s_demean[t_itr->v0]);	// v2_v0
		edges.push_back(coords_s_demean[t_itr->v2] - coords_s_demean[t_itr->v1]);	// v2_v1

		for (int ai = 0; ai < 9; ++ai)		axes[ai + 4] = CC3D::cross(edges[ai / 3], axes[ai % 3]);


		// ##############################################################
		// 2.1 identify the cubes which intersect with the bounding box of triangle
		AABC.at<int>(0, 0) = floor((std::min(std::min(coords_s_demean[t_itr->v1].x, coords_s_demean[t_itr->v0].x), coords_s_demean[t_itr->v2].x) - m_Xmin) / m_cubeSize[0]);
		AABC.at<int>(1, 0) = floor((std::max(std::max(coords_s_demean[t_itr->v1].x, coords_s_demean[t_itr->v0].x), coords_s_demean[t_itr->v2].x) - m_Xmin) / m_cubeSize[0]);
		AABC.at<int>(0, 1) = floor((std::min(std::min(coords_s_demean[t_itr->v1].y, coords_s_demean[t_itr->v0].y), coords_s_demean[t_itr->v2].y) - m_Ymin) / m_cubeSize[1]);
		AABC.at<int>(1, 1) = floor((std::max(std::max(coords_s_demean[t_itr->v1].y, coords_s_demean[t_itr->v0].y), coords_s_demean[t_itr->v2].y) - m_Ymin) / m_cubeSize[1]);
		AABC.at<int>(0, 2) = floor((std::min(std::min(coords_s_demean[t_itr->v1].z, coords_s_demean[t_itr->v0].z), coords_s_demean[t_itr->v2].z) - m_Zmin) / m_cubeSize[2]);
		AABC.at<int>(1, 2) = floor((std::max(std::max(coords_s_demean[t_itr->v1].z, coords_s_demean[t_itr->v0].z), coords_s_demean[t_itr->v2].z) - m_Zmin) / m_cubeSize[2]);



		std::vector<cv::Point3i> AABB_intersec_cubes;

		for (int xi = AABC.at<int>(0, 0); xi <= AABC.at<int>(1, 0); ++xi){
			for (int yi = AABC.at<int>(0, 1); yi <= AABC.at<int>(1, 1); ++yi){
				for (int zi = AABC.at<int>(0, 2); zi <= AABC.at<int>(1, 2); ++zi){
					cv::Point3i tmp(xi, yi, zi);
					AABB_intersec_cubes.push_back(tmp);
				}
			}
		}

		std::vector<float3>	corners_coords(8);

		// ##############################################################
		// 2.2 for all those cubes, perform SAT (Separating Axis Theorem) test 
		BOOST_FOREACH(const cv::Point3i& AABB_intersec_cubes_i, AABB_intersec_cubes){
			getCornerPointsCoords(AABB_intersec_cubes_i, m_cubeSize, m_Xmin, m_Ymin, m_Zmin, corners_coords);

			float	proj_tri[3];		float	proj_cub[8];
			int count = 0;
			BOOST_FOREACH(const float3& axis_i, axes){
				proj_tri[0] = dot(axis_i, coords_s_demean[t_itr->v0]);
				proj_tri[1] = dot(axis_i, coords_s_demean[t_itr->v1]);
				proj_tri[2] = dot(axis_i, coords_s_demean[t_itr->v2]);

				//#pragma omp parallel for		// don't do parallel computing here. it doesn't help at all, even makes it (much) slower 
				for (int cri = 0; cri < 8; ++cri)
					proj_cub[cri] = dot(axis_i, corners_coords[cri]);

				auto MinMaxTri = std::minmax_element(proj_tri, proj_tri + 3);
				auto MinMaxCub = std::minmax_element(proj_cub, proj_cub + 8);

				if ((MinMaxCub.first < MinMaxTri.second) || (MinMaxTri.first < MinMaxCub.second)){
					if (AABB_intersec_cubes_i.z*m_WxL + AABB_intersec_cubes_i.x*m_W + AABB_intersec_cubes_i.y < normals_Vol.size()){
						normals_Vol[AABB_intersec_cubes_i.z*m_WxL + AABB_intersec_cubes_i.x*m_W + AABB_intersec_cubes_i.y] += axes[3];
						occ_slices[AABB_intersec_cubes_i.z].at<unsigned char>(AABB_intersec_cubes_i.x, AABB_intersec_cubes_i.y) = 0;
						break;
					}
				}

			}
		}
	}

	assert(tnl_itr == normals_tri.end());

	// ##############################################################
	// 3. normalization of normals
#pragma omp parallel for 
	for (int ni = 0; ni < m_NUM_Vxls; ++ni){
		float norm = norm2(normals_Vol[ni]);
		if (norm != 0)	normals_Vol[ni] /= norm;
	}

	if (m_SHOW_FLAG){
		cv::namedWindow("A", CV_WINDOW_AUTOSIZE);
		cvNamedWindow("A", 0);
		for (int zi = 0; zi < m_H; ++zi){
			cv::imshow("A", occ_slices[zi]);
			cv::waitKey(30);
		}
	}

	// ##############################################################
	// 4. distinguish between inside (-2) and outside (2), by check the connected compoment (a bit ad-hoc)
#pragma omp parallel for 
	for (int zi = 0; zi < m_H; ++zi){
		IplImage imgtmp = IplImage(occ_slices[zi]);
		int* labels = new int[m_WxL];

		int conncomp = bwlabel(&imgtmp, 4, labels);

		int zixWxL = zi*m_WxL;
		int count = 0;
		const int* const begin = labels;
		const int* const end = labels + m_WxL;

		for (const int* label = begin; label != end; ++label){
			if (*label == conncomp){ normals_Vol[zixWxL + count] = make_float3(2.0f, 2.0f, 2.0f); }
			else if ((*label < conncomp) && (*label > 0)){ normals_Vol[zixWxL + count] = make_float3(-2.0f, -2.0f, -2.0f); }
			else{}
			count++;
		}

		/*for (int ri = 0; ri < L; ++ri){
		for (int ci = 0; ci < W; ++ci){
		int label = labels[ri*W+ci];
		if (label==conncomp){ normals_Vol[zixWxL + ri*W + ci] = make_float3(2.0f, 2.0f, 2.0f);}
		else if((label < conncomp)&&(label > 0)){ normals_Vol[zixWxL + ri*W + ci] =  make_float3(-2.0f, -2.0f, -2.0f);}
		else{}
		}
		}*/

		delete[] labels;
	}

	/*std::ofstream f_out("normalx_afterConnComp");
	for (int ri = 0; ri < m_L; ++ri){
	for (int ci = 0; ci < m_W; ++ci){
	f_out << normals_Vol[94*m_WxL + ri*m_L + ci].x << " ";
	}
	f_out << "\n";
	}
	f_out.close();*/

	//		// ##############################################################
	//		// 5. when needed, compute cubeSum (gotta replace this part with integral volume)
	//#pragma omp parallel for 
	//		for (int ni = 0; ni < m_NUM_Vxls; ++ni){
	//			
	//			int curr_x	   = std::get<0>(m_ArrayIdx2VolCoords[ni]);
	//			int curr_y	   = std::get<1>(m_ArrayIdx2VolCoords[ni]);
	//			int curr_z	   = std::get<2>(m_ArrayIdx2VolCoords[ni]);
	//			/*int xAux_begin = curr_x-CUBE_SIZE_FEATURE;		int xAux_end = curr_x+CUBE_SIZE_FEATURE;
	//			int yAux_begin = curr_y-CUBE_SIZE_FEATURE;		int yAux_end = curr_y+CUBE_SIZE_FEATURE;
	//			int zAux_begin = curr_z-CUBE_SIZE_FEATURE;		int zAux_end = curr_z+CUBE_SIZE_FEATURE;
	//			IndexedMesh3D::float3 sum = IndexedMesh3D::make_float3(0,0,0);
	//			for(int xAux = xAux_begin; xAux <= xAux_end; xAux++) {
	//				for(int yAux = yAux_begin; yAux <= yAux_end; yAux++) {
	//					for(int zAux = zAux_begin; zAux <= zAux_end; zAux++) {
	//						if(xAux >= 0 && xAux < m_W && yAux >= 0 && yAux < m_L && zAux >= 0 && zAux < m_H) {
	//							int voxel =  zAux*m_WxL + xAux*m_W + yAux;
	//							const float3& currentNormal = normals_Vol[voxel];
	//							if(currentNormal.x == 2 && currentNormal.y == 2 && currentNormal.z == 2) {
	//								sum += IndexedMesh3D::make_float3(1, 1, 1);
	//							} else if(currentNormal.x == -2 && currentNormal.y == -2 && currentNormal.z == -2) {
	//								sum += IndexedMesh3D::make_float3(-1, -1, -1);
	//							} else {
	//								sum += currentNormal;
	//							}
	//						}
	//					}
	//				}
	//			}
	//			sumCube[ni] = sum;*/
	//
	//			float3 sum2 = make_float3(0,0,0);
	//			std::vector<unsigned int>& cubeNeighbor_idx = m_cube_neighbors_idx_VolumeSub[ni];
	//			//BOOST_FOREACH(const int neighbor_idx_i, volumeNeighbor_idx){ const float3& currentNormal = normals_Vol[neighbor_idx_i];
	//			for (int nbi = 0; nbi < cubeNeighbor_idx.size(); ++nbi){								
	//				const float3& currentNormal = normals_Vol[cubeNeighbor_idx[nbi]];
	//				if(currentNormal.x == 2) {
	//					sum2 += IndexedMesh3D::make_float3(1, 1, 1);
	//				} else if(currentNormal.x == -2) {
	//					sum2 += IndexedMesh3D::make_float3(-1, -1, -1);
	//				} else {
	//					sum2 += currentNormal;
	//				}
	//			}
	//			sumCube[ni] = sum2;				
	//			
	//		}	

	//std::cout << "LCF part" << std::endl;

	// ##############################################################
	// 6. LCF part
#pragma omp parallel for 
	for (int vi = 0; vi < numVertex; ++vi){

		int ni = Vtx2Vxl[vi];
		int curr_x = std::get<0>(m_ArrayIdx2VolCoords[ni]);
		int curr_y = std::get<1>(m_ArrayIdx2VolCoords[ni]);
		int curr_z = std::get<2>(m_ArrayIdx2VolCoords[ni]);

		const float3& normal = normals_vtx[vi];
		float3 ref_frame_vector = make_float3(0, 0, 0);

		float NdotPoint = dot(normal, make_float3(curr_x, curr_y, curr_z));
		float maxValue = -STATS_HUGE_NUMBER;
		cv::Mat data(1, 3, CV_32F);


		for (std::vector<Neighbors>::const_iterator n_itr = m_cube_neighbors_idx_LCF.begin(); n_itr != m_cube_neighbors_idx_LCF.end(); n_itr++){
			int x_n = curr_x + std::get<0>(*n_itr);
			int y_n = curr_y + std::get<1>(*n_itr);
			int z_n = curr_z + std::get<2>(*n_itr);

			if (x_n >= 0 && x_n < m_L && y_n >= 0 && y_n < m_W && z_n >= 0 && z_n < m_H) {
				int voxel = z_n*m_WxL + x_n*m_W + y_n;
				float3 currentNormal = normals_Vol[voxel];

				if (currentNormal.x != 2 && currentNormal.x != -2) {
					if (m_LCF_FLAG == LCF::SIGN){

						float3 tempVoxel = make_float3(x_n, y_n, z_n);
						float dot_prod = dot(tempVoxel, normal);

						if ((dot_prod - NdotPoint) > maxValue) {			// compute the signed dist. 					
							maxValue = dot_prod - NdotPoint;
							ref_frame_vector = make_float3(std::get<0>(*n_itr), std::get<1>(*n_itr), std::get<2>(*n_itr));
						}
					}
					else if (m_LCF_FLAG == LCF::EVD){

						float3 tempVoxel = make_float3(x_n, y_n, z_n);
						float dot_prod = dot(tempVoxel, normal);
						float k = NdotPoint - dot_prod;
						float3 proj = tempVoxel + k*normal;

						data.at<float>(data.rows - 1, 0) = proj.x; data.at<float>(data.rows - 1, 1) = proj.y; data.at<float>(data.rows - 1, 2) = proj.z;
						cv::vconcat(data, cv::Mat::zeros(1, 3, CV_32F), data);


					}
					else{ // m_LCF_FLAG==LCF_MY
						if (abs(dot(std::get<4>(*n_itr), normal)) < 0.4){
							float dot_prod = dot(currentNormal, normal);
							if ((std::get<3>(*n_itr) + dot_prod) > maxValue) {
								maxValue = std::get<3>(*n_itr) + dot_prod;
								ref_frame_vector = make_float3(std::get<0>(*n_itr), std::get<1>(*n_itr), std::get<2>(*n_itr));
							}
						}
					}


				}
			}
		}

		if (m_LCF_FLAG == LCF::EVD){
			cv::Mat center = (cv::Mat_<float>(1, 3) << curr_x, curr_y, curr_z);

			data.pop_back();		//		std::cout<< "rows: " << data.rows << ", cols: " << data.cols << std::endl;
			//cv::PCA pca(data, Mat(), CV_PCA_DATA_AS_ROW, 1);		//		std::cout<< "pca done" << std::endl;
			cv::PCA pca(data, center, CV_PCA_DATA_AS_ROW, 1);		//		std::cout<< "pca done" << std::endl;
			ref_frame_vector = make_float3(pca.eigenvectors.at<float>(0, 0), pca.eigenvectors.at<float>(0, 1), pca.eigenvectors.at<float>(0, 2));
			/*std::cout<< "x: " << pca.eigenvectors << ", y: " << pca.eigenvectors.at<float>(0,1) << ", z: " << pca.eigenvectors.at<float>(0,2) << std::endl;
			system("pause");*/
		}
		ref_frame_vectors[vi] = ref_frame_vector;
	}

	/*std::cout << "(" << ref_frame_vectors[1].x << " ," << ref_frame_vectors[1].y << " ," << ref_frame_vectors[1].z << ")" << std::endl;
	system("pause");*/
}
#endif

void Voxelizer::voxelization(	const std::vector<float3>& coords_s_demean, const std::vector<Triangle>&	triangles_s,  
									const std::vector<float3>& normals_vtx,		const std::vector<float3>&		normals_tri, 
										std::vector<int>&	Vtx2Vxl, std::vector<float3>& normals_Vol, std::vector<float3>& sumCube, std::vector<float3>& ref_frame_vectors, 
										std::vector<float3>& limits){
		
		size_t numVertex = coords_s_demean.size();

		// ##############################################################
		// 1. create vertice-voxels mapping list
#pragma omp parallel for 
		for (int vi = 0; vi < numVertex; ++vi){
			int idxX = floor((coords_s_demean[vi].x-m_Xmin)/m_cubeSize[0]);		// start counting at 0.
			int idxY = floor((coords_s_demean[vi].y-m_Ymin)/m_cubeSize[1]);
			int idxZ = floor((coords_s_demean[vi].z-m_Zmin)/m_cubeSize[2]);
			
			assert(idxX < m_L);
			assert(idxY < m_W);
			assert(idxZ < m_H);
			if(idxX < limits[0].x) {
				limits[0].x = idxX;
			}
			if(idxX > limits[1].x) {
				limits[1].x = idxX;
			}
			if(idxY < limits[0].y) {
				limits[0].y = idxY;
			}
			if(idxY > limits[1].y) {
				limits[1].y = idxY;
			}
			if(idxZ < limits[0].z) {
				limits[0].z = idxZ;
			}
			if(idxZ > limits[1].z) {
				limits[1].z = idxZ;
			}

			Vtx2Vxl[vi] = idxZ*m_WxL + idxX*m_W + idxY;
		}


		// ##############################################################
		// 2. create shell by checking triangle intersection
		std::vector<float3> axes(13);
		axes[0] = make_float3(1.0f, 0.0f, 0.0f);		axes[1] = make_float3(0.0f, 1.0f, 0.0f);	axes[2] = make_float3(0.0f, 0.0f, 1.0f);
	
		std::vector<float3>::const_iterator tnl_itr = normals_tri.begin();
		cv::Mat AABC(2, 3, CV_32SC1);
		
		std::vector<cv::Mat>	occ_slices(m_H);	
#pragma omp parallel for 
		for (int zi = 0; zi < m_H; ++zi){
			occ_slices[zi].create(m_L,m_W,CV_8UC1);
			occ_slices[zi].setTo(255);
		}

		for(std::vector<Triangle>::const_iterator t_itr = triangles_s.begin(); t_itr != triangles_s.end(); ++t_itr){
			
			
			axes[3] = *tnl_itr++;

			std::vector<float3>	edges(3);
			edges.push_back(coords_s_demean[t_itr->v1] - coords_s_demean[t_itr->v0]);	// v1_v0		
			edges.push_back(coords_s_demean[t_itr->v2] - coords_s_demean[t_itr->v0]);	// v2_v0
			edges.push_back(coords_s_demean[t_itr->v2] - coords_s_demean[t_itr->v1]);	// v2_v1
		
			for (int ai = 0; ai < 9; ++ai)		axes[ai+4] = CC3D::cross(edges[ai/3], axes[ai%3]);
		

			// ##############################################################
			// 2.1 identify the cubes which intersect with the bounding box of triangle
			AABC.at<int>(0,0) = floor((std::min(std::min(coords_s_demean[t_itr->v1].x, coords_s_demean[t_itr->v0].x), coords_s_demean[t_itr->v2].x)-m_Xmin)/m_cubeSize[0]);
			AABC.at<int>(1,0) = floor((std::max(std::max(coords_s_demean[t_itr->v1].x, coords_s_demean[t_itr->v0].x), coords_s_demean[t_itr->v2].x)-m_Xmin)/m_cubeSize[0]);
			AABC.at<int>(0,1) = floor((std::min(std::min(coords_s_demean[t_itr->v1].y, coords_s_demean[t_itr->v0].y), coords_s_demean[t_itr->v2].y)-m_Ymin)/m_cubeSize[1]);
			AABC.at<int>(1,1) = floor((std::max(std::max(coords_s_demean[t_itr->v1].y, coords_s_demean[t_itr->v0].y), coords_s_demean[t_itr->v2].y)-m_Ymin)/m_cubeSize[1]);
			AABC.at<int>(0,2) = floor((std::min(std::min(coords_s_demean[t_itr->v1].z, coords_s_demean[t_itr->v0].z), coords_s_demean[t_itr->v2].z)-m_Zmin)/m_cubeSize[2]);
			AABC.at<int>(1,2) = floor((std::max(std::max(coords_s_demean[t_itr->v1].z, coords_s_demean[t_itr->v0].z), coords_s_demean[t_itr->v2].z)-m_Zmin)/m_cubeSize[2]);

			

			std::vector<cv::Point3i> AABB_intersec_cubes;

			for (int xi = AABC.at<int>(0,0); xi <= AABC.at<int>(1,0); ++xi){
				for (int yi = AABC.at<int>(0,1); yi <= AABC.at<int>(1,1); ++yi){
					for (int zi = AABC.at<int>(0,2); zi <= AABC.at<int>(1,2); ++zi){	
						cv::Point3i tmp(xi,yi,zi);
						AABB_intersec_cubes.push_back(tmp);
					}
				}
			}
		
			std::vector<float3>	corners_coords(8);
						
			// ##############################################################
			// 2.2 for all those cubes, perform SAT (Separating Axis Theorem) test 
			BOOST_FOREACH(const cv::Point3i& AABB_intersec_cubes_i, AABB_intersec_cubes){
				getCornerPointsCoords(AABB_intersec_cubes_i, m_cubeSize, m_Xmin, m_Ymin, m_Zmin, corners_coords);	

				float	proj_tri[3];		float	proj_cub[8];
				int count = 0;
				BOOST_FOREACH(const float3& axis_i, axes){
					proj_tri[0] = dot(axis_i, coords_s_demean[t_itr->v0]);
					proj_tri[1] = dot(axis_i, coords_s_demean[t_itr->v1]);
					proj_tri[2] = dot(axis_i, coords_s_demean[t_itr->v2]);
				
					//#pragma omp parallel for		// don't do parallel computing here. it doesn't help at all, even makes it (much) slower 
					for (int cri = 0; cri < 8; ++cri)
						proj_cub[cri] = dot(axis_i, corners_coords[cri]);
				
					auto MinMaxTri = std::minmax_element (proj_tri, proj_tri+3);
					auto MinMaxCub = std::minmax_element (proj_cub, proj_cub+8);

					if ((MinMaxCub.first < MinMaxTri.second) || (MinMaxTri.first < MinMaxCub.second)){
						if(AABB_intersec_cubes_i.z*m_WxL + AABB_intersec_cubes_i.x*m_W + AABB_intersec_cubes_i.y < normals_Vol.size() ){
							normals_Vol[AABB_intersec_cubes_i.z*m_WxL + AABB_intersec_cubes_i.x*m_W + AABB_intersec_cubes_i.y] += axes[3];
							occ_slices[AABB_intersec_cubes_i.z].at<unsigned char>(AABB_intersec_cubes_i.x, AABB_intersec_cubes_i.y) = 0;
							break;
						}
					}			

				}	
			}	
		}

		assert(tnl_itr==normals_tri.end());
		
		// ##############################################################
		// 3. normalization of normals
		#pragma omp parallel for 
		for (int ni = 0; ni < m_NUM_Vxls; ++ni){
			float norm = norm2(normals_Vol[ni]);
			if (norm!=0)	normals_Vol[ni] /=norm;			
		}
	
		if(m_SHOW_FLAG){
			cv::namedWindow("A", CV_WINDOW_AUTOSIZE);
			cvNamedWindow("A", 0);
			for (int zi = 0; zi < m_H; ++zi){
				cv::imshow("A", occ_slices[zi]);
				cv::waitKey(30);
			}
		}
	
		// ##############################################################
		// 4. distinguish between inside (-2) and outside (2), by check the connected compoment (a bit ad-hoc)
#pragma omp parallel for 
		for (int zi = 0; zi < m_H; ++zi){
			IplImage imgtmp = IplImage(occ_slices[zi]);
			int* labels= new int[m_WxL];

			int conncomp = bwlabel(&imgtmp, 4, labels);

			int zixWxL				= zi*m_WxL;
			int count				= 0;
			const int* const begin	= labels;
			const int* const end	= labels + m_WxL;

			for (const int* label = begin; label != end; ++label){		
				if (*label==conncomp){ normals_Vol[zixWxL + count] = make_float3(2.0f, 2.0f, 2.0f);}
				else if((*label < conncomp)&&(*label > 0)){ normals_Vol[zixWxL + count] =  make_float3(-2.0f, -2.0f, -2.0f);}
				else{}
				count++;
			}

			/*for (int ri = 0; ri < L; ++ri){
				for (int ci = 0; ci < W; ++ci){
					int label = labels[ri*W+ci];
					if (label==conncomp){ normals_Vol[zixWxL + ri*W + ci] = make_float3(2.0f, 2.0f, 2.0f);}
					else if((label < conncomp)&&(label > 0)){ normals_Vol[zixWxL + ri*W + ci] =  make_float3(-2.0f, -2.0f, -2.0f);}
					else{}			
				}		
			}*/

			delete[] labels;
		}
		
		/*std::ofstream f_out("normalx_afterConnComp");	
		for (int ri = 0; ri < m_L; ++ri){
			for (int ci = 0; ci < m_W; ++ci){
				f_out << normals_Vol[94*m_WxL + ri*m_L + ci].x << " ";
			}
			f_out << "\n";
		}
		f_out.close();*/

//		// ##############################################################
//		// 5. when needed, compute cubeSum (gotta replace this part with integral volume)
//#pragma omp parallel for 
//		for (int ni = 0; ni < m_NUM_Vxls; ++ni){
//			
//			int curr_x	   = std::get<0>(m_ArrayIdx2VolCoords[ni]);
//			int curr_y	   = std::get<1>(m_ArrayIdx2VolCoords[ni]);
//			int curr_z	   = std::get<2>(m_ArrayIdx2VolCoords[ni]);
//			/*int xAux_begin = curr_x-CUBE_SIZE_FEATURE;		int xAux_end = curr_x+CUBE_SIZE_FEATURE;
//			int yAux_begin = curr_y-CUBE_SIZE_FEATURE;		int yAux_end = curr_y+CUBE_SIZE_FEATURE;
//			int zAux_begin = curr_z-CUBE_SIZE_FEATURE;		int zAux_end = curr_z+CUBE_SIZE_FEATURE;
//			IndexedMesh3D::float3 sum = IndexedMesh3D::make_float3(0,0,0);
//			for(int xAux = xAux_begin; xAux <= xAux_end; xAux++) {
//				for(int yAux = yAux_begin; yAux <= yAux_end; yAux++) {
//					for(int zAux = zAux_begin; zAux <= zAux_end; zAux++) {
//						if(xAux >= 0 && xAux < m_W && yAux >= 0 && yAux < m_L && zAux >= 0 && zAux < m_H) {
//							int voxel =  zAux*m_WxL + xAux*m_W + yAux;
//							const float3& currentNormal = normals_Vol[voxel];
//							if(currentNormal.x == 2 && currentNormal.y == 2 && currentNormal.z == 2) {
//								sum += IndexedMesh3D::make_float3(1, 1, 1);
//							} else if(currentNormal.x == -2 && currentNormal.y == -2 && currentNormal.z == -2) {
//								sum += IndexedMesh3D::make_float3(-1, -1, -1);
//							} else {
//								sum += currentNormal;
//							}
//						}
//					}
//				}
//			}
//			sumCube[ni] = sum;*/
//
//			float3 sum2 = make_float3(0,0,0);
//			std::vector<unsigned int>& cubeNeighbor_idx = m_cube_neighbors_idx_VolumeSub[ni];
//			//BOOST_FOREACH(const int neighbor_idx_i, volumeNeighbor_idx){ const float3& currentNormal = normals_Vol[neighbor_idx_i];
//			for (int nbi = 0; nbi < cubeNeighbor_idx.size(); ++nbi){								
//				const float3& currentNormal = normals_Vol[cubeNeighbor_idx[nbi]];
//				if(currentNormal.x == 2) {
//					sum2 += IndexedMesh3D::make_float3(1, 1, 1);
//				} else if(currentNormal.x == -2) {
//					sum2 += IndexedMesh3D::make_float3(-1, -1, -1);
//				} else {
//					sum2 += currentNormal;
//				}
//			}
//			sumCube[ni] = sum2;				
//			
//		}	

		//std::cout << "LCF part" << std::endl;

		// ##############################################################
		// 6. LCF part
#pragma omp parallel for 
		for (int vi = 0; vi < numVertex; ++vi){

			int ni		   = Vtx2Vxl[vi];			
			int curr_x	   = std::get<0>(m_ArrayIdx2VolCoords[ni]);		
			int curr_y	   = std::get<1>(m_ArrayIdx2VolCoords[ni]);
			int curr_z	   = std::get<2>(m_ArrayIdx2VolCoords[ni]);

			const float3& normal	= normals_vtx[vi];	
			float3 ref_frame_vector = make_float3(0,0,0);

			float NdotPoint = dot(normal, make_float3(curr_x, curr_y, curr_z));
			float maxValue = -STATS_HUGE_NUMBER;
			cv::Mat data(1,3,CV_32F);
			

			for (std::vector<Neighbors>::const_iterator n_itr = m_cube_neighbors_idx_LCF.begin(); n_itr != m_cube_neighbors_idx_LCF.end(); n_itr++){
				int x_n = curr_x + std::get<0>(*n_itr);
				int y_n = curr_y + std::get<1>(*n_itr);
				int z_n = curr_z + std::get<2>(*n_itr);

				if(x_n >= 0 && x_n < m_L && y_n >= 0 && y_n < m_W && z_n >= 0 && z_n < m_H) {
					int voxel =  z_n*m_WxL + x_n*m_W + y_n;
					float3 currentNormal = normals_Vol[voxel];

					if(currentNormal.x!=2 && currentNormal.x!=-2 ) {
						if (m_LCF_FLAG == LCF::SIGN){							

							float3 tempVoxel = make_float3(x_n, y_n, z_n);
							float dot_prod = dot(tempVoxel, normal);

							if((dot_prod - NdotPoint) > maxValue) {			// compute the signed dist. 					
								maxValue = dot_prod - NdotPoint;
								ref_frame_vector = make_float3( std::get<0>(*n_itr), std::get<1>(*n_itr), std::get<2>(*n_itr));			
							}
						}
						else if (m_LCF_FLAG == LCF::EVD){

							float3 tempVoxel = make_float3(x_n, y_n, z_n);
							float dot_prod = dot(tempVoxel, normal);
							float k = NdotPoint - dot_prod;
							float3 proj = tempVoxel + k*normal;
																
							data.at<float>(data.rows-1,0) = proj.x; data.at<float>(data.rows-1,1) = proj.y; data.at<float>(data.rows-1,2) = proj.z;
							cv::vconcat(data, cv::Mat::zeros(1, 3, CV_32F), data);		

						
						}
						else{ // m_LCF_FLAG==LCF_MY
							if(abs(dot(std::get<4>(*n_itr), normal)) < 0.4){
								float dot_prod = dot(currentNormal, normal);
								if((std::get<3>(*n_itr) + dot_prod) > maxValue) {								
									maxValue = std::get<3>(*n_itr) + dot_prod;
									ref_frame_vector = make_float3( std::get<0>(*n_itr), std::get<1>(*n_itr), std::get<2>(*n_itr));									
								}
							}
						}


					}
				}
			}

			if (m_LCF_FLAG == LCF::EVD){
				cv::Mat center = (cv::Mat_<float>(1, 3) << curr_x, curr_y, curr_z);

				data.pop_back();		//		std::cout<< "rows: " << data.rows << ", cols: " << data.cols << std::endl;
				//cv::PCA pca(data, Mat(), CV_PCA_DATA_AS_ROW, 1);		//		std::cout<< "pca done" << std::endl;
				cv::PCA pca(data, center, CV_PCA_DATA_AS_ROW, 1);		//		std::cout<< "pca done" << std::endl;
				ref_frame_vector = make_float3(pca.eigenvectors.at<float>(0,0), pca.eigenvectors.at<float>(0,1), pca.eigenvectors.at<float>(0,2));		
				/*std::cout<< "x: " << pca.eigenvectors << ", y: " << pca.eigenvectors.at<float>(0,1) << ", z: " << pca.eigenvectors.at<float>(0,2) << std::endl;
				system("pause");*/
			}
			ref_frame_vectors[vi] = ref_frame_vector;															
		}				

		/*std::cout << "(" << ref_frame_vectors[1].x << " ," << ref_frame_vectors[1].y << " ," << ref_frame_vectors[1].z << ")" << std::endl;
		system("pause");*/
	}

void Voxelizer::voxelization(const std::vector<float3>& coords_s_demean, const std::vector<Triangle>&	triangles_s,
		const std::vector<float3>& normals_Sur, const std::vector<float3>&		normals_tri,
		std::vector<int>&	Vtx2Vxl, std::vector<float3>& normals_Vol, std::vector<float3>& sumCube, std::vector<std::vector<float3>>& ref_frame_vectors,
		std::vector<float3>& limits){

		size_t numVertex = coords_s_demean.size();

		// ##############################################################
		// 1. create vertice-voxels mapping list
#pragma omp parallel for 
		for (int vi = 0; vi < numVertex; ++vi){
			int idxX = floor((coords_s_demean[vi].x - m_Xmin) / m_cubeSize[0]);		// start counting at 0.
			int idxY = floor((coords_s_demean[vi].y - m_Ymin) / m_cubeSize[1]);
			int idxZ = floor((coords_s_demean[vi].z - m_Zmin) / m_cubeSize[2]);

			assert(idxX < m_L);
			assert(idxY < m_W);
			assert(idxZ < m_H);
			if (idxX < limits[0].x) {
				limits[0].x = idxX;
			}
			if (idxX > limits[1].x) {
				limits[1].x = idxX;
			}
			if (idxY < limits[0].y) {
				limits[0].y = idxY;
			}
			if (idxY > limits[1].y) {
				limits[1].y = idxY;
			}
			if (idxZ < limits[0].z) {
				limits[0].z = idxZ;
			}
			if (idxZ > limits[1].z) {
				limits[1].z = idxZ;
			}

			Vtx2Vxl[vi] = idxZ*m_WxL + idxX*m_W + idxY;
		}


		// ##############################################################
		// 2. create shell by checking triangle intersection
		std::vector<float3> axes(13);
		axes[0] = make_float3(1.0f, 0.0f, 0.0f);		axes[1] = make_float3(0.0f, 1.0f, 0.0f);	axes[2] = make_float3(0.0f, 0.0f, 1.0f);

		std::vector<float3>::const_iterator tnl_itr = normals_tri.begin();
		cv::Mat AABC(2, 3, CV_32SC1);

		std::vector<cv::Mat>	occ_slices(m_H);
#pragma omp parallel for 
		for (int zi = 0; zi < m_H; ++zi){
			occ_slices[zi].create(m_L, m_W, CV_8UC1);
			occ_slices[zi].setTo(255);
		}

		for (std::vector<Triangle>::const_iterator t_itr = triangles_s.begin(); t_itr != triangles_s.end(); ++t_itr){


			axes[3] = *tnl_itr++;

			std::vector<float3>	edges(3);
			edges.push_back(coords_s_demean[t_itr->v1] - coords_s_demean[t_itr->v0]);	// v1_v0		
			edges.push_back(coords_s_demean[t_itr->v2] - coords_s_demean[t_itr->v0]);	// v2_v0
			edges.push_back(coords_s_demean[t_itr->v2] - coords_s_demean[t_itr->v1]);	// v2_v1

			for (int ai = 0; ai < 9; ++ai)		axes[ai + 4] = CC3D::cross(edges[ai / 3], axes[ai % 3]);


			// ##############################################################
			// 2.1 identify the cubes which intersect with the bounding box of triangle
			AABC.at<int>(0, 0) = floor((std::min(std::min(coords_s_demean[t_itr->v1].x, coords_s_demean[t_itr->v0].x), coords_s_demean[t_itr->v2].x) - m_Xmin) / m_cubeSize[0]);
			AABC.at<int>(1, 0) = floor((std::max(std::max(coords_s_demean[t_itr->v1].x, coords_s_demean[t_itr->v0].x), coords_s_demean[t_itr->v2].x) - m_Xmin) / m_cubeSize[0]);
			AABC.at<int>(0, 1) = floor((std::min(std::min(coords_s_demean[t_itr->v1].y, coords_s_demean[t_itr->v0].y), coords_s_demean[t_itr->v2].y) - m_Ymin) / m_cubeSize[1]);
			AABC.at<int>(1, 1) = floor((std::max(std::max(coords_s_demean[t_itr->v1].y, coords_s_demean[t_itr->v0].y), coords_s_demean[t_itr->v2].y) - m_Ymin) / m_cubeSize[1]);
			AABC.at<int>(0, 2) = floor((std::min(std::min(coords_s_demean[t_itr->v1].z, coords_s_demean[t_itr->v0].z), coords_s_demean[t_itr->v2].z) - m_Zmin) / m_cubeSize[2]);
			AABC.at<int>(1, 2) = floor((std::max(std::max(coords_s_demean[t_itr->v1].z, coords_s_demean[t_itr->v0].z), coords_s_demean[t_itr->v2].z) - m_Zmin) / m_cubeSize[2]);



			std::vector<cv::Point3i> AABB_intersec_cubes;

			for (int xi = AABC.at<int>(0, 0); xi <= AABC.at<int>(1, 0); ++xi){
				for (int yi = AABC.at<int>(0, 1); yi <= AABC.at<int>(1, 1); ++yi){
					for (int zi = AABC.at<int>(0, 2); zi <= AABC.at<int>(1, 2); ++zi){
						cv::Point3i tmp(xi, yi, zi);
						AABB_intersec_cubes.push_back(tmp);
					}
				}
			}

			std::vector<float3>	corners_coords(8);

			// ##############################################################
			// 2.2 for all those cubes, perform SAT (Separating Axis Theorem) test 
			BOOST_FOREACH(const cv::Point3i& AABB_intersec_cubes_i, AABB_intersec_cubes){
				getCornerPointsCoords(AABB_intersec_cubes_i, m_cubeSize, m_Xmin, m_Ymin, m_Zmin, corners_coords);

				float	proj_tri[3];		float	proj_cub[8];
				int count = 0;
				BOOST_FOREACH(const float3& axis_i, axes){
					proj_tri[0] = dot(axis_i, coords_s_demean[t_itr->v0]);
					proj_tri[1] = dot(axis_i, coords_s_demean[t_itr->v1]);
					proj_tri[2] = dot(axis_i, coords_s_demean[t_itr->v2]);

					//#pragma omp parallel for		// don't do parallel computing here. it doesn't help at all, even makes it (much) slower 
					for (int cri = 0; cri < 8; ++cri)
						proj_cub[cri] = dot(axis_i, corners_coords[cri]);

					auto MinMaxTri = std::minmax_element(proj_tri, proj_tri + 3);
					auto MinMaxCub = std::minmax_element(proj_cub, proj_cub + 8);

					if ((MinMaxCub.first < MinMaxTri.second) || (MinMaxTri.first < MinMaxCub.second)){
						if (AABB_intersec_cubes_i.z*m_WxL + AABB_intersec_cubes_i.x*m_W + AABB_intersec_cubes_i.y < normals_Vol.size()){
							normals_Vol[AABB_intersec_cubes_i.z*m_WxL + AABB_intersec_cubes_i.x*m_W + AABB_intersec_cubes_i.y] += axes[3];
							occ_slices[AABB_intersec_cubes_i.z].at<unsigned char>(AABB_intersec_cubes_i.x, AABB_intersec_cubes_i.y) = 0;
							break;
						}
					}

				}
			}
		}

		assert(tnl_itr == normals_tri.end());

		// ##############################################################
		// 3. normalization of normals
#pragma omp parallel for 
		for (int ni = 0; ni < m_NUM_Vxls; ++ni){
			float norm = norm2(normals_Vol[ni]);
			if (norm != 0)	normals_Vol[ni] /= norm;
		}

		if (m_SHOW_FLAG){
			cv::namedWindow("A", CV_WINDOW_AUTOSIZE);
			cvNamedWindow("A", 0);
			for (int zi = 0; zi < m_H; ++zi){
				cv::imshow("A", occ_slices[zi]);
				cv::waitKey(30);
			}
		}

		// ##############################################################
		// 4. distinguish between inside (-2) and outside (2), by check the connected compoment (a bit ad-hoc)
#pragma omp parallel for 
		for (int zi = 0; zi < m_H; ++zi){
			IplImage imgtmp = IplImage(occ_slices[zi]);
			int* labels = new int[m_WxL];

			int conncomp = bwlabel(&imgtmp, 4, labels);

			int zixWxL = zi*m_WxL;
			int count = 0;
			const int* const begin = labels;
			const int* const end = labels + m_WxL;

			for (const int* label = begin; label != end; ++label){
				if (*label == conncomp){ normals_Vol[zixWxL + count] = make_float3(2.0f, 2.0f, 2.0f); }
				else if ((*label < conncomp) && (*label > 0)){ normals_Vol[zixWxL + count] = make_float3(-2.0f, -2.0f, -2.0f); }
				else{}
				count++;
			}

			/*for (int ri = 0; ri < L; ++ri){
				for (int ci = 0; ci < W; ++ci){
					int label = labels[ri*W+ci];
					if (label==conncomp){ normals_Vol[zixWxL + ri*W + ci] = make_float3(2.0f, 2.0f, 2.0f);}
					else if((label < conncomp)&&(label > 0)){ normals_Vol[zixWxL + ri*W + ci] =  make_float3(-2.0f, -2.0f, -2.0f);}
					else{}
				}
			}*/

			delete[] labels;
		}

		/*std::ofstream f_out("normalx_afterConnComp");
		for (int ri = 0; ri < m_L; ++ri){
		for (int ci = 0; ci < m_W; ++ci){
		f_out << normals_Vol[94*m_WxL + ri*m_L + ci].x << " ";
		}
		f_out << "\n";
		}
		f_out.close();*/

//		// ##############################################################
//		// 5. when needed, compute cubeSum (gotta replace this part with integral volume)
//#pragma omp parallel for 
//		for (int ni = 0; ni < m_NUM_Vxls; ++ni){
//
//			int curr_x = std::get<0>(m_ArrayIdx2VolCoords[ni]);
//			int curr_y = std::get<1>(m_ArrayIdx2VolCoords[ni]);
//			int curr_z = std::get<2>(m_ArrayIdx2VolCoords[ni]);
//			/*int xAux_begin = curr_x-CUBE_SIZE_FEATURE;		int xAux_end = curr_x+CUBE_SIZE_FEATURE;
//			int yAux_begin = curr_y-CUBE_SIZE_FEATURE;		int yAux_end = curr_y+CUBE_SIZE_FEATURE;
//			int zAux_begin = curr_z-CUBE_SIZE_FEATURE;		int zAux_end = curr_z+CUBE_SIZE_FEATURE;
//			IndexedMesh3D::float3 sum = IndexedMesh3D::make_float3(0,0,0);
//			for(int xAux = xAux_begin; xAux <= xAux_end; xAux++) {
//			for(int yAux = yAux_begin; yAux <= yAux_end; yAux++) {
//			for(int zAux = zAux_begin; zAux <= zAux_end; zAux++) {
//			if(xAux >= 0 && xAux < m_W && yAux >= 0 && yAux < m_L && zAux >= 0 && zAux < m_H) {
//			int voxel =  zAux*m_WxL + xAux*m_W + yAux;
//			const float3& currentNormal = normals_Vol[voxel];
//			if(currentNormal.x == 2 && currentNormal.y == 2 && currentNormal.z == 2) {
//			sum += IndexedMesh3D::make_float3(1, 1, 1);
//			} else if(currentNormal.x == -2 && currentNormal.y == -2 && currentNormal.z == -2) {
//			sum += IndexedMesh3D::make_float3(-1, -1, -1);
//			} else {
//			sum += currentNormal;
//			}
//			}
//			}
//			}
//			}
//			sumCube[ni] = sum;*/
//
//			float3 sum2 = make_float3(0, 0, 0);
//			std::vector<unsigned int>& cubeNeighbor_idx = m_cube_neighbors_idx_VolumeSub[ni];
//			//BOOST_FOREACH(const int neighbor_idx_i, volumeNeighbor_idx){ const float3& currentNormal = normals_Vol[neighbor_idx_i];
//			for (int nbi = 0; nbi < cubeNeighbor_idx.size(); ++nbi){
//				const float3& currentNormal = normals_Vol[cubeNeighbor_idx[nbi]];
//				if (currentNormal.x == 2) {
//					sum2 += IndexedMesh3D::make_float3(1, 1, 1);
//				}
//				else if (currentNormal.x == -2) {
//					sum2 += IndexedMesh3D::make_float3(-1, -1, -1);
//				}
//				else {
//					sum2 += currentNormal;
//				}
//			}
//			sumCube[ni] = sum2;
//
//		}


		// ##############################################################
		// 6. LCF part

		//std::ofstream f_out("local profile7");
#pragma omp parallel for 		
		for (int vi = 0; vi < numVertex; ++vi){

			int ni = Vtx2Vxl[vi];
			int curr_x = std::get<0>(m_ArrayIdx2VolCoords[ni]);
			int curr_y = std::get<1>(m_ArrayIdx2VolCoords[ni]);
			int curr_z = std::get<2>(m_ArrayIdx2VolCoords[ni]);

			const float3& normal = normals_Sur[vi];
			float3 ref_v_1 = make_float3(normal.y, -normal.x, 0);		normalize(ref_v_1);
			float3 ref_v_2 = cross(normal, ref_v_1);

			std::vector<std::tuple<float, float, bool, int, int, int>>	projAngle_SignDist_isMax_offsets;


			float3 ref_frame_vector = make_float3(0, 0, 0);

			const float NdotPoint = dot(normal, make_float3(curr_x, curr_y, curr_z));
			float maxValue = -STATS_HUGE_NUMBER;
			float minValue = STATS_HUGE_NUMBER;
			cv::Mat data(1, 3, CV_32F);
			
			for (std::vector<Neighbors>::const_iterator n_itr = m_cube_neighbors_idx_LCF.begin(); n_itr != m_cube_neighbors_idx_LCF.end(); n_itr++){
				int x_n = curr_x + std::get<0>(*n_itr);
				int y_n = curr_y + std::get<1>(*n_itr);
				int z_n = curr_z + std::get<2>(*n_itr);

				if (x_n >= 0 && x_n < m_L && y_n >= 0 && y_n < m_W && z_n >= 0 && z_n < m_H) {
					int voxel = z_n*m_WxL + x_n*m_W + y_n;
					float3 currentNormal = normals_Vol[voxel];

					if (currentNormal.x != 2 && currentNormal.x != -2) {						

						if (m_LCF_FLAG <= LCF::SIGN){  // for signDist, x_n, y_n, z_n lie on the border of the cuboids		

							float3 tempVoxel = make_float3(x_n, y_n, z_n);
							
							float x_proj = dot(tempVoxel, ref_v_1);
							float y_proj = dot(tempVoxel, ref_v_2);

							float dot_prod = dot(tempVoxel, normal);
							float signDist_tmp = dot_prod - NdotPoint;
							
							projAngle_SignDist_isMax_offsets.push_back(std::make_tuple(atan2(y_proj, x_proj), signDist_tmp, true, std::get<0>(*n_itr), std::get<1>(*n_itr), std::get<2>(*n_itr)));

							if (signDist_tmp > maxValue) {			
								maxValue = signDist_tmp;
								ref_frame_vector = make_float3(std::get<0>(*n_itr), std::get<1>(*n_itr), std::get<2>(*n_itr));
							}

							if (signDist_tmp < minValue) 			minValue = signDist_tmp;
							
						}
						else if (m_LCF_FLAG == LCF::EVD){

							float3 tempVoxel = make_float3(x_n, y_n, z_n);
							float dot_prod = dot(tempVoxel, normal);
							float k = NdotPoint - dot_prod;
							float3 proj = tempVoxel + k*normal;

							data.at<float>(data.rows - 1, 0) = proj.x; data.at<float>(data.rows - 1, 1) = proj.y; data.at<float>(data.rows - 1, 2) = proj.z;
							cv::vconcat(data, cv::Mat::zeros(1, 3, CV_32F), data);


						}
						else{ // m_LCF_FLAG==LCF_MY, my old CVPR approaches
							if (abs(dot(std::get<4>(*n_itr), normal)) < 0.4){
								float dot_prod = dot(currentNormal, normal);
								if ((std::get<3>(*n_itr) + dot_prod) > maxValue) {
									maxValue = std::get<3>(*n_itr) + dot_prod;
									ref_frame_vector = make_float3(std::get<0>(*n_itr), std::get<1>(*n_itr), std::get<2>(*n_itr));
								}
							}
						}
					}
				}
			}

			

			if (m_LCF_FLAG == LCF::EVD){
				cv::Mat center = (cv::Mat_<float>(1, 3) << curr_x, curr_y, curr_z);

				data.pop_back();						
				cv::PCA pca(data, center, CV_PCA_DATA_AS_ROW, 1);		
				ref_frame_vector = make_float3(pca.eigenvectors.at<float>(0, 0), pca.eigenvectors.at<float>(0, 1), pca.eigenvectors.at<float>(0, 2));
			}
			ref_frame_vectors[vi].push_back(ref_frame_vector);

			if (m_LCF_FLAG <= LCF::SIGN){ // adding more local max, as ref. vectors

				std::sort(projAngle_SignDist_isMax_offsets.begin(), projAngle_SignDist_isMax_offsets.end());

				int numSurV = projAngle_SignDist_isMax_offsets.size();		//f_out << numSurV << " ";
				int radiusNMS = cv::max(4, (int)(numSurV*0.075));		// radius: adaptive to the length of this 1D array, but no less than 4
				//std::cout << "radius: " << radiusNMS << std::endl;		

				double th = maxValue - 0.05*(maxValue-minValue);		// thresholds for the bandwidths


				for (int i = 0; i < numSurV; ++i){					

					float signDist_i = std::get<1>(projAngle_SignDist_isMax_offsets[i]);
					//f_out << std::get<0>(projAngle_SignDist_isMax_offsets[i]) << " " << signDist_i << " ";

					int beginNMS = cv::max(i - radiusNMS, 0);
					int endNMS = cv::min(i + radiusNMS, numSurV);

					if (signDist_i > th){
						for (int j = beginNMS; j < endNMS; ++j){	// non-maximum suppression
							if (std::get<1>(projAngle_SignDist_isMax_offsets[j]) > signDist_i){
								std::get<2>(projAngle_SignDist_isMax_offsets[i]) = false;
								break;
							}
						}
					}
					else{
						std::get<2>(projAngle_SignDist_isMax_offsets[i]) = false;
					}
				}
				//f_out << "\n";

				for (std::vector<std::tuple<float, float, bool, int, int, int>>::const_iterator itr = projAngle_SignDist_isMax_offsets.begin(); itr != projAngle_SignDist_isMax_offsets.end(); ++itr){
					if (std::get<2>(*itr) && (std::get<1>(*itr) != maxValue))
						ref_frame_vectors[vi].push_back(make_float3(std::get<3>(*itr), std::get<4>(*itr), std::get<5>(*itr)));
					
					if (ref_frame_vectors[vi].size() == NUM_LOCALMAX) break;
				}

			}
		}		
		
		//f_out.close();
		//system("pause");
	}
	

	inline void Voxelizer::getCornerPointsCoords(const cv::Point3i& cubi, const float* const cubeSize, 
										float Xmin, float Ymin, float Zmin, 
										std::vector<float3>& corners_coords){
	
		corners_coords[0].x = Xmin + cubi.x*cubeSize[0];				corners_coords[0].y = Ymin + cubi.y*cubeSize[1];				corners_coords[0].z = Zmin + cubi.z*cubeSize[2];				// 0, 0, 0
		corners_coords[1].x = Xmin + cubi.x*cubeSize[0] + cubeSize[0];	corners_coords[1].y = Ymin + cubi.y*cubeSize[1];				corners_coords[1].z = Zmin + cubi.z*cubeSize[2];				// 1, 0, 0
		corners_coords[2].x = Xmin + cubi.x*cubeSize[0] + cubeSize[0];	corners_coords[2].y = Ymin + cubi.y*cubeSize[1] + cubeSize[1];	corners_coords[2].z = Zmin + cubi.z*cubeSize[2];				// 1, 1, 0
		corners_coords[3].x = Xmin + cubi.x*cubeSize[0];				corners_coords[3].y = Ymin + cubi.y*cubeSize[1] + cubeSize[1];	corners_coords[3].z = Zmin + cubi.z*cubeSize[2];				// 0, 1, 0
		corners_coords[4].x = Xmin + cubi.x*cubeSize[0];				corners_coords[4].y = Ymin + cubi.y*cubeSize[1];				corners_coords[4].z = Zmin + cubi.z*cubeSize[2] + cubeSize[2];	// 0, 0, 1
		corners_coords[5].x = Xmin + cubi.x*cubeSize[0] + cubeSize[0];	corners_coords[5].y = Ymin + cubi.y*cubeSize[1];				corners_coords[5].z = Zmin + cubi.z*cubeSize[2] + cubeSize[2];	// 1, 0, 1
		corners_coords[6].x = Xmin + cubi.x*cubeSize[0] + cubeSize[0];	corners_coords[6].y = Ymin + cubi.y*cubeSize[1] + cubeSize[1];	corners_coords[6].z = Zmin + cubi.z*cubeSize[2] + cubeSize[2];	// 1, 1, 1
		corners_coords[7].x = Xmin + cubi.x*cubeSize[0];				corners_coords[7].y = Ymin + cubi.y*cubeSize[1] + cubeSize[1];	corners_coords[7].z = Zmin + cubi.z*cubeSize[2] + cubeSize[2];	// 0, 1, 1
	}

	int Voxelizer::bwlabel(const IplImage* const img, int n, int* labels)  
	{  /* C++ equivalent of bwlabel in MATLAB, do know if there's a bug.... at least so far it serves me well....
		source: http://blog.csdn.net/yysdsyl/article/details/5343075 */

		if(n != 4 && n != 8)  n = 4;  
		int nr = img->height;  
		int nc = img->width;  
		int total = nr * nc;  
		// results  
		memset(labels, 0, total * sizeof(int));  
		int nobj = 0;                               // number of objects found in image  
		// other variables                               
		int* lset = new int[total];   // label table  
		memset(lset, 0, total * sizeof(int));  
		int ntable = 0;  
		for( int r = 0; r < nr; r++ ) {  
			for( int c = 0; c < nc; c++ )   
			{              
				if ( ELEM(img, r, c) )   // if A is an object  
				{                 
					// get the neighboring pixels B, C, D, and E  
					int B, C, D, E;  
					if ( c == 0 )   
						B = 0;   
					else   
						B = find( lset, ONETWO(labels, r, c - 1, nc) );  
					if ( r == 0 )   
						C = 0;   
					else   
						C = find( lset, ONETWO(labels, r - 1, c, nc) );  
					if ( r == 0 || c == 0 )   
						D = 0;   
					else   
						D = find( lset, ONETWO(labels, r - 1, c - 1, nc) );  
					if ( r == 0 || c == nc - 1 )   
						E = 0;  
					else   
						E = find( lset, ONETWO(labels, r - 1, c + 1, nc) );  
					if ( n == 4 )   
					{  
						// apply 4 connectedness  
						if ( B && C )   
						{        // B and C are labeled  
							if ( B == C )  
								ONETWO(labels, r, c, nc) = B;  
							else {  
								lset[C] = B;  
								ONETWO(labels, r, c, nc) = B;  
							}  
						}   
						else if ( B )             // B is object but C is not  
							ONETWO(labels, r, c, nc) = B;  
						else if ( C )               // C is object but B is not  
							ONETWO(labels, r, c, nc) = C;  
						else   
						{                      // B, C, D not object - new object  
							//   label and put into table  
							ntable++;  
							ONETWO(labels, r, c, nc) = lset[ ntable ] = ntable;  
						}  
					}   
					else if ( n == 6 )   
					{  
						// apply 6 connected ness  
						if ( D )                    // D object, copy label and move on  
							ONETWO(labels, r, c, nc) = D;  
						else if ( B && C )   
						{        // B and C are labeled  
							if ( B == C )  
								ONETWO(labels, r, c, nc) = B;  
							else   
							{  
								int tlabel = MIN(B,C);  
								lset[B] = tlabel;  
								lset[C] = tlabel;  
								ONETWO(labels, r, c, nc) = tlabel;  
							}  
						}   
						else if ( B )             // B is object but C is not  
							ONETWO(labels, r, c, nc) = B;  
						else if ( C )               // C is object but B is not  
							ONETWO(labels, r, c, nc) = C;  
						else   
						{                      // B, C, D not object - new object  
							//   label and put into table  
							ntable++;  
							ONETWO(labels, r, c, nc) = lset[ ntable ] = ntable;  
						}  
					}  
					else if ( n == 8 )   
					{  
						// apply 8 connectedness  
						if ( B || C || D || E )   
						{  
							int tlabel = B;  
							if ( B )   
								tlabel = B;  
							else if ( C )   
								tlabel = C;  
							else if ( D )   
								tlabel = D;  
							else if ( E )   
								tlabel = E;  
							ONETWO(labels, r, c, nc) = tlabel;  
							if ( B && B != tlabel )   
								lset[B] = tlabel;  
							if ( C && C != tlabel )   
								lset[C] = tlabel;  
							if ( D && D != tlabel )   
								lset[D] = tlabel;  
							if ( E && E != tlabel )   
								lset[E] = tlabel;  
						}   
						else   
						{  
							//   label and put into table  
							ntable++;  
							ONETWO(labels, r, c, nc) = lset[ ntable ] = ntable;  
						}  
					}  
				}   
				else   
				{  
					ONETWO(labels, r, c, nc) = NO_OBJECT;      // A is not an object so leave it  
				}  
			}  
		}  
		// consolidate component table  
		for( int i = 0; i <= ntable; i++ )  
			lset[i] = find( lset, i );                                                                                                   
		// run image through the look-up table  
		for( int r = 0; r < nr; r++ )  
			for( int c = 0; c < nc; c++ )  
				ONETWO(labels, r, c, nc) = lset[ ONETWO(labels, r, c, nc) ];  
		// count up the objects in the image  
		for( int i = 0; i <= ntable; i++ )  
			lset[i] = 0;  
		for( int r = 0; r < nr; r++ ){  
			for( int c = 0; c < nc; c++ ) { 
				lset[ ONETWO(labels, r, c, nc) ]++;  
			}
		}
		// number the objects from 1 through n objects  
		nobj = 0;  
		lset[0] = 0;  
		for( int i = 1; i <= ntable; i++ )  {
			if ( lset[i] > 0 )  {
				lset[i] = ++nobj;
			}
		}
		// run through the look-up table again  

		for( int r = 0; r < nr; r++ )  {
			for( int c = 0; c < nc; c++ )  {
				ONETWO(labels, r, c, nc) = lset[ ONETWO(labels, r, c, nc) ];  
			}
		}
		//  
		delete[] lset;  
		return nobj;  
	}  	
}