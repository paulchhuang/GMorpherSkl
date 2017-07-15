/* *************************************************
 * Copyright (2015) : Paul Huang
 * Implementing Ijaz's joint angle prior model
 * *************************************************/
#include <JointAnglePrior.h>


#define NUM_JOINT_SUB 16

namespace GMorpher 
{
	void gramschmidt(const float3& u, const float3& v, const float3& w, std::vector<float>& R){

		float norm_tmp = norm2(u);
		float3 B0 = make_float3(u.x / norm_tmp, u.y / norm_tmp, u.z / norm_tmp);

		float3 v_tmp = v;
		std::vector<float3> U(1, B0);

		double pc = dot(U[0], v_tmp);

		float3 p = pc*U[0];

		v_tmp = v_tmp - p;		normalize(v_tmp);



		float3 B1 = v_tmp;

		v_tmp = w;
		U.push_back(B1);

		pc = dot(U[0], v_tmp);
		p = pc*U[0];

		pc = dot(U[1], v_tmp);
		p += pc*U[1];

		v_tmp = v_tmp - p;		normalize(v_tmp);


		R[0] = B0.x;		R[3] = B1.x;		R[6] = v_tmp.x;
		R[1] = B0.y;		R[4] = B1.y;		R[7] = v_tmp.y;
		R[2] = B0.z;		R[5] = B1.z;		R[8] = v_tmp.z;
	}

	void findClosestValidPoint(const std::vector<double>& boundary, const std::tuple<double, double>& pt, std::tuple<double, double>& pt2){
		int num_Pts = boundary.size() / 2;

		std::vector<double> dist(num_Pts, 0);

		for (int pi = 0; pi < num_Pts; ++pi){
			double diff1 = boundary[pi * 2 + 0] - std::get<0>(pt);
			double diff2 = boundary[pi * 2 + 1] - std::get<1>(pt);

			dist[pi] = std::pow(diff1, 2) + std::pow(diff2, 2);
		}

		int ind = std::min_element(dist.begin(), dist.end()) - dist.begin();

		std::get<0>(pt2) = boundary[ind * 2 + 0];
		std::get<1>(pt2) = boundary[ind * 2 + 1];
	}

	inline void pointPlane2Sphere(const float3& Xh, const float3& v, const double ri, float3& X){
		double a = sqrt(1 - dot(Xh, Xh));
		X = ri*(Xh - v*a);
	}

	void estimateZ(const std::vector<float3>& dZ, const std::vector<double>& edges, const float3& Z1, std::vector<float3>& Z){
		int P = dZ.size() + 1;
		int numEdge = edges.size() / 2;

		Z.resize(P);

		Z[0] = Z1;

		for (int i = 0; i < numEdge; i++){		
			Z[edges[i + numEdge] - 1] = Z[edges[i] - 1] - dZ[i];
		}
	}


	JointAnglePrior::JointAnglePrior(const char* const fileName){

		m_BPName.insert(std::make_pair(1, "back-bone"));
		m_BPName.insert(std::make_pair(2, "R-shldr"));
		m_BPName.insert(std::make_pair(3, "R-Uarm"));
		m_BPName.insert(std::make_pair(4, "R-Larm"));
		m_BPName.insert(std::make_pair(5, "L-shldr"));
		m_BPName.insert(std::make_pair(6, "L-Uarm"));
		m_BPName.insert(std::make_pair(7, "L-Larm"));
		m_BPName.insert(std::make_pair(8, "head"));
		m_BPName.insert(std::make_pair(9, "R-hip"));
		m_BPName.insert(std::make_pair(10, "R-Uleg"));
		m_BPName.insert(std::make_pair(11, "R-Lleg"));
		m_BPName.insert(std::make_pair(12, "R-feet"));
		m_BPName.insert(std::make_pair(13, "L-hip"));
		m_BPName.insert(std::make_pair(14, "L-Uleg"));
		m_BPName.insert(std::make_pair(15, "L-Lleg"));
		m_BPName.insert(std::make_pair(16, "L-feet"));

		m_arrMap.insert(std::make_pair("edges", 0));
		m_arrMap.insert(std::make_pair("a", 1));
		m_arrMap.insert(std::make_pair("di", 2));
		m_arrMap.insert(std::make_pair("jmp", 3));
		m_arrMap.insert(std::make_pair("chlds", 4));
		m_arrMap.insert(std::make_pair("prnts", 5));

		m_BcellMap.insert(std::make_pair("angleSprd", 0));

		m_cellMap.insert(std::make_pair("sepPlane", 0));
		m_cellMap.insert(std::make_pair("bounds", 1));
		m_cellMap.insert(std::make_pair("E2", 2));
		m_cellMap.insert(std::make_pair("boundries", 3));

		m_arrayData.resize(m_arrMap.size());
		m_cellData.resize(m_cellMap.size());
		m_BcellData.resize(m_BcellMap.size());

		m_chldsT.push_back(3);
		m_chldsT.push_back(6);
		m_chldsT.push_back(8);
		m_chldsT.push_back(10);
		m_chldsT.push_back(14);

		m_jointMappingIjaz2Paul.insert(std::make_pair(1, 1));			// belly-torso
		m_jointMappingIjaz2Paul.insert(std::make_pair(10, 2));			// RHip-LHip
		m_jointMappingIjaz2Paul.insert(std::make_pair(11, 3));			// RKnee-LKnee
		m_jointMappingIjaz2Paul.insert(std::make_pair(12, 4));			// RAnkle-LFeet
		m_jointMappingIjaz2Paul.insert(std::make_pair(14, 5));			// LHip-RHip
		m_jointMappingIjaz2Paul.insert(std::make_pair(15, 6));			// LKnee-RKnee
		m_jointMappingIjaz2Paul.insert(std::make_pair(16, 7));			// LAnkle-RFeet
		m_jointMappingIjaz2Paul.insert(std::make_pair(2, 8));			// Neck-Neck
		m_jointMappingIjaz2Paul.insert(std::make_pair(9, 9));			// Face-Head
		m_jointMappingIjaz2Paul.insert(std::make_pair(3, 10));			// RShldr-LShldr
		m_jointMappingIjaz2Paul.insert(std::make_pair(4, 11));			// RElbow-LElbow
		m_jointMappingIjaz2Paul.insert(std::make_pair(5, 12));			// RWrist-LWrist
		m_jointMappingIjaz2Paul.insert(std::make_pair(6, 13));			// LShldr-RShldr
		m_jointMappingIjaz2Paul.insert(std::make_pair(7, 14));			// LElbow-RElbow
		m_jointMappingIjaz2Paul.insert(std::make_pair(8, 15));			// LWrist-RWrist


		const char** varNames;
		int			 numVar;

		// open MAT-file
		MATFile *pmat = matOpen(fileName, "r");
		if (pmat == NULL) { printf("Error opening file %s\n", fileName); return; }
		varNames = (const char **)matGetDir(pmat, &numVar);
		std::cout << "there are " << numVar << " variables:" << std::endl;
		//for (int i=0; i < numVar; i++)	printf("%s\n",varNames[i]);


		mxArray *arr;
		mxClassID category;

		for (int vi = 0; vi < numVar; ++vi){
			arr = matGetVariable(pmat, varNames[vi]);
			if (arr == NULL) {
				printf("Error reading in file %s\n", varNames[vi]);
				return;
			}
			//if(vi==3) continue;
			category = mxGetClassID(arr);
			mwSize num = (mwSize)mxGetNumberOfElements(arr);

			switch (category) {
			case mxCELL_CLASS:
				printf("variable %s is a cell, with %d elements\n", varNames[vi], num);
				mxArray *ePtr;

				if (vi == 3)		// the 4th variable is a logical cell, this is hard coded.
					m_BcellData[m_BcellMap[std::string(varNames[vi])]].resize(num);
				else
					m_cellData[m_cellMap[std::string(varNames[vi])]].resize(num);

				for (mwIndex ei = 0; ei < num; ei++){
					ePtr = mxGetCell(arr, ei);
					if (ePtr == NULL) {
						//printf("\tEmpty Cell\n");							
					}

					else {
						mwSize num_ele = (mwSize)mxGetNumberOfElements(ePtr);
						//printf("\tnumeric class with number of elements: %d\n", num_ele);

						if (mxGetClassID(ePtr) == mxLOGICAL_CLASS){
							bool *pr = mxGetLogicals(ePtr);
							if (pr != NULL) {
								m_BcellData[m_BcellMap[std::string(varNames[vi])]][ei].resize(num_ele);
								m_BcellData[m_BcellMap[std::string(varNames[vi])]][ei].assign(pr, pr + num_ele);
							}
						}
						else{

							double *pr = mxGetPr(ePtr);
							if (pr != NULL) {
								m_cellData[m_cellMap[std::string(varNames[vi])]][ei].resize(num_ele);
								m_cellData[m_cellMap[std::string(varNames[vi])]][ei].assign(pr, pr + num_ele);
							}
						}
					}
					//mxDestroyArray(ePtr);
				}

				break;
			default:
				printf("variable %s is a double array, with %d elements\n", varNames[vi], num);
				double *pr = mxGetPr(arr);
				if (pr != NULL) {
					m_arrayData[m_arrMap[std::string(varNames[vi])]].resize(num);
					m_arrayData[m_arrMap[std::string(varNames[vi])]].assign(pr, pr + num);
				}

				break;

			}

			mxDestroyArray(arr);
		}

		matClose(pmat);
	}

	JointAnglePrior::~JointAnglePrior(){}

	void JointAnglePrior::global2local(const std::vector<double>& S, std::vector<float3>& dSl){

		const std::vector<double>& edges = m_arrayData[m_arrMap["edge"]];
		int numEdge = edges.size() / 2;

		std::vector<float3> dS(numEdge, make_float3(.0, .0, .0));

		for (int i = 0; i < numEdge; i++){
			dS[i].x = S[(edges[i] - 1) * 3] - S[(edges[i + numEdge] - 1) * 3];
			dS[i].y = S[(edges[i] - 1) * 3 + 1] - S[(edges[i + numEdge] - 1) * 3 + 1];
			dS[i].z = S[(edges[i] - 1) * 3 + 2] - S[(edges[i + numEdge] - 1) * 3 + 2];
		}

		const std::vector<double>& chlds = m_arrayData[m_arrMap["chlds"]];
		const std::vector<double>& prnts = m_arrayData[m_arrMap["prnts"]];
		const std::vector<double>& a = m_arrayData[m_arrMap["a"]];
		const std::vector<double>& di = m_arrayData[m_arrMap["di"]];

		int nprts = chlds.size();

		float3 shldr = dS[4] - dS[1];	// std::cout << "(" << shldr.x << " ," << shldr.y << " ," << shldr.z <<  ")" << std::endl;
		float3 hip = dS[12] - dS[8];

		dSl = dS;

		std::vector<float>	R(9, .0);

		for (int i = 0; i < nprts; i++){
			int current_child = chlds[i];		// MATLAB index
			int current_prnt = prnts[i];		// MATLAB index

			float3 u, v, w;

			if (std::find(m_chldsT.begin(), m_chldsT.end(), current_child) != m_chldsT.end()){
				if ((i == 0) || (i == 2) || (i == 4))	u = shldr;
				else u = hip;
				normalize(u);

				v = dS[0];					normalize(v);
			}
			else{

				u = dS[current_prnt - 1];		normalize(u);

				float3 Ra = make_float3((R[0] * a[0] + R[3] * a[1] + R[6] * a[2]),
					(R[1] * a[0] + R[4] * a[1] + R[7] * a[2]),
					(R[2] * a[0] + R[5] * a[1] + R[8] * a[2]));

				float3 di_tmp = make_float3(di[3 * (current_prnt - 1)], di[3 * (current_prnt - 1) + 1], di[3 * (current_prnt - 1) + 2]);

				float3 Rdi = make_float3((R[0] * di_tmp.x + R[3] * di_tmp.y + R[6] * di_tmp.z),
					(R[1] * di_tmp.x + R[4] * di_tmp.y + R[7] * di_tmp.z),
					(R[2] * di_tmp.x + R[5] * di_tmp.y + R[8] * di_tmp.z));

				if ((norm2(u - Ra) < 1e-4) || (norm2(u + Ra) < 1e-4))	v = cross(u, Rdi);
				else				v = cross(Ra, u);

				normalize(v);
			}

			w = cross(u, v);				normalize(w);

			gramschmidt(u, v, w, R);

			dSl[current_child - 1].x = dot(make_float3(R[0], R[1], R[2]), dS[current_child - 1]);
			dSl[current_child - 1].y = dot(make_float3(R[3], R[4], R[5]), dS[current_child - 1]);
			dSl[current_child - 1].z = dot(make_float3(R[6], R[7], R[8]), dS[current_child - 1]);

		}

	}

	void JointAnglePrior::local2global(const std::vector<float3>& dSl, std::vector<float3>& dS){

		const std::vector<double>& edges = m_arrayData[m_arrMap["edge"]];
		int numEdge = edges.size() / 2;

		const std::vector<double>& chlds = m_arrayData[m_arrMap["chlds"]];
		const std::vector<double>& prnts = m_arrayData[m_arrMap["prnts"]];
		const std::vector<double>& a = m_arrayData[m_arrMap["a"]];
		const std::vector<double>& di = m_arrayData[m_arrMap["di"]];

		int nprts = chlds.size();

		float3 shldr = dSl[4] - dSl[1];	// std::cout << "(" << shldr.x << " ," << shldr.y << " ," << shldr.z <<  ")" << std::endl;
		float3 hip = dSl[12] - dSl[8];

		dS[1 - 1] = dSl[1 - 1];
		dS[2 - 1] = dSl[2 - 1];
		dS[5 - 1] = dSl[5 - 1];
		dS[9 - 1] = dSl[9 - 1];
		dS[13- 1] = dSl[13- 1];

		std::vector<float>	R(9, .0);

		for (int i = 0; i < nprts; i++){
			int current_child = chlds[i];		// MATLAB index
			int current_prnt = prnts[i];		// MATLAB index

			float3 u, v, w;

			if (std::find(m_chldsT.begin(), m_chldsT.end(), current_child) != m_chldsT.end()){
				if ((i == 0) || (i == 2) || (i == 4))	u = shldr;
				else u = hip;
				normalize(u);

				v = dS[0];					normalize(v);
			}
			else{

				u = dS[current_prnt - 1];		normalize(u);

				float3 Ra = make_float3((R[0] * a[0] + R[3] * a[1] + R[6] * a[2]),
					(R[1] * a[0] + R[4] * a[1] + R[7] * a[2]),
					(R[2] * a[0] + R[5] * a[1] + R[8] * a[2]));

				float3 di_tmp = make_float3(di[3 * (current_prnt - 1)], di[3 * (current_prnt - 1) + 1], di[3 * (current_prnt - 1) + 2]);

				float3 Rdi = make_float3((R[0] * di_tmp.x + R[3] * di_tmp.y + R[6] * di_tmp.z),
										(R[1] * di_tmp.x + R[4] * di_tmp.y + R[7] * di_tmp.z),
										(R[2] * di_tmp.x + R[5] * di_tmp.y + R[8] * di_tmp.z));

				if ((norm2(u - Ra) < 1e-4) || (norm2(u + Ra) < 1e-4))	v = cross(u, Rdi);
				else				v = cross(Ra, u);

				normalize(v);
			}

			w = cross(u, v);				normalize(w);

			gramschmidt(u, v, w, R);

			dS[current_child - 1].x = dot(make_float3(R[0], R[3], R[6]), dSl[current_child - 1]);
			dS[current_child - 1].y = dot(make_float3(R[1], R[4], R[7]), dSl[current_child - 1]);
			dS[current_child - 1].z = dot(make_float3(R[2], R[5], R[8]), dSl[current_child - 1]);

		}

	}

	void JointAnglePrior::isValid(const std::vector<double>& S, std::vector<bool>& flag){
				

		const std::vector<double>&				edges = m_arrayData[m_arrMap["edges"]];		int numEdges = edges.size() / 2;
		const std::vector<double>&				jmp = m_arrayData[m_arrMap["jmp"]];
		const std::vector<double>&				chlds = m_arrayData[m_arrMap["chlds"]];		int nprts = chlds.size();
		const std::vector<double>&				prnts = m_arrayData[m_arrMap["prnts"]];

		const std::vector<std::vector<double>>&	sepPlane = m_cellData[m_cellMap["sepPlane"]];
		const std::vector<std::vector<double>>&	bounds = m_cellData[m_cellMap["bounds"]];
		const std::vector<std::vector<double>>&	E2 = m_cellData[m_cellMap["E2"]];

		const std::vector<std::vector<bool>>&	angleSprd = m_BcellData[m_BcellMap["angleSprd"]];

		std::vector<float3> dSl(numEdges, make_float3(.0, .0, .0));

		global2local(S, dSl);

		std::vector<std::tuple<int, int>> angle(dSl.size(), std::make_tuple(0, 0));

		for (int i = 0; i < nprts; ++i){
			int current_child = chlds[i];		// MATLAB index
			int current_prnt = prnts[i];		// MATLAB index
			float3 chldB = dSl[current_child - 1];

			float3 chldB_sph = cart2sph(chldB, "DEG");			// r, 90-phi, th	
			chldB_sph.y = 90 - chldB_sph.y;

			chldB /= chldB_sph.x;

			//std::cout << "(" << chldB.x << " ," << chldB.y << " ," << chldB.z <<  ")" << std::endl;

			int t_j = floor((chldB_sph.z + 180) / jmp[0] + 1);
			int p_j = floor((chldB_sph.y + 90) / jmp[0] + 1);

			angle[current_child - 1] = std::make_tuple(t_j, p_j);

			if (std::find(m_chldsT.begin(), m_chldsT.end(), current_child) != m_chldsT.end()){
				int idx = (p_j - 1) * 121 + t_j - 1;
				if (!angleSprd[i][idx]){
					flag[current_child - 1] = false;
					//std::cout << m_BPName[current_child] << " is invalid due to criteria 0." << std::endl;
				}
			}
			else{
				int t_p = std::get<0>(angle[current_prnt - 1]);
				int p_p = std::get<1>(angle[current_prnt - 1]);

				int idx = (p_p - 1) * 121 + t_p - 1;

				float3 v_tmp = make_float3(sepPlane[i][idx], sepPlane[i][idx + 7381], sepPlane[i][idx + 2 * 7381]);	normalize(v_tmp);
				float v4 = sepPlane[i][idx + 3 * 7381] / norm2(v_tmp);

				if (_isnanf(v_tmp.x) || _isnanf(v_tmp.y) || _isnanf(v_tmp.z) || _isnanf(v4) || ((dot(v_tmp, chldB) + v4) > 0)){
					flag[current_child - 1] = false;
					//std::cout << m_BPName[current_child] << " is invalid due to criteria 1." << std::endl;
				}
				else{
					float3 e2 = make_float3(E2[i][idx], E2[i][idx + 7381], E2[i][idx + 2 * 7381]);
					std::vector<float>	T(9, .0);

					gramschmidt(v_tmp, e2, cross(v_tmp, e2), T);

					//std::cout << "(" << T[3] << " ," << T[4] << " ," << T[5] <<  ")" << std::endl;

					double bnd1 = bounds[i][idx];
					double bnd2 = bounds[i][idx + 7381];
					double bnd3 = bounds[i][idx + 2 * 7381];
					double bnd4 = bounds[i][idx + 3 * 7381];

					double u1 = T[3] * chldB.x + T[4] * chldB.y + T[5] * chldB.z;
					double u2 = T[6] * chldB.x + T[7] * chldB.y + T[8] * chldB.z;

					if ((u1 < bnd1) || (u1 > bnd2) || (u2 < bnd3) || (u2 > bnd4)){
						flag[current_child - 1] = false;
						//std::cout << m_BPName[current_child] << " is invalid due to criteria 2." << std::endl;
					}
				}


				/*std::cout << "(" << t_p << " ," << p_p << ")" << std::endl;
				std::cout << "(" << v_tmp.x << " ," << v_tmp.y << " ," << v_tmp.z <<  ")" << std::endl;
				system("pause");*/
			}



		}
	}

	void JointAnglePrior::isValid(const std::vector<double>& S, std::vector<bool>& flag, std::vector<float3>& S2){

		const std::vector<double>&				edges = m_arrayData[m_arrMap["edges"]];		int numEdges = edges.size() / 2;
		const std::vector<double>&				jmp = m_arrayData[m_arrMap["jmp"]];
		const std::vector<double>&				chlds = m_arrayData[m_arrMap["chlds"]];		int nprts = chlds.size();
		const std::vector<double>&				prnts = m_arrayData[m_arrMap["prnts"]];

		const std::vector<std::vector<double>>&	sepPlane = m_cellData[m_cellMap["sepPlane"]];
		const std::vector<std::vector<double>>&	bounds = m_cellData[m_cellMap["bounds"]];
		const std::vector<std::vector<double>>&	E2 = m_cellData[m_cellMap["E2"]];

		const std::vector<std::vector<bool>>&	angleSprd = m_BcellData[m_BcellMap["angleSprd"]];

		std::vector<float3> dSl(numEdges, make_float3(.0, .0, .0));

		global2local(S, dSl);

		std::vector<std::tuple<int, int>> angle(dSl.size(), std::make_tuple(0, 0));

		for (int i = 0; i < nprts; ++i){
			int current_child = chlds[i];		// MATLAB index
			int current_prnt = prnts[i];		// MATLAB index
			float3 chldB = dSl[current_child - 1];

			float3 chldB_sph = cart2sph(chldB, "DEG");			// r, 90-phi, th	
			chldB_sph.y = 90 - chldB_sph.y;

			chldB /= chldB_sph.x;

			//std::cout << "(" << chldB.x << " ," << chldB.y << " ," << chldB.z <<  ")" << std::endl;

			int t_j = floor((chldB_sph.z + 180) / jmp[0] + 1);
			int p_j = floor((chldB_sph.y + 90) / jmp[0] + 1);

			angle[current_child - 1] = std::make_tuple(t_j, p_j);

			if (std::find(m_chldsT.begin(), m_chldsT.end(), current_child) != m_chldsT.end()){
				int idx = (p_j - 1) * 121 + t_j - 1;
				if (!angleSprd[i][idx]){
					flag[current_child - 1] = false;
					//std::cout << m_BPName[current_child] << " is invalid due to criteria 0." << std::endl;
				}
			}
			else{
				int t_p = std::get<0>(angle[current_prnt - 1]);
				int p_p = std::get<1>(angle[current_prnt - 1]);

				int idx = (p_p - 1) * 121 + t_p - 1;

				float3 v_tmp = make_float3(sepPlane[i][idx], sepPlane[i][idx + 7381], sepPlane[i][idx + 2 * 7381]);	normalize(v_tmp);
				float v4 = sepPlane[i][idx + 3 * 7381] / norm2(v_tmp);

				if (_isnanf(v_tmp.x) || _isnanf(v_tmp.y) || _isnanf(v_tmp.z) || _isnanf(v4) || ((dot(v_tmp, chldB) + v4) > 0)){
					flag[current_child - 1] = false;
					//std::cout << m_BPName[current_child] << " is invalid due to criteria 1." << std::endl;
				}
				else{
					float3 e2 = make_float3(E2[i][idx], E2[i][idx + 7381], E2[i][idx + 2 * 7381]);
					std::vector<float>	T(9, .0);

					gramschmidt(v_tmp, e2, cross(v_tmp, e2), T);

					//std::cout << "(" << T[3] << " ," << T[4] << " ," << T[5] <<  ")" << std::endl;

					double bnd1 = bounds[i][idx];
					double bnd2 = bounds[i][idx + 7381];
					double bnd3 = bounds[i][idx + 2 * 7381];
					double bnd4 = bounds[i][idx + 3 * 7381];

					double u1 = T[3] * chldB.x + T[4] * chldB.y + T[5] * chldB.z;
					double u2 = T[6] * chldB.x + T[7] * chldB.y + T[8] * chldB.z;

					if ((u1 < bnd1) || (u1 > bnd2) || (u2 < bnd3) || (u2 > bnd4)){
						flag[current_child - 1] = false;
						//std::cout << m_BPName[current_child] << " is invalid due to criteria 2." << std::endl;
					}
				}
			}
		}

		/*for (size_t ai = 0; ai < angle.size(); ++ai){
			std::cout << std::get<0>(angle[ai]) << ", " << std::get<1>(angle[ai]) << std::endl;
		}*/

		const std::vector<std::vector<double>>&	boundries = m_cellData[m_cellMap["boundries"]];

		std::vector<double> thEdge;// (1 + 360 / jmp[0], .0);
		for (int thi = -180; thi <= 180; thi+=jmp[0]){
			thEdge.push_back(thi);
		}
		assert(thEdge.size() == (1 + 360 / jmp[0]));

		std::vector<double> phEdge;// (1 + 180 / jmp[0], .0);

		for (int phi = -90; phi <= 90; phi += jmp[0]){
			phEdge.push_back(phi);
		}
		assert(phEdge.size() == (1 + 180 / jmp[0]));


		std::vector<float3> dSl2 = dSl;
		std::vector<std::tuple<int, int>> angle2 = angle;

		std::vector<double> epsilon;
		epsilon.push_back(-0.06);
		epsilon.push_back(0.06);
		epsilon.push_back(-0.06);
		epsilon.push_back(0.06);


		for (int i = 0; i < nprts; i++){
			int current_child = chlds[i];		// MATLAB index
			int current_prnt = prnts[i];		// MATLAB index
			float3 chldB = dSl[current_child - 1];		
			float ri = norm2(chldB);		chldB /= ri;


			if (std::find(m_chldsT.begin(), m_chldsT.end(), current_child) != m_chldsT.end()){
				if (~flag[current_child - 1]){
					std::tuple<double, double> angi2;

					findClosestValidPoint(boundries[i], angle2[current_child - 1], angi2);

					angle2[current_child - 1] = angi2;

					double thetas = thEdge[int(std::get<0>(angi2)) - 1] * M_PI / 180;
					double phis = phEdge[int(std::get<1>(angi2)) - 1] * M_PI / 180;

					/*be super careful, we don't know if the input of cos/sin is within the same range*/
					double Xi = ri*cos(phis)*cos(thetas);
					double Yi = ri*cos(phis)*sin(thetas);
					double Zi = ri*sin(phis);
										
					dSl2[current_child - 1] = make_float3(Xi, Yi, Zi);

				}
			}
			else{
				if (~flag[current_child - 1]){
					int t_p = std::get<0>(angle2[current_prnt - 1]);
					int p_p = std::get<1>(angle2[current_prnt - 1]);

					int idx = (p_p - 1) * 121 + t_p - 1;

					float3 v_tmp = make_float3(sepPlane[i][idx], sepPlane[i][idx + 7381], sepPlane[i][idx + 2 * 7381]);	normalize(v_tmp);
					float v4 = sepPlane[i][idx + 3 * 7381] / norm2(v_tmp);

					float3 e2 = make_float3(E2[i][idx], E2[i][idx + 7381], E2[i][idx + 2 * 7381]);
					std::vector<float>	T(9, .0);

					gramschmidt(v_tmp, e2, cross(v_tmp, e2), T);

					//std::cout << "(" << T[3] << " ," << T[4] << " ," << T[5] <<  ")" << std::endl;

					double bnd1 = bounds[i][idx]			- epsilon[0];
					double bnd2 = bounds[i][idx + 7381]		- epsilon[1];
					double bnd3 = bounds[i][idx + 2 * 7381] - epsilon[2];;
					double bnd4 = bounds[i][idx + 3 * 7381] - epsilon[3];;


					if ((dot(v_tmp, chldB) + v4) > 0){
						double u3 = T[3] * chldB.x + T[4] * chldB.y + T[5] * chldB.z;	
						double u4 = T[6] * chldB.x + T[7] * chldB.y + T[8] * chldB.z;	

						chldB.x = T[3] * u3 + T[6] * u4;		chldB.y = T[4] * u3 + T[7] * u4;		chldB.z = T[5] * u3 + T[8] * u4;
					}

					double u1 = T[3] * chldB.x + T[4] * chldB.y + T[5] * chldB.z;		u1 = std::max(std::min(u1, bnd2), bnd1);
					double u2 = T[6] * chldB.x + T[7] * chldB.y + T[8] * chldB.z;		u2 = std::max(std::min(u2, bnd4), bnd3);

					float3 Xh, X2;

					Xh.x = T[3] * u1 + T[6] * u2;		Xh.y = T[4] * u1 + T[7] * u2;		Xh.z = T[5] * u1 + T[8] * u2;

					pointPlane2Sphere(Xh, v_tmp, ri, X2);

					float3 tmp = cart2sph(X2, "DEG");			// r, 90-phi, th	

					tmp.y = 90 - tmp.y;

					int t_j = floor((tmp.z + 180) / jmp[0] + 1);
					int p_j = floor((tmp.y + 90) / jmp[0] + 1);
					idx = (p_j - 1) * 121 + t_j - 1;

					if (((i == 6) || (i == 9)) && (!angleSprd[i][idx])){
						const std::vector<double>& bndry = boundries[i];

						int length_bndry = bndry.size() / 2;

						std::vector<bool>	ind(length_bndry, false);
						std::vector<float3>	bndryPts(length_bndry, make_float3(.0, .0, .0));
						std::vector<double> bndryU;


						for (int bi = 0; bi < length_bndry; bi++){
							double theta = thEdge[bndry[bi * 2] - 1] * M_PI / 180;
							double phi = phEdge[bndry[bi * 2 + 1] - 1] * M_PI / 180;

							double Xi = ri*cos(phi)*cos(theta);
							double Yi = ri*cos(phi)*sin(theta);
							double Zi = ri*sin(phi);

							bndryPts[bi] = make_float3(Xi, Yi, Zi);

							ind[bi] = ((dot(v_tmp, bndryPts[bi]) + v4) < 0);

							if (ind[bi]){
								double u1_tmp = T[3] * bndryPts[bi].x + T[4] * bndryPts[bi].y + T[5] * bndryPts[bi].z;		u1_tmp = std::max(std::min(u1_tmp, bnd2), bnd1);
								double u2_tmp = T[6] * bndryPts[bi].x + T[7] * bndryPts[bi].y + T[8] * bndryPts[bi].z;		u2_tmp = std::max(std::min(u2_tmp, bnd4), bnd3);

								bndryU.push_back(u1_tmp);
								bndryU.push_back(u2_tmp);
							}


						}

						std::tuple<double, double> u2_tuple;

						findClosestValidPoint(bndryU, std::make_tuple(u1, u2), u2_tuple);

						Xh.x = T[3] * std::get<0>(u2_tuple) +T[6] * std::get<1>(u2_tuple);		
						Xh.y = T[4] * std::get<0>(u2_tuple) +T[7] * std::get<1>(u2_tuple);		
						Xh.z = T[5] * std::get<0>(u2_tuple) +T[8] * std::get<1>(u2_tuple);

						pointPlane2Sphere(Xh, v_tmp, ri, X2);

						float3 tmp = cart2sph(X2, "DEG");			// r, 90-phi, th	

						tmp.y = 90 - tmp.y;

						t_j = floor((tmp.z + 180) / jmp[0] + 1);
						p_j = floor((tmp.y + 90) / jmp[0] + 1);
						idx = (p_j - 1) * 121 + t_j - 1;

					}


					dSl2[current_child - 1] = X2;
					angle2[current_child - 1] = std::make_tuple(t_j, p_j);
				}
			}

		}
	
		/*int peek_idx = 2;
		std::cout << dSl2[peek_idx].x << ", " << dSl2[peek_idx].y << ", " << dSl2[peek_idx].z << std::endl;*/

		std::vector<float3> dS2(numEdges, make_float3(.0, .0, .0));
		local2global(dSl2, dS2);
		//std::cout << dS2[peek_idx].x << ", " << dS2[peek_idx].y << ", " << dS2[peek_idx].z << std::endl;

		estimateZ(dS2, edges, make_float3(S[0], S[1], S[2]), S2);
		//std::cout << S2[peek_idx].x << ", " << S2[peek_idx].y << ", " << S2[peek_idx].z << std::endl;
	}


} // end namespace GMorpher

/* ##############################################################################*/
/*put this part in main function when one wants to test sth with JointAnglePrior*/
/*

std::vector<double> S;
if (true)
{
std::map<int, int>	jointMappingIjaz2Paul;

jointMappingIjaz2Paul.insert(std::make_pair(1, 1));			// belly-torso
jointMappingIjaz2Paul.insert(std::make_pair(10, 2));			// RHip-LHip
jointMappingIjaz2Paul.insert(std::make_pair(11, 3));			// RKnee-LKnee
jointMappingIjaz2Paul.insert(std::make_pair(12, 4));			// RAnkle-LFeet
jointMappingIjaz2Paul.insert(std::make_pair(14, 5));			// LHip-RHip
jointMappingIjaz2Paul.insert(std::make_pair(15, 6));			// LKnee-RKnee
jointMappingIjaz2Paul.insert(std::make_pair(16, 7));			// LAnkle-RFeet
jointMappingIjaz2Paul.insert(std::make_pair(2, 8));			// Neck-Neck
jointMappingIjaz2Paul.insert(std::make_pair(9, 9));			// Face-Head
jointMappingIjaz2Paul.insert(std::make_pair(3, 10));			// RShldr-LShldr
jointMappingIjaz2Paul.insert(std::make_pair(4, 11));			// RElbow-LElbow
jointMappingIjaz2Paul.insert(std::make_pair(5, 12));			// RWrist-LWrist
jointMappingIjaz2Paul.insert(std::make_pair(6, 13));			// LShldr-RShldr
jointMappingIjaz2Paul.insert(std::make_pair(7, 14));			// LElbow-RElbow
jointMappingIjaz2Paul.insert(std::make_pair(8, 15));			// LWrist-RWrist

std::ifstream f_in("C:\\Users\\Paul\\Documents\\Visual Studio 2010\\Projects\\JointPS_forest\\bin\\Release\\output\\ghandstand\\joints\\joints_0403");

std::vector<float3> JX(15, make_float3(.0, .0, .0));
for (size_t ji = 0; ji < 15; ji++){
f_in >> JX[ji].x >> JX[ji].z >> JX[ji].y;
}
f_in.close();		//		std::cout << "loading joints successfully" << std::endl;
S.resize(17*3);

for (size_t ji = 0; ji < 17; ji++){
if ((ji == 16) || (ji==12))
{
S[ji * 3 + 0] = std::numeric_limits<double>::quiet_NaN();
S[ji * 3 + 1] = std::numeric_limits<double>::quiet_NaN();
S[ji * 3 + 2] = std::numeric_limits<double>::quiet_NaN();
}
else{
S[ji * 3 + 0] = JX[jointMappingIjaz2Paul[ji + 1] - 1].x;
S[ji * 3 + 1] = JX[jointMappingIjaz2Paul[ji + 1] - 1].y;
S[ji * 3 + 2] = JX[jointMappingIjaz2Paul[ji + 1] - 1].z;
}
}		//		std::cout << "mapping joints successfully" << std::endl;
}
else{
string testPose = "C:\\Users\\Paul\\Desktop\\Ijaz\\pose-conditioned-prior3\\testPose.mat";
MATFile *pmat = matOpen(testPose.c_str(), "r");
mxArray *pa;

if (pmat == NULL) { printf("Error opening file %s\n", testPose.c_str()); return(EXIT_FAILURE); }
pa = matGetVariable(pmat, "Sd");
if (pa == NULL) {
printf("Error reading existing matrix LocalDouble\n");
return(EXIT_FAILURE);
}

if (matClose(pmat) != 0) {
printf("Error closing file %s\n", testPose.c_str());
return(EXIT_FAILURE);
}

mwSize num_ele = (mwSize)mxGetNumberOfElements(pa);
double *pr = mxGetPr(pa);

S.resize(num_ele);
S.assign(pr, pr + num_ele);

}
// prior model
string priorName = "C:\\Users\\Paul\\Desktop\\Ijaz\\pose-conditioned-prior3\\jointAngleModel.mat";
GMorpher::JointAnglePrior prior(priorName.c_str());

std::vector<bool>	flag(16, true);
std::vector<float3> S2(17, make_float3(.0, .0, .0));
prior.isValid(S, flag, S2);
//prior.isValid(S, flag);

for (size_t i = 0; i < flag.size(); i++){
std::cout << flag[i] << ", ";
}
std::cout << std::endl;

return 0;

*/