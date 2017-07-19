/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#include <MVLib.h>
#include <fstream>
#include <iostream>

namespace MVLib
{



	void loadCameras(const char*            cameraCalibBaseName,
	                 const std::list<int>&  indexList,
	                 std::map<int, Camera>& cams, const bool KinovisFlag/*=false*/)
	{
		cams.clear();

		if (!KinovisFlag)
		{
			for (std::list<int>::const_iterator i_itr = indexList.begin(); i_itr != indexList.end(); ++i_itr)
			{
				Matrix3x4 P = loadPMatrix(buildFilename(cameraCalibBaseName, *i_itr).c_str());
				Camera cam(P);
				cams.insert(std::pair<int, Camera>(*i_itr, cam));
			}
		}
		else
		{
			std::ifstream fin(cameraCalibBaseName);
			if (!fin.is_open()) throw std::runtime_error(std::string("could not open :") + std::string(cameraCalibBaseName));
			for (std::list<int>::const_iterator i_itr = indexList.begin(); i_itr != indexList.end(); ++i_itr)
			{
				double tempvec[12];
				for (int i = 0; i<12; ++i) fin >> tempvec[i];
				//Matrix3x4 res(tempvec, tempvec+12);
				Matrix3x4 P;
				P << tempvec[0], tempvec[1], tempvec[2], tempvec[3],
					tempvec[4], tempvec[5], tempvec[6], tempvec[7],
					tempvec[8], tempvec[9], tempvec[10], tempvec[11];
				
				Camera cam(P);
				cams.insert(std::pair<int, Camera>(*i_itr, cam));
			}

			fin.close();

		}
	}


	void getNearFar(const std::map<int, Camera>& cams,
	                const double xmin, const double ymin, const double zmin,
	                const double xmax, const double ymax, const double zmax,
	                double& znear, double& zfar)
	{
		znear = std::numeric_limits<double>::max();
		zfar = std::numeric_limits<double>::min();

		double corners[3*8];
		corners[0*3+0] = xmin;corners[0*3+1] = ymin;corners[0*3+2] = zmin;
		corners[1*3+0] = xmax;corners[1*3+1] = ymin;corners[1*3+2] = zmin;
		corners[2*3+0] = xmax;corners[2*3+1] = ymax;corners[2*3+2] = zmin;
		corners[3*3+0] = xmin;corners[3*3+1] = ymax;corners[3*3+2] = zmin;
		corners[4*3+0] = xmin;corners[4*3+1] = ymin;corners[4*3+2] = zmax;
		corners[5*3+0] = xmax;corners[5*3+1] = ymin;corners[5*3+2] = zmax;
		corners[6*3+0] = xmax;corners[6*3+1] = ymax;corners[6*3+2] = zmax;
		corners[7*3+0] = xmin;corners[7*3+1] = ymax;corners[7*3+2] = zmax;

		for( std::map<int,Camera>::const_iterator cam_itr = cams.begin(); cam_itr != cams.end(); ++ cam_itr)
		{
			const Matrix3& Rw2c = cam_itr->second.Rw2c();
			const Vector3& Pos  = cam_itr->second.Pos();
			Vector3 pointPosW;
			for (int i=0;i<8;++i)
			{
				/*pointPosW = corners[i*3+0], corners[i*3+1], corners[i*3+2];
				Vector3 pointPosCam(prod(Rw2c, pointPosW - Pos));*/
				pointPosW << corners[i*3+0], corners[i*3+1], corners[i*3+2];
				Vector3 pointPosCam(Rw2c*(pointPosW - Pos));

				znear = std::min(znear, pointPosCam(2));
				zfar = std::max(zfar, pointPosCam(2));
			}
		}
	}

}
