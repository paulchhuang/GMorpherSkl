/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#include <RigidT.h>
#include <Quaternions.h>





namespace GMorpher
{

void RigidTUpdate( const std::vector<rigidT>& RT_old, const double* X, 
                   const double scale, const std::vector<int>& patch_comp, 
				   std::vector<rigidT>& RT)
{
	int numPatches = RT_old.size();
	RT.resize(numPatches);

	#pragma omp parallel for
	for(int pi=0;pi<numPatches;++pi )
	{
		float w[3];
		float v[3];
		for(int i=0;i<3;++i) w[i] = X[6*pi+ 3*patch_comp[pi]*15 +i]   * scale;
		for(int i=0;i<3;++i) v[i] = X[6*pi+ 3*patch_comp[pi]*15 +3+i] * scale;

		// ------------------
		// do quaternion magic
		float theta = sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
		float qupdate[4];
		qupdate[0] = cos(theta/2.);
		float scaler = 0.;
		if( fabs(theta) < 0.000001 ) scaler = 0.5;
		else scaler =  sin(theta/2.) / theta;
		qupdate[1] = scaler * w[0];
		qupdate[2] = scaler * w[1];
		qupdate[3] = scaler * w[2];

		// ----------------------
		// write output
		const rigidT& RTi_old = RT_old[pi];
		rigidT RTi_new;
		Quat_multiply<float>( qupdate, (const float*)( &(RTi_old.m_q) ) , (float*)( &(RTi_new.m_q) ) );
		Quat_normalise<float>((float*)( &(RTi_new.m_q)) );
		RTi_new.m_t.x = RTi_old.m_t.x + v[0];
		RTi_new.m_t.y = RTi_old.m_t.y + v[1];
		RTi_new.m_t.z = RTi_old.m_t.z + v[2];

		RT[pi] = RTi_new;
	}
}


} // end namespace GMorpher
