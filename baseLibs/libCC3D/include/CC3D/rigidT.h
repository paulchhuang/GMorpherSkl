/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef CC3D_RIGIDT_H_DEFINED
#define CC3D_RIGIDT_H_DEFINED


#include "quaternions.h"
#include <vector>

namespace CC3D 
{


struct rigidT {
	quaternion m_q;
	float3     m_t;
};


inline void RigidTFlatten( const rigidT& RT, float* RTf)
{
	Quat_to_Rot<float>( (const float*)( &(RT.m_q)), RTf );
	RTf[9]  = RT.m_t.x;
	RTf[10] = RT.m_t.y;
	RTf[11] = RT.m_t.z;
}

inline void RigidTFlatten( const std::vector<rigidT>& RT, float* RTf ){
	int numPatches = RT.size();

	#pragma omp parallel for
	for(int pi=0;pi<numPatches;++pi) {
		RigidTFlatten( RT[pi], RTf+12*pi);
		//RTf += 12;
	}
}


//inline RigidT compose( const RigidT& RT1, const RigidT& RT2 ) {
//	RigidT res;
//	float R[9];
//	Quat_multiply<float>( (const float*)&(RT1.m_q), (const float*)&(RT2.m_q), (float*)&(res.m_q));
//	Quat_to_Rot<float>( (const float*)&(RT1.m_q), R);
//	res.m_t = prod3( R, RT2.m_t) + RT1.m_t;

//	return res;
//}


}



#endif
