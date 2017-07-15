/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef RIGIDT_H_DEFINED
#define RIGIDT_H_DEFINED

//#include "Quaternions.h"
#include <cmath>
#include <vector>
#include <iostream>

#include <CC3D/float3.h>
#include <CC3D/rigidT.h>


namespace GMorpher
{

typedef CC3D::float3 float3;
typedef CC3D::rigidT rigidT;

void RigidTUpdate( const std::vector<rigidT>& RT_old, const double* X, 
                   const double scale, const std::vector<int>& patch_comp, std::vector<rigidT>& RT);


/*
struct float3     { float x,y,z; };
struct quaternion { float w,x,y,z; };

struct RigidT {
	quaternion m_q;
	float3     m_t;
};


void RigidTFlatten( const RigidT& RT, float* RTf);
void RigidTFlatten( const std::vector<RigidT>& RT, float* RTf );


inline float3 make_float3( const float& x, const float& y, const float& z ) {
	float3 res;
	res.x = x; res.y = y; res.z = z;
	return res;
}

inline quaternion make_quaternion( const float& w, const float& x, const float& y, const float& z ) {
	quaternion res;
	res.w = w; res.x = x; res.y = y; res.z = z;
	return res;
}


// ####################
// VECTOR GROUP
inline float3 operator - ( const float3& a ) {
	return make_float3(-a.x, -a.y, -a.z);
}

inline float3 operator + ( const float3& a, const float3& b ) {
	return make_float3( a.x+b.x, a.y+b.y, a.z+b.z );
}

inline float3 operator - ( const float3& a, const float3& b ) {
	return make_float3( a.x-b.x, a.y-b.y, a.z-b.z) ;
}

inline float3& operator += (float3& a, const float3& b) {
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
	return a;
}

inline float3& operator -= (float3& a, const float3& b) {
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
	return a;
}

inline float3 prod3( const float* A, const float3& b )
{
	return make_float3( A[0*3+0] * b.x +  A[0*3+1] * b.y + A[0*3+2] * b.z,
	                    A[1*3+0] * b.x +  A[1*3+1] * b.y + A[1*3+2] * b.z,
	                    A[2*3+0] * b.x +  A[2*3+1] * b.y + A[2*3+2] * b.z);
}

inline float3 prodtranspose3( const float* A, const float3& b )
{
	return make_float3( A[0*3+0] * b.x +  A[1*3+0] * b.y + A[2*3+0] * b.z,
	                    A[0*3+1] * b.x +  A[1*3+1] * b.y + A[2*3+1] * b.z,
	                    A[0*3+2] * b.x +  A[1*3+2] * b.y + A[2*3+2] * b.z);
}

// ####################
// VECTOR SPACE
inline float3 operator * ( const float3& a, const float& s ) {
	return make_float3( a.x*s, a.y*s, a.z*s );
}

inline float3 operator * ( const float& s, const float3& a ) {
	return make_float3( a.x*s, a.y*s, a.z*s );
}

inline float3 operator / ( const float3& a, const float& s ) {
	return make_float3( a.x/s, a.y/s, a.z/s );
}

inline float dot( const float3& a, const float3& b ) {
	return a.x*b.x + a.y*b.y+a.z*b.z;
}

// ####################
// VECTOR ALGEBRA
inline float3 cross ( const float3& a, const float3& b ) {
	return make_float3( a.y*b.z - a.z*b.y,
	                    a.z*b.x - a.x*b.z,
	                    a.x*b.y - a.y*b.x );
}


// ####################
// Normalization
inline float norm2( const float3& a ) {
	return sqrtf( dot(a,a) );
}

inline void normalize( float3& a ) {
	float norm = norm2(a);
	a.x /= norm;
	a.y /= norm;
	a.z /= norm;
}


// ####################
// RigidT Compose
inline RigidT compose( const RigidT& RT1, const RigidT& RT2 ) {
	RigidT res;
	float R[9];
	Quat_multiply<float>( (const float*)&(RT1.m_q), (const float*)&(RT2.m_q), (float*)&(res.m_q));
	Quat_to_Rot<float>( (const float*)&(RT1.m_q), R);
	res.m_t = prod3( R, RT2.m_t) + RT1.m_t;

	return res;
}


// ####################
// IO
inline std::ostream& operator << ( std::ostream& os, const float3& a ) {
	os << "["<< a.x<<", "<<a.y<<", "<<a.z<<"]";
	return os;
}

inline void printMat3( std::ostream& os, const float* A ) {
	os<<"[";
	for(int i=0;i<3;++i) {
		os<<"\n[";
		for(int j=0;j<3;++j) {
			os<<A[3*i+j]<<", ";
		}
		os<<"\b\b],";
	}
	os<<"\b]\n";
}

inline std::istream& operator >> ( std::istream&is, RigidT& RT ) {
	is >> RT.m_q.w >> RT.m_q.x >> RT.m_q.y >> RT.m_q.z
	   >> RT.m_t.x >> RT.m_t.y >> RT.m_t.z;
	return is;
}

inline std::ostream& operator << ( std::ostream&os, const RigidT& RT ) {
	os << RT.m_q.w <<" "<< RT.m_q.x <<" "<< RT.m_q.y <<" "<< RT.m_q.z
	   <<" "<< RT.m_t.x <<" "<< RT.m_t.y <<" "<< RT.m_t.z<<"\n";

	return os;
}
*/

}









#endif
