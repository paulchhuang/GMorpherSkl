/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef CC3D_FLOAT3_H_DEFINED
#define CC3D_FLOAT3_H_DEFINED
#pragma warning(disable: 4267) 
#pragma warning(disable: 4244) 
#define _USE_MATH_DEFINES
#include <math.h> // for sqrtf
#include <iostream>

namespace CC3D
{

struct float3     { float x,y,z; };

// ####################
// CONSTRUCTION
inline float3 make_float3( const float& x, const float& y, const float& z ) {
	float3 res;
	res.x = x; res.y = y; res.z = z;
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


// ####################
// MATRIX ALGEBRA
inline float3 prod3( const float* A, const float3& b )
{
	return make_float3( A[0*3+0] * b.x +  A[0*3+1] * b.y + A[0*3+2] * b.z,
	                    A[1*3+0] * b.x +  A[1*3+1] * b.y + A[1*3+2] * b.z,
	                    A[2*3+0] * b.x +  A[2*3+1] * b.y + A[2*3+2] * b.z);
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

inline float3& operator *= ( float3&a, const float s) {
	a.x *= s;
	a.y *= s;
	a.z *= s;
	return a;
}

inline float3& operator /= ( float3&a, const float s) {
	a.x /= s;
	a.y /= s;
	a.z /= s;
	return a;
}

inline float dot( const float3& a, const float3& b ) {
	return a.x*b.x + a.y*b.y + a.z*b.z;
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

inline float norm_squared( const float3& a ) {
	return dot(a,a);
}

inline void normalize( float3& a ) {
	float norm = norm2(a);
	a.x /= norm;
	a.y /= norm;
	a.z /= norm;
}

// ###################
// cart2sph
inline float3 cart2sph(const float3& cart, const std::string& degree_flag){
	float r = sqrtf(dot(cart,cart));
	float theta;
	float phi;

	/******************************//**
	* \brief conversion from cartesian to spherical coordinates (if the
	* resulting radius is 0 also the azimutal and inclination
	* angles get set to 0).
	* theta(y): [0 180]
	* phi(z):[-180 180]	
	******************************/

	if(r!=0){
		theta = (!degree_flag.compare("DEG")) ? (acos( cart.z / r )*180/M_PI):acos( cart.z / r );
		phi =	(!degree_flag.compare("DEG")) ? (atan2( cart.y, cart.x )*180/M_PI):atan2( cart.y, cart.x );
	}
	else{
		theta = phi = 0.;
	}

	return make_float3(r,theta,phi);
}


}








#endif
