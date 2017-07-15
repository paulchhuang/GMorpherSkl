/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef QUATERNIONS_H_DEFINED
#define QUATERNIONS_H_DEFINED

#include <cmath>

namespace GMorpher
{

// q is wxyz
template <typename T>
inline void Quat_to_Rot( const T* q, T* R)
{
	R[0*3+0] = 1.  - (2*q[2]*q[2] + 2*q[3]*q[3]); // 1 -2yy - 2zz
	R[1*3+0] =        2*q[1]*q[2] + 2*q[3]*q[0];  //    2xy + 2zw
	R[2*3+0] =        2*q[1]*q[3] - 2*q[2]*q[0];  //    2xz - 2yw

	R[0*3+1] =        2*q[1]*q[2] - 2*q[3]*q[0];  //    2xy - 2zw
	R[1*3+1] = 1.  - (2*q[1]*q[1] + 2*q[3]*q[3]); // 1 - 2xx - 2zz
	R[2*3+1] =        2*q[2]*q[3] + 2*q[1]*q[0];  //     2yz + 2xww

	R[0*3+2] =        2*q[1]*q[3] + 2*q[2]*q[0];  //    2xz + 2yw
	R[1*3+2] =        2*q[2]*q[3] - 2*q[1]*q[0];  //    2yz - 2xw
	R[2*3+2] = 1.  - (2*q[1]*q[1] + 2*q[2]*q[2]); // 1 - 2xx - 2yy
}

template <typename T>
inline void Quat_from_rot( const T* R, T* q)
{
	T trace = 1. + R[0*3+0] + R[1*3+1] + R[2*3+2];

	if ( trace > 0.00000001 ) {
		double s = 0.5 / sqrt(trace);
		q[0] = 0.25 / s;
		q[1] = ( R[2*3+1] - R[1*3+2] ) * s;
		q[2] = ( R[0*3+2] - R[2*3+0] ) * s;
		q[3] = ( R[1*3+0] - R[0*3+1] ) * s;
	}

	else {
		if ( R[0*3+0] > R[1*3+1] && R[0*3+0] > R[2*3+2] ) {
			double s = 2.0 * sqrt( 1.0 + R[0*3+0] - R[1*3+1] - R[2*3+2]);
			q[0] = (R[2*3+1] - R[1*3+2] ) / s;
			q[1] = 0.25 * s;
			q[2] = (R[0*3+1] + R[1*3+0] ) / s;
			q[3] = (R[0*3+2] + R[2*3+0] ) / s;
		} else if (R[1*3+1] > R[2*3+2]) {
			double s = 2.0 * sqrt( 1.0 + R[1*3+1] - R[0*3+0] - R[2*3+2]);
			q[0] = (R[0*3+2] - R[2*3+0] ) / s;
			q[1] = (R[0*3+1] + R[1*3+0] ) / s;
			q[2] = 0.25 * s;
			q[3] = (R[1*3+2] + R[2*3+1] ) / s;
		} else {
			double s = 2.0 * sqrt( 1.0 + R[2*3+2] - R[0*3+0] - R[1*3+1] );
			q[0] = (R[1*3+0] - R[0*3+1] ) / s;
			q[1] = (R[0*3+2] + R[2*3+0] ) / s;
			q[2] = (R[1*3+2] + R[2*3+1] ) / s;
			q[3] = 0.25 * s;
		}
	}
}

template <typename T>
inline void Quat_normalise( T* q)
{
	T norm = sqrt( q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] );
	q[0] /= norm;
	q[1] /= norm;
	q[2] /= norm;
	q[3] /= norm;
}

template <typename T>
inline void Quat_multiply( const T* q1, const T* q2, T* q)
{
	// w1.w2 - v1.v2
	q[0] = q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3];
	// w1.w2 - v1.v2
	q[1] = q1[0]*q2[1] + q2[0]*q1[1] + ( q1[2]*q2[3] - q1[3]*q2[2] );
	q[2] = q1[0]*q2[2] + q2[0]*q1[2] + ( q1[3]*q2[1] - q1[1]*q2[3] );
	q[3] = q1[0]*q2[3] + q2[0]*q1[3] + ( q1[1]*q2[2] - q1[2]*q2[1] );
}

} // end namespace GMorpher

#endif
