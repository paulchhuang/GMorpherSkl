/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef ICP_FUNCTION_H_DEFINED
#define ICP_FUNCTION_H_DEFINED

#include <PatchedCloud.h>


inline void accumEnergy( const double  w,
                         const SklRgedPatchedCloud::float3& xo,
                         const SklRgedPatchedCloud::float3& no,
                         const SklRgedPatchedCloud::float3& X,   // the transformed version
                         const SklRgedPatchedCloud::float3& N,   // the rotated version
                         double&       E )
{
	SklRgedPatchedCloud::float3 d;
	d.x = X.x - xo.x;
	d.y = X.y - xo.y;
	d.z = X.z - xo.z;
	E += w *( d.x*d.x + d.y*d.y + d.z*d.z);
}


inline void accumJacobian( const double  w,
                           const SklRgedPatchedCloud::float3& xo,
                           const SklRgedPatchedCloud::float3& no,
                           const SklRgedPatchedCloud::float3& dX,  // the rotated version
                           const SklRgedPatchedCloud::float3& X,   // the transformed version
                           const SklRgedPatchedCloud::float3& N,   // the rotated version
                           double*       JTJ,               // 6x6
                           double*       JTb )              // 6x1
{
	SklRgedPatchedCloud::float3 d;
	d.x = X.x - xo.x;
	d.y = X.y - xo.y;
	d.z = X.z - xo.z;

	SklRgedPatchedCloud::float3 wdX;
	wdX.x = w * dX.x;
	wdX.y = w * dX.y;
	wdX.z = w * dX.z;


	// top left 3x3
	JTJ[0*6+0] += wdX.y*dX.y + wdX.z*dX.z;    JTJ[0*6+1] -= wdX.x*dX.y;               JTJ[0*6+2] -= wdX.x*dX.z;
	/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/  JTJ[1*6+1] += wdX.x*dX.x + wdX.z*dX.z;  JTJ[1*6+2] -= wdX.y*dX.z;
	/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/ /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/ JTJ[2*6+2] += wdX.x*dX.x + wdX.y*dX.y;
	// top right 3x3
	/*xxxxxxxxxxxxxxxxxxx*/ JTJ[0*6+4] -= wdX.z;   JTJ[0*6+5] += wdX.y;
	JTJ[1*6+3] += wdX.z;   /*xxxxxxxxxxxxxxxxxxx*/ JTJ[1*6+5] -= wdX.x;
	JTJ[2*6+3] -= wdX.y;   JTJ[2*6+4] += wdX.x;   /*xxxxxxxxxxxxxxxxxxx*/
	// bottom right 3x3
	JTJ[3*6+3] += w;
	JTJ[4*6+4] += w;
	JTJ[5*6+5] += w;




	JTb[0] -= w*(- dX.z*d.y + dX.y*d.z ); // -dXz * y + dXy * z
	JTb[1] -= w*(  dX.z*d.x - dX.x*d.z);  //  dXz * x - dXx * z
	JTb[2] -= w*(- dX.y*d.x + dX.x*d.y ); // -dXy * x + dXx * y
	JTb[3] -= w*(d.x);
	JTb[4] -= w*(d.y);
	JTb[5] -= w*(d.z);
}





//inline void accumEnergy( const double  w,
//                         const PatchedCloud::float3& xo,
//                         const PatchedCloud::float3& no,
//                         const PatchedCloud::float3& X,   // the transformed version
//                         const PatchedCloud::float3& N,   // the rotated version
//                         double&       E )
//{
//	PatchedCloud::float3 Ns    = no + N;
//	PatchedCloud::float3 delta = X - xo;
//	double poinoPlane = dot(Ns,delta);
//	E += w*(poinoPlane*poinoPlane);
//}




//inline void accumJacobian( const double  w,
//                           const PatchedCloud::float3& xo,
//                           const PatchedCloud::float3& no,
//                           const PatchedCloud::float3& dX,  // the rotated version
//                           const PatchedCloud::float3& X,   // the transformed version
//                           const PatchedCloud::float3& N,   // the rotated version
//                           double*       JTJn,               // 6x6
//                           double*       JTbn )              // 6x1
//{
//	PatchedCloud::float3 Ns    = no + N;
//	PatchedCloud::float3 delta = X - xo;
//	double poinoPlane = dot(Ns,delta);

//	/*
//	f  = w N.( X2 - X1 )
//	   = w N.( [w] dX2 + v + c - X1 )
//	gf = - w N [dX2] , w N
//	     (1x3)x(3x3), (1x3)
//	*/
//	double gn[6]; // the normalized gradieno
//	gn[0] = (dX.y * Ns.z - dX.z * Ns.y); // D poinoPlane / DW1
//	gn[1] = (dX.z * Ns.x - dX.x * Ns.z); // D poinoPlane / DW2
//	gn[2] = (dX.x * Ns.y - dX.y * Ns.x); // D poinoPlane / DW3
//	gn[3] = (Ns.x);
//	gn[4] = (Ns.y);
//	gn[5] = (Ns.z);

//	// TODO remove the bottom right part
//	for(ino i=0;i<6;++i ) {
//		for(ino j=0;j<6;++j) JTJn[6*i+j] += w * gn[i] * gn[j];
//	}
//	for( ino i=0;i<6;++i ) JTbn[i] -= w * poinoPlane * gn[i];
//}



#endif
