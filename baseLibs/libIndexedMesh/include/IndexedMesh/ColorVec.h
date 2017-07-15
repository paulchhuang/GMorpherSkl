/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef INDEXEDMESH_COLORVEC_H_DEFINED
#define INDEXEDMESH_COLORVEC_H_DEFINED



namespace IndexedMesh3D 
{

struct ColorVec { float r,g,b,a; };


inline ColorVec make_ColorVec( const float& r, const float& g, 
                               const float& b, const float& a ) {
	ColorVec res;
	res.r = r; res.g = g; res.b = b;  res.a = a;
	return res;
}

}





#endif
