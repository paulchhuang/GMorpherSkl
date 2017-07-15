/* *************************************************
 * Copyright (2017) : Paul Huang, Cedric Cagniart
 * *************************************************/
#pragma warning(disable: 4267) 
#pragma warning(disable: 4244) 
#include <IndexedMesh/IndexedMesh3D.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

#include <boost\foreach.hpp>
#include <cassert>
#include <numeric>  // for std::accumulate
#include <queue>

#define VERBOSE false



namespace IndexedMesh3D
{

	// ##################################################################
	// AUX IO FUNCTIONS
	// ##################################################################
	void loadTriangles(std::istream& s,
	                   std::vector<Triangle>& triangles)
	{
		int numTriangles = triangles.size();
		
		BOOST_FOREACH(Triangle& t_i, triangles){
			int faceVertexCount;
			s >> faceVertexCount;
			if(faceVertexCount != 3) throw(std::ios_base::failure("this is not a triangular mesh"));
			s >> t_i.v0 >> t_i.v1>> t_i.v2;
		}

	}

	void loadCoords(std::istream& s,
	                std::vector<float3>& coords)
	{
		int numVertices = coords.size();

		BOOST_FOREACH(float3& coords_i, coords){			
			s >> coords_i.x >> coords_i.y >> coords_i.z;
		}
	}

	void loadCoordsColors(std::istream& s,
	                      std::vector<float3>&   coords,
	                      std::vector<ColorVec>& colors)
	{
		int numVertices = coords.size();
		for(int vi = 0; vi < numVertices; ++vi)
		{
			float3  & p = coords[vi];
			ColorVec& c = colors[vi];
			s >> p.x >> p.y >> p.z >> c.r >> c.g >> c.b >> c.a;
		}
	}

	void loadCoordsNormals( std::istream& s,
	                        std::vector<float3>& coords )
	{
		int numVertices = coords.size();
		double a;
		for(int vi = 0; vi < numVertices; ++vi)
		{
			float3& p = coords[vi];
			s >> p.x >> p.y >> p.z >> a >> a >> a;
		}
	}

	void parseObj(	const std::string& meshFile,
					std::vector<float3>&   coords,
					std::vector<ColorVec>& colors,
					std::vector<Triangle>& triangles)
	{
		std::ifstream file_in(meshFile);	//this line reads the .off file
		if (!file_in.is_open())throw(std::ios_base::failure(meshFile));

		char line[1024];
		std::string op;
		std::istringstream line_in;

		while (file_in.good()){
			file_in.getline(line, 1023);
			line_in.clear();
			line_in.str(line);

			if (!(line_in >> op))
				continue;
			if (op == "v"){
				float3 point_i;
				line_in >> point_i.x >> point_i.y >> point_i.z;
				coords.push_back(point_i);
			}
			else if (op == "vt"){
				/*To do: parse texture map*/
			}
				
			else if (op == "vn"){
				/*To do: parse normals*/
			}
				
			else if (op == "g"){
				/*groups.clear();
				while (line_in >> groups);
				groups.insert("default");*/
			}
			else if (op == "f") {
				Triangle tri_i;
				line_in >> tri_i.v0 >> tri_i.v1 >> tri_i.v2;
				triangles.push_back(tri_i);
			}
		}
		file_in.close();
	}





	// ##################################################################
	// IO
	// ##################################################################
	void loadIndexedMesh3D( const std::string& meshFile,
	                        std::vector<float3>&   coords,
	                        std::vector<ColorVec>& colors,
	                        std::vector<Triangle>& triangles)
	{
		std::size_t found = meshFile.find_last_of(".");
		if ((meshFile.substr(found + 1) == "obj") || (meshFile.substr(found + 1) == "OBJ")){
			parseObj(meshFile, coords, colors, triangles);
			return;
		}


		enum OFFType {TYPE_OFF, TYPE_COFF, TYPE_NOFF};

		OFFType offtype;
		int nVertices, nFaces, nEdges;
		std::string buffer;
		std::cout << meshFile << std::endl;
		std::ifstream file_in(meshFile);	//this line reads the .off file
		if (!file_in.is_open())throw(std::ios_base::failure(meshFile));

		std::getline(file_in,buffer);	//reads the first line of .off file

		if (buffer.find("OFF") != std::string::npos)  { offtype = TYPE_OFF; }
		else throw(std::ios_base::failure("not an off file"));
		if (buffer.find("COFF") != std::string::npos) { offtype = TYPE_COFF; }
		if (buffer.find("NOFF") != std::string::npos) { offtype = TYPE_NOFF; }


		file_in >> nVertices >> nFaces >> nEdges;

		coords.resize(nVertices);
		colors.resize(nVertices);
		triangles.resize(nFaces);
		
		switch(offtype)
		{
			case TYPE_OFF :
			{				
				loadCoords(file_in, coords);
				ColorVec zerocolor = make_ColorVec(0,0,0,255);
				std::fill(colors.begin(), colors.end(), zerocolor);		// initialize colors with all 0.
				break;
			}
			case TYPE_COFF :
			{
				loadCoordsColors(file_in, coords, colors);
				break;
			}
			case TYPE_NOFF :
			{
				loadCoordsNormals(file_in, coords);
				break;
			}
			default :
			{
				throw(std::ios_base::failure("could not infer the subtype of off file"));
				break;
			}
		}	
		
		loadTriangles(file_in, triangles);
		file_in.close();		
		if (VERBOSE) std::cout << "-- MESHIO : loaded a mesh with verts/tris : " << nVertices << "/" << nFaces << std::endl;
	}




	void saveIndexedMesh3D( const char* filename,
	                        const std::vector<float3>&   coords,
	                        const std::vector<ColorVec>& colors,
	                        const std::vector<Triangle>& triangles)
	{
		int numVertices = coords.size();
		int numTriangles = triangles.size();

		assert(int(colors.size()) == numVertices);

		std::ofstream f_out(filename);
		if(!f_out.is_open())throw(std::ios_base::failure(filename));

		f_out<<"COFF\n"<<numVertices<<" "<<numTriangles<<" 0"<<std::endl;

		// write down vertices
		for (int vi=0;vi<numVertices;++vi)
		{
			const float3&   p = coords[vi];
			const ColorVec& c = colors[vi];
			f_out<<p.x<<" "<<p.y<<" "<<p.z<<" "<<c.r<<" "<<c.g<<" "<<c.b<<" "<<c.a<<"\n";
		}
		// write down triangles

		BOOST_FOREACH(const Triangle& t, triangles){
			f_out<<"3 "<<t.v0<<" "<<t.v1<<" "<<t.v2<<"\n";
		}
		/*for(int  ti=0;ti<numTriangles;++ti)
		{
			const Triangle& t = triangles[ti];
			f_out<<"3 "<<t.v0<<" "<<t.v1<<" "<<t.v2<<"\n";
		}*/

		f_out.close();

		if(VERBOSE) std::cout<<"-- MESHIO : saved a mesh with verts/tris : "<<numVertices<<"/"<<numTriangles<<std::endl;
	}






	// ##################################################################
	// count triangles
	// ##################################################################
	int getNumTriangles(const std::vector<Edge>& CWEdges,
	                    const std::vector<int>&  CWEdgesBoundaries)
	{
		int numVertices = CWEdgesBoundaries.size() - 1;
		int numTriangles = 0;
		for(int vi = 0; vi< numVertices; ++vi)
		{
			std::vector<Edge>::const_iterator eBegin = CWEdges.begin() + CWEdgesBoundaries[vi];
			std::vector<Edge>::const_iterator eEnd = CWEdges.begin() + CWEdgesBoundaries[vi+1];

			for(std::vector<Edge>::const_iterator e_itr = eBegin; e_itr != eEnd; ++e_itr) {
				if( e_itr->hasNext && e_itr->v < vi && e_itr->vnext < vi) numTriangles++;
			}
		}

		return numTriangles;
	}





	// ##################################################################
	// Extract the Normals
	// ##################################################################
	void extractNormals(const std::vector<Edge>&     CWEdges,
	                    const std::vector<int>&      CWEdgesBoundaries,
	                    const std::vector<float3>&   coords,
	                    std::vector<float3>&         normals)
	{
		typedef std::vector<Edge>::const_iterator EdgeIterator;

		int numVertices = CWEdgesBoundaries.size() - 1;
		assert(numVertices == int(coords.size()));

		normals.resize(numVertices);

		for(int vi = 0; vi < numVertices; ++vi)
		{
			EdgeIterator eBegin = CWEdges.begin() + CWEdgesBoundaries[vi];
			EdgeIterator eEnd= CWEdges.begin() + CWEdgesBoundaries[vi+1];

			const float3& x   = coords[vi];
			int numNeighbours = std::distance(eBegin, eEnd);
			float3 AccumNormal = make_float3(0,0,0);
			double sumNorm = 0.0;

			for(EdgeIterator e_itr = eBegin; e_itr != eEnd; ++e_itr)
			{
				if( e_itr->hasNext )
				{
					float3 v;
					float3 w;

					v = coords[e_itr->v] - x;
					w = coords[e_itr->vnext] - x;

					float3 localNorm;
					localNorm = cross(w,v);

					// HACK : in opengl CCW faces are front facing so we adopt this convention
					AccumNormal -= localNorm;
					sumNorm     += norm2(localNorm);
				}
			}

			normals[vi] = AccumNormal / norm2(AccumNormal);
		}
	}

	void averageNormals( const std::vector<Triangle>& triangles,
	                     const std::vector<float3>&   coords,
	                     std::vector<float3>&         normals, 
						 std::vector<float3>&		  normals_tri)
	{
		int numVertices = coords.size();
		normals.resize(numVertices);
		normals_tri.resize(triangles.size());
		std::fill( normals.begin(), normals.end(), make_float3(0,0,0));
		std::fill( normals_tri.begin(), normals_tri.end(), make_float3(0,0,0));
		std::vector<float3>::iterator tnrl_itr = normals_tri.begin();
		
		for(std::vector<Triangle>::const_iterator t_itr = triangles.begin(); t_itr != triangles.end(); ++t_itr)
		{
			const float3& v0 = coords[t_itr->v0];
			const float3& v1 = coords[t_itr->v1];
			const float3& v2 = coords[t_itr->v2];

			float3 localNorm( cross(v1-v0, v2-v0) );
			
			*tnrl_itr++ = localNorm;
			//tnrl_itr++;

			normals[t_itr->v0] += localNorm;
			normals[t_itr->v1] += localNorm;
			normals[t_itr->v2] += localNorm;			
		}
		
		/*for( std::vector<float3>::iterator n_itr = normals.begin(); n_itr != normals.end(); ++n_itr){
			double l = norm2(*n_itr);
			if( l != 0.0) *n_itr /= l;
		}*/

		#pragma omp parallel for 
		for (int ni = 0; ni < numVertices; ++ni){
			double l = norm2(normals[ni]);
			if( l != 0.0) normals[ni] /= l;
		}
	}

	void averageNormalsFliped( const std::vector<Triangle>& triangles,
								 const std::vector<float3>&   coords,
								 std::vector<float3>&         normals, 
								 std::vector<float3>&		  normals_tri)
	{
		int numVertices = coords.size();
		normals.resize(numVertices);
		normals_tri.resize(triangles.size());
		std::fill( normals.begin(), normals.end(), make_float3(0,0,0));
		std::fill( normals_tri.begin(), normals_tri.end(), make_float3(0,0,0));
		std::vector<float3>::iterator tnrl_itr = normals_tri.begin();
		
		for(std::vector<Triangle>::const_iterator t_itr = triangles.begin(); t_itr != triangles.end(); ++t_itr)
		{
			const float3& v0 = coords[t_itr->v0];
			const float3& v1 = coords[t_itr->v1];
			const float3& v2 = coords[t_itr->v2];

			float3 localNorm( -cross(v1-v0, v2-v0) );
			
			*tnrl_itr++ = localNorm;
			//tnrl_itr++;

			normals[t_itr->v0] += localNorm;
			normals[t_itr->v1] += localNorm;
			normals[t_itr->v2] += localNorm;			
		}
		
		/*for( std::vector<float3>::iterator n_itr = normals.begin(); n_itr != normals.end(); ++n_itr){
			double l = norm2(*n_itr);
			if( l != 0.0) *n_itr /= l;
		}*/

		#pragma omp parallel for 
		for (int ni = 0; ni < numVertices; ++ni){
			double l = norm2(normals[ni]);
			if( l != 0.0) normals[ni] /= l;
		}
	}


	void averageNormals( const std::vector<Triangle>& triangles,
	                     const std::vector<float3>&   coords,
	                     std::vector<float3>&         normals)
	{
		int numVertices = coords.size();
		normals.resize(numVertices);

		std::fill( normals.begin(), normals.end(), make_float3(0,0,0));

		/*for(std::vector<Triangle>::const_iterator t_itr = triangles.begin(); t_itr != triangles.end(); ++t_itr)
		{
			const float3& v0 = coords[t_itr->v0];
			const float3& v1 = coords[t_itr->v1];
			const float3& v2 = coords[t_itr->v2];

			float3 localNorm( cross(v1-v0, v2-v0) );
			normals[t_itr->v0] += localNorm;
			normals[t_itr->v1] += localNorm;
			normals[t_itr->v2] += localNorm;
		}*/
		
		BOOST_FOREACH(const Triangle& t_itr, triangles){
			const float3& v0 = coords[t_itr.v0];
			const float3& v1 = coords[t_itr.v1];
			const float3& v2 = coords[t_itr.v2];

			float3 localNorm( cross(v1-v0, v2-v0) );
			normals[t_itr.v0] += localNorm;
			normals[t_itr.v1] += localNorm;
			normals[t_itr.v2] += localNorm;
		}
		
		/*for( std::vector<float3>::iterator n_itr = normals.begin(); n_itr != normals.end(); ++n_itr){
			double l = norm2(*n_itr);
			if( l != 0.0) *n_itr /= l;
		}*/
		
		#pragma omp parallel for 
		for (int ni = 0; ni < numVertices; ++ni){
			double l = norm2(normals[ni]);
			if( l != 0.0) normals[ni] /= l;
		}
	}

	void averageNormalsFliped( const	std::vector<Triangle>& triangles,
								const	std::vector<float3>&   coords,
										std::vector<float3>&	normals)
	{	//Paul
		int numVertices = coords.size();
		normals.resize(numVertices);

		std::fill( normals.begin(), normals.end(), make_float3(0,0,0));

		/*for(std::vector<Triangle>::const_iterator t_itr = triangles.begin(); t_itr != triangles.end(); ++t_itr)
		{
			const float3& v0 = coords[t_itr->v0];
			const float3& v1 = coords[t_itr->v1];
			const float3& v2 = coords[t_itr->v2];

			float3 localNorm( cross(v1-v0, v2-v0) );
			normals[t_itr->v0] -= localNorm;
			normals[t_itr->v1] -= localNorm;
			normals[t_itr->v2] -= localNorm;
		}*/

		
		BOOST_FOREACH(const Triangle& t_itr, triangles){
			const float3& v0 = coords[t_itr.v0];
			const float3& v1 = coords[t_itr.v1];
			const float3& v2 = coords[t_itr.v2];

			float3 localNorm( cross(v1-v0, v2-v0) );
			normals[t_itr.v0] -= localNorm;
			normals[t_itr.v1] -= localNorm;
			normals[t_itr.v2] -= localNorm;
		}

		/*for( std::vector<float3>::iterator n_itr = normals.begin(); n_itr != normals.end(); ++n_itr){
			double l = norm2(*n_itr);
			if( l != 0.0) *n_itr /= l;
		}*/

		#pragma omp parallel for 
		for (int ni = 0; ni < numVertices; ++ni){
			double l = norm2(normals[ni]);
			if( l != 0.0) normals[ni] /= l;
		}

	}


	// ##################################################################
	// GetMeanEdgeLength
	// ##################################################################
	double getMeanEdgeLength(const std::vector<Triangle>& triangles,
	                         const std::vector<float3>&   coords)
	{
		if( triangles.size() == 0 ) throw( std::runtime_error(" no triangle in the list" ) );

		double meanEdge = 0.0;
		for( std::vector<Triangle>::const_iterator t_itr = triangles.begin(); t_itr != triangles.end(); ++ t_itr )
		{
			double t1 = norm2( coords[ t_itr->v0] - coords[ t_itr->v1] );
			double t2 = norm2( coords[ t_itr->v1] - coords[ t_itr->v2] );
			double t3 = norm2( coords[ t_itr->v2] - coords[ t_itr->v0] );
			meanEdge += t1+t2+t3;
		}

		return meanEdge / double( 3*triangles.size() );
	}

	double getMeanEdgeLength(const std::vector<Edge>&     CWEdges,
	                         const std::vector<int>&      CWEdgesBoundaries,
	                         const std::vector<float3>&   coords)
	{
		int numVertices = CWEdgesBoundaries.size() - 1;
		if( CWEdges.size() == 0 ) throw( std::runtime_error(" no edges in the list" ) );
		if ( int(coords.size()) != numVertices ) throw( std::runtime_error(" mismatch between the graph and the coords" ) );

		double meanEdge = 0.0;
		for( int vi =0; vi<numVertices; ++vi )
		{
			const float3& ci = coords[vi];
			const std::vector<Edge>::const_iterator eBegin = CWEdges.begin() + CWEdgesBoundaries[vi];
			const std::vector<Edge>::const_iterator eEnd = CWEdges.begin() + CWEdgesBoundaries[vi+1];

			for( std::vector<Edge>::const_iterator e_itr = eBegin; e_itr != eEnd; ++ e_itr ) meanEdge += norm2( ci - coords[ e_itr->v] );
		}

		return meanEdge / double( CWEdges.size() );
	}


	// ##################################################################
	// GetBBox
	// ##################################################################
	void getBBox(const std::vector<float3>& coords, 
	             double& xmin, double& ymin, double& zmin, 
	             double& xmax, double& ymax, double& zmax)
	{
		xmin = ymin = zmin = FLT_MAX;// std::numeric_limits<double>::max();
		xmax = ymax = zmax = FLT_MIN;// std::numeric_limits<double>::min();

		for(std::vector<float3>::const_iterator v_itr = coords.begin(); v_itr != coords.end(); ++ v_itr)
		{
			const float3& v = *v_itr;
			if( xmin > v.x ) xmin = v.x;
			if( xmax < v.x ) xmax = v.x;
			if( ymin > v.y ) ymin = v.y;
			if( ymax < v.y ) ymax = v.y;
			if( zmin > v.z ) zmin = v.z;
			if( zmax < v.z ) zmax = v.z;
		}
	}










	// ##################################################################
	// ajdacency to triangles
	// ##################################################################
	void convertEdgeGraphToTriangles( const std::vector<Edge>& CWEdges,
	                                  const std::vector<int>&  CWEdgesBoundaries,
	                                  std::vector<Triangle>&   triangles)
	{
		typedef std::vector<Edge>::const_iterator EdgeIterator;

		int numVertices = CWEdgesBoundaries.size() - 1;
		int numTriangles = getNumTriangles(CWEdges, CWEdgesBoundaries);
		triangles.clear();
		triangles.reserve(numTriangles);

		for(int vi=0;vi<numVertices;++vi)
		{
			std::vector<Edge>::const_iterator eBegin = CWEdges.begin() + CWEdgesBoundaries[vi];
			std::vector<Edge>::const_iterator eEnd = CWEdges.begin() + CWEdgesBoundaries[vi+1];

			for(std::vector<Edge>::const_iterator e_itr = eBegin; e_itr != eEnd; ++e_itr) {
				if( e_itr->hasNext && e_itr->v < vi && e_itr->vnext < vi) {
					Triangle t;
					t.v0 = vi;
					t.v1 = e_itr->v;
					t.v2 = e_itr->vnext;
					triangles.push_back(t);
				}
			}
		}
		assert( int(triangles.size()) == numTriangles );
	}









	// ##################################################################
	// reindex the mesh
	// It would for sure be easier to do this at the Triangle stage... but graph algorithms operate on adjacency graphs
	// ##################################################################
	void reindexAdjacency( const std::vector<Edge>& CWEdges0,
	                       const std::vector<int>&  CWEdgesBoundaries0,
	                       const std::vector<int>&  oldIndices, // the new vertex at position i had the indice oldIndices[i] in A0
								 std::vector<Edge>& CWEdges,
								 std::vector<int>&  CWEdgesBoundaries)
	{
		typedef std::vector<Edge>::const_iterator EdgeIterator;

		// Size verifications
		int numVertices = oldIndices.size();
		assert( int(CWEdgesBoundaries0.size() - 1) == numVertices );

		// WS
		std::vector<int> newIndices(numVertices);
		for(int vi=0;vi<numVertices;++vi) {
			newIndices[oldIndices[vi]] = vi;
		}


		// Allocation
		CWEdgesBoundaries.resize(numVertices + 1);
		CWEdges.clear();
		CWEdges.reserve( CWEdges0.size() );


		// And lets write the thing
		int offset = 0;
		for(int vi=0;vi<numVertices;++vi)
		{
			CWEdgesBoundaries[vi] = offset;
			const EdgeIterator eBegin_old = CWEdges0.begin() + CWEdgesBoundaries0[oldIndices[vi]];
			const EdgeIterator eEnd_old   = CWEdges0.begin() + CWEdgesBoundaries0[oldIndices[vi]+1];

			for(EdgeIterator e_itr = eBegin_old; e_itr != eEnd_old; ++e_itr) {
				Edge e = *e_itr;
				e.v = newIndices[e.v];
				if( e.hasNext ) e.vnext = newIndices[e.vnext];
				CWEdges.push_back(e);
				offset++;
			}
		}
		CWEdgesBoundaries[numVertices] = offset;
	}










	// ##################################################################
	// triangles to adjacency
	// ##################################################################
	void convertTrianglesToEdgeGraph(const int numVertices,
	                                 const std::vector<Triangle>& triangles,
	                                 std::vector<Edge>&           CWEdges,
	                                 std::vector<int>&            CWEdgesBoundaries)
	{
		int numTriangles = triangles.size();

		// 1 - start by counting how many triangles each vertex will have attached to it
		std::vector<int> vertex_int(numVertices, 0);	
		for(std::vector<Triangle>::const_iterator t_itr = triangles.begin(); t_itr != triangles.end(); ++t_itr)
		{
			vertex_int[t_itr->v0]++;
			vertex_int[t_itr->v1]++;
			vertex_int[t_itr->v2]++;
		}

		// 2 - create the triangle /vertex adjacency vector
		std::vector<Triangle> vertex_Triangles(3*numTriangles); // each vertex has a list to the incident triangles
		std::vector<int>      vertex_Triangles_boundaries(numVertices + 1);
			//2.a -build the boundary vector
		{
			vertex_Triangles_boundaries.resize(numVertices+1);
			int offset = 0;
			for(int vi = 0; vi < numVertices; ++vi) {
				vertex_Triangles_boundaries[vi] = offset;
				offset += vertex_int[vi];
			}
			vertex_Triangles_boundaries[numVertices] = offset;
		}
			// 2.b - write the vertex triangle adj vector
			// this holds 3 times the original triangles, each replicated for its three vertices
		std::fill(vertex_int.begin(), vertex_int.end(), 0); // we are reusing vertex_int to have offsets
		for(std::vector<Triangle>::const_iterator t_itr = triangles.begin(); t_itr != triangles.end(); ++t_itr)
		{
			// we write direct rotations of the triangles
			{
				int offset = vertex_Triangles_boundaries[t_itr->v0] + vertex_int[t_itr->v0];
				vertex_Triangles[offset].v0 = t_itr->v0;
				vertex_Triangles[offset].v1 = t_itr->v1;
				vertex_Triangles[offset].v2 = t_itr->v2;
				vertex_int[t_itr->v0]++;
			}
			{
				int offset = vertex_Triangles_boundaries[t_itr->v1] + vertex_int[t_itr->v1];
				vertex_Triangles[offset].v0 = t_itr->v1;
				vertex_Triangles[offset].v1 = t_itr->v2;
				vertex_Triangles[offset].v2 = t_itr->v0;
				vertex_int[t_itr->v1]++;
			}
			{
				int offset = vertex_Triangles_boundaries[t_itr->v2] + vertex_int[t_itr->v2];
				vertex_Triangles[offset].v0 = t_itr->v2;
				vertex_Triangles[offset].v1 = t_itr->v0;
				vertex_Triangles[offset].v2 = t_itr->v1;
				vertex_int[t_itr->v2]++;
			}
		}






		// 3 - fill the Edge vector
			// 3.a - compute the number of edges
			// we count the number of adjacent vertices, which is actually how many triangles shared one vertex.
		CWEdgesBoundaries.resize(numVertices+1);
		{
			int offset = 0;
			for(int vi=0; vi<numVertices; ++vi)
			{
				CWEdgesBoundaries[vi] = offset;
				std::fill( vertex_int.begin(), vertex_int.end(), 0 );
				std::vector<Triangle>::const_iterator tbegin = vertex_Triangles.begin() + vertex_Triangles_boundaries[vi];
				std::vector<Triangle>::const_iterator tend   = vertex_Triangles.begin() + vertex_Triangles_boundaries[vi+1];
				for(std::vector<Triangle>::const_iterator t_itr = tbegin; t_itr != tend; ++t_itr) {
					vertex_int[t_itr->v1] = 1;
					vertex_int[t_itr->v2] = 1;
				}
				offset += std::accumulate(vertex_int.begin(), vertex_int.end(), 0);
			}
			CWEdgesBoundaries[numVertices] = offset;
		}

		//debug
		std::cout<<"We'll have "<<CWEdgesBoundaries[numVertices]<<" Edges"<<std::endl;
		//end debug


			// 3.b - fill the edge vector
		CWEdges.resize(CWEdgesBoundaries[numVertices]);
		std::vector<std::pair<int, int> > vertex_int_int(numVertices);
		for(int vi=0; vi<numVertices; ++vi)
		{
			std::fill( vertex_int_int.begin(), vertex_int_int.end(), std::make_pair<int,int>(-1,-1) );

			std::vector<Triangle>::const_iterator tbegin = vertex_Triangles.begin() + vertex_Triangles_boundaries[vi];
			std::vector<Triangle>::const_iterator tend   = vertex_Triangles.begin() + vertex_Triangles_boundaries[vi+1];

			int nTris = std::distance(tbegin, tend);
			if( nTris )
			{
				// 4.1 - for each vertex we store the incoming and outgoing edges
				for(std::vector<Triangle>::const_iterator t_itr = tbegin; t_itr != tend; ++t_itr)
				{
					if( vertex_int_int[t_itr->v1].second != -1 ) throw ( std::runtime_error("Not an oriented Manifold.. two edges from the same 1-ring point to the same vertex") );
					vertex_int_int[t_itr->v1].second = t_itr->v2;
					vertex_int_int[t_itr->v2].first = t_itr->v1;
				}

				std::vector<Edge>::iterator e_itr = CWEdges.begin() + CWEdgesBoundaries[vi];
				std::vector<Edge>::iterator eEnd  = CWEdges.begin() + CWEdgesBoundaries[vi+1];
				bool closedFace = false; // if this is true... we should only go one time through the following loop
				int numComponents = 0;
				while( e_itr != eEnd)
				{
					Edge e;
					// 2.a - we find a vertex and try to find the start of the loop (if there is one)
					int vstart0 = 0;
					while( vertex_int_int[vstart0].second == -1 ) vstart0++;
					int vstart = vstart0;
					while( (vertex_int_int[vstart].first != -1) && (vertex_int_int[vstart].first != vstart0) ) vstart = vertex_int_int[vstart].first; //rewind

					closedFace |= (vertex_int_int[vstart].first == vstart0); // iif we looped back... the face has a closed loop
					numComponents++;

					// 2.b - we now unroll the loop
					e.v = vstart;
					while( vertex_int_int[e.v].second != -1 ) // while there is a next triangle
					{
						e.vnext   = vertex_int_int[e.v].second;
						e.hasNext = 1;

						vertex_int_int[e.v].second = -1;  // break loop
						*e_itr++ = e;                     // add to edgelist
						e.v = e.vnext;                    // move forward
					}

					if( !closedFace ) {
						e.vnext = -1;
						e.hasNext = 0;
						*e_itr++ = e;
					}
				}

				if ( closedFace && numComponents > 1 ) throw ( std::runtime_error(" Not an oriented Manifold... a closed loop + other shit") );
			}
		}
	}






	// ##################################################################
	// GetConnectedComponents
	// ##################################################################
	int getConnectedComponents( const std::vector<Edge>& CWEdges,
	                            const std::vector<int>&  CWEdgesBoundaries,
	                            std::vector<int>&       vertex_ccomponent )
	{
		std::cout<<"(W) untested function "<<__FUNCTION__<<std::endl;
		int numVertices = CWEdgesBoundaries.size()-1;

		// resize the output vector
		vertex_ccomponent.resize(numVertices);
		std::fill( vertex_ccomponent.begin(), vertex_ccomponent.end(), -1 );


		int ccId = 0;
		int numUnassigned = numVertices;
		std::queue<int> processQueue;
		while( numUnassigned > 0 )
		{
			int seedIndex = 0;
			while( vertex_ccomponent[seedIndex] != -1 ) seedIndex++;

			vertex_ccomponent[seedIndex] = ccId;
			processQueue.push(seedIndex);
			numUnassigned--;

			while( ! processQueue.empty() )
			{
				int visitedId = processQueue.front();
				processQueue.pop();

				const std::vector<Edge>::const_iterator nBegin = CWEdges.begin() + CWEdgesBoundaries[visitedId];
				const std::vector<Edge>::const_iterator nEnd = CWEdges.begin() + CWEdgesBoundaries[visitedId+1];
				for( std::vector<Edge>::const_iterator n_itr = nBegin; n_itr != nEnd; ++n_itr)
				{
					if( vertex_ccomponent[n_itr->v] == -1 ) {
						vertex_ccomponent[n_itr->v] = ccId;
						processQueue.push(n_itr->v);
						numUnassigned--;
					}
					// TODO (or not) assert on the fact that it is either -1 or ccId
				}
			}

			ccId ++; // move to next connected component id
		}

		return ccId;
	}







}
