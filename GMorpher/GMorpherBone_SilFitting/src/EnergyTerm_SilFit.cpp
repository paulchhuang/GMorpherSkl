/* *************************************************
 * Copyright (2017) : Chun-Hao Paul Huang
 * *************************************************/
//#define DEBUG_SILRENDER
//#define DEBUG_SILCORR

#define NOMINMAX
#include <EnergyTerm_SilFit.h>
#include <EnergyTerm_2DPixelConstraint_ceres.h>
#include <LineRasterizer.h>
/*
 glut is used to create window and initialize context. it replaces OpenGLContextMS and is used only here.
 if ones needs to declare opengl variables, e.g., GLfloat, he needs to include only glew.h
*/
#include <gl\glut.h>		
// misc
#include <OptionParser.h>
#include<fstream>


using CC3D::float3;
using namespace cv;

/**
 * Auxilliary functor...
 * takes an image as input, and computes the gradient where it is requested
 */

void EnergyTerm_SilFit::addResidualBlock(ceres::Problem& problem, double* const RT, double* const JX)
{
	
	int numVertices = m_pv_bounds[m_numPatches];
	// this is a buffer which contains info on whether if a vertex is in a border or not
	std::vector<unsigned char>       hasTarget(numVertices);
	std::vector<std::pair<int, int> > target2D(numVertices);

	
	/*transform to the current shape-pose*/
	// TO-DO: transform dX0_smooth not only dX0
	float3* vGL_itr = (float3*)m_vBuffer;
	for (int pi = 0; pi<m_numPatches; ++pi) {
		const float3*       dX0_itr = m_dX0 + m_pv_bounds[pi];
		const float3* const dX0_end = m_dX0 + m_pv_bounds[pi + 1];

		float rt_32f[6];

		rt_32f[0] = (float)(*(RT + 6 * pi + 0));
		rt_32f[1] = (float)(*(RT + 6 * pi + 1));
		rt_32f[2] = (float)(*(RT + 6 * pi + 2));
		rt_32f[3] = (float)(*(RT + 6 * pi + 3));
		rt_32f[4] = (float)(*(RT + 6 * pi + 4));
		rt_32f[5] = (float)(*(RT + 6 * pi + 5));
	
		while (dX0_itr != dX0_end) {
			// coord
			ceres::AngleAxisRotatePoint(rt_32f, (float*)dX0_itr, (float*)vGL_itr);
			vGL_itr->x += rt_32f[3];
			vGL_itr->y += rt_32f[4];
			vGL_itr->z += rt_32f[5];
			
			vGL_itr++;
			dX0_itr++;
		}
	}
		
	
	// compute the error for all cameras
	for( std::map<int, CamDataPtr>::const_iterator c_itr = m_PhotoData.begin(); c_itr != m_PhotoData.end(); ++c_itr)
	{
		std::fill(hasTarget.begin(), hasTarget.end(), 0x00);
		std::fill(target2D.begin(), target2D.end(), std::pair<int, int>(0, 0));

		//std::cout << "cam: " << c_itr->first << std::endl;
		SilhouetteData& camData = *(c_itr->second);
		int w = camData.w();
		int h = camData.h();
		unsigned char* sil = camData.sil().data;
		const double const* P3x4 = camData.P3x4();
		int umax = w - 1;
		int vmax = h - 1;
	
		OffScreenRendering::SilGLRenderJob job(w, h, camData.GLPMat(), camData.GLMMat(), camData.GLTMat(), m_vBuffer, m_iBuffer, m_iBufferSize);
		m_glRenderer->draw(&job, 0, 0, w, h, m_sil_rendered, m_depth);
			
#ifdef DEBUG_SILRENDER				
		cv::Mat sil_debug = cv::Mat::zeros(camData.h(), camData.w(), CV_8UC1);
		#pragma omp parallel for
		for (int i = 0; i < camData.h()*camData.w(); ++i)	sil_debug.data[i] = (int)m_sil_rendered[i];			//sil.data[i] = (m_sil_border[i]) ? 255: 0;				
				

		cv::namedWindow("Rendered sil", CV_WINDOW_NORMAL);		cv::resizeWindow("Rendered sil", camData.w()*0.5,camData.h()*0.5);		cv::moveWindow("Rendered sil", 200, 100);
		cv::setMouseCallback("Rendered sil", onMouseGRAY, &sil_debug);
		cv::imshow("Rendered sil", sil_debug);

		cv::namedWindow("Observed sil", CV_WINDOW_NORMAL);		cv::resizeWindow("Observed sil", camData.w()*0.5, camData.h()*0.5);		cv::moveWindow("Observed sil", 200, 100);
		cv::setMouseCallback("Observed sil", onMouseGRAY, &camData.sil());
		cv::imshow("Observed sil", camData.sil());

		cv::waitKey(0);							
#endif // DEBUG_SILRENDER
		
		/*int size = w*h;
		int numDiff = 0;
		for( int i=0; i<size; ++i){
			if ( m_sil_rendered[i] ^ sil[i] ) numDiff++;
		}
		E += double(numDiff);		*/

		MVLib::findEdges(w, h, m_sil_rendered, m_sil_border);
		// create the gradientcompute functors
		MVLib::GradientCompute<unsigned char> silrendered_grad(m_sil_rendered, w);
		MVLib::GradientCompute<unsigned char> sil_grad(sil, w);


#ifdef DEBUG_SILCORR
		// ################################################
		// DEBUG
		cv::Mat debug_corr = cv::Mat::zeros(camData.h(), camData.w(), CV_8UC3);
		unsigned char *outputImage = debug_corr.data;
		std::fill(outputImage, outputImage+w*h*3, 0);
		{
			unsigned char* temp = outputImage;
			unsigned char* temp2 = sil;
			unsigned char* temp3 = m_sil_border;
			for(int i=0; i<w*h; ++i) {
				temp[0] = temp[1] = temp[2] = *temp2++ / 2;
				if(*temp3++) {temp[0] = temp[1] = 0; temp[2] = 0xFF; }
				temp+=3;
			}
		}
		// END DEBUG
		// ################################################
#endif //DEBUG_SILCORR

		// --------------------------
		// 1 - find targets
//#pragma omp parallel for
		for (int vi = 0; vi<numVertices; ++vi)
		{
			float3 v = *(const float3*)(m_vBuffer + 3 * vi);
			double xcam = P3x4[0 * 4 + 0] * v.x	+ P3x4[0 * 4 + 1] * v.y	+ P3x4[0 * 4 + 2] * v.z	+ P3x4[0 * 4 + 3];
			double ycam = P3x4[1 * 4 + 0] * v.x	+ P3x4[1 * 4 + 1] * v.y	+ P3x4[1 * 4 + 2] * v.z	+ P3x4[1 * 4 + 3];
			double zcam = P3x4[2 * 4 + 0] * v.x	+ P3x4[2 * 4 + 1] * v.y	+ P3x4[2 * 4 + 2] * v.z	+ P3x4[2 * 4 + 3];

			int u_i = int(xcam / zcam + 0.5);
			int v_i = int(ycam / zcam + 0.5);

			//if querying in the textures is useless, abort
			if ((u_i <= 1) || (u_i >= umax) || (v_i <= 1) || (v_i >= vmax))
				continue;

			int px_coord = u_i + v_i*w;


			if (m_sil_border[px_coord])
			{
				//we know we have at least one pixel border
				double gx, gy;
				silrendered_grad(px_coord, gx, gy); // query the gradient value

				int u2 = u_i + int(gx * 10000);
				int v2 = v_i + int(gy * 10000);

				int maxXsteps = std::min(u_i, w - u_i);
				int maxYsteps = std::min(v_i, h - v_i);
				int maxSteps = std::min(maxXsteps, maxYsteps);
				maxSteps = std::min(maxSteps, 20); // maximum 20 pix distance

				int utarget = 0;
				int vtarget = 0;

				// todo ... try to find a target
				if (fabs(gx) > fabs(gy)){
					LineRasterizerDXsupDY<unsigned char> lr(sil, u_i, v_i, u2, v2, w, 1);
					LineRasterizerDXsupDY<unsigned char> lr2(sil, u_i, v_i, u2, v2, w, 1);
					for (int s = 0; s<maxSteps; ++s)
					{
						double ogx, ogy;
						++lr;
						--lr2;
						if (*lr) {
							int offset = lr.offset();
							sil_grad(offset, ogx, ogy);
							if ((ogx*gx + ogy*gy)>0){ utarget = offset%w; vtarget = offset / w; hasTarget[vi] = 1; break; }
						}
						if (*lr2) {
							int offset = lr2.offset();
							sil_grad(offset, ogx, ogy);
							if ((ogx*gx + ogy*gy)>0){ utarget = offset%w; vtarget = offset / w; hasTarget[vi] = 1; break; }
						}
					}
				}
				else{
					LineRasterizerDYsupDX<unsigned char> lr(sil, u_i, v_i, u2, v2, w, 1);
					LineRasterizerDYsupDX<unsigned char> lr2(sil, u_i, v_i, u2, v2, w, 1);
					for (int s = 0; s<maxSteps; ++s)
					{
						double ogx, ogy;
						++lr;
						--lr2;
						if (*lr) {
							int offset = lr.offset();
							sil_grad(offset, ogx, ogy);
							if ((ogx*gx + ogy*gy)>0){ utarget = offset%w; vtarget = offset / w; hasTarget[vi] = 1; break; }
						}
						if (*lr2) {
							int offset = lr2.offset();
							sil_grad(offset, ogx, ogy);
							if ((ogx*gx + ogy*gy)>0){ utarget = offset%w; vtarget = offset / w; hasTarget[vi] = 1; break; }
						}
					}
				}

				// if we found a target
				if (hasTarget[vi]){
					
					target2D[vi] = std::pair<int, int>(utarget, vtarget);
					ceres::CostFunction* cost_function = GMorpher::Pixel2DConstraint::Create(m_dX0 + vi, target2D[vi], P3x4);
					problem.AddResidualBlock(cost_function,
						new ceres::ScaledLoss(NULL, m_alpha, ceres::TAKE_OWNERSHIP), RT + 6 * m_vp_reidxed[vi]);

#ifdef DEBUG_SILCORR
					// ################################################
					// DEBUG
					int du = abs(utarget - u_i);
					int dv = abs(vtarget - v_i);
					if(du > dv){
						LineRasterizerDXsupDY<unsigned char> lr(outputImage, u_i, v_i, utarget, vtarget, 3*w, 3);
						lr()[0] = 0xFF;
						for(int i =0;i<du;++i) { ++lr; lr()[1] =  0xFF;}
	//						lr()[2] =  0xFF;
					}
					else {
						LineRasterizerDYsupDX<unsigned char> lr(outputImage, u_i, v_i, utarget, vtarget, 3*w, 3);
	//						lr()[0] = 0xFF;
						for(int i =0;i<dv;++i) { ++lr; lr()[1] =  0xFF;}
	//						lr()[2] =  0xFF;
					}
					// END DEBUG
					// ################################################
#endif // DEBUG_SILCORR
				}
			}
		}
#ifdef DEBUG_SILCORR
		cv::namedWindow("corrs", CV_WINDOW_NORMAL);		cv::resizeWindow("corrs", camData.w()*0.5, camData.h()*0.5);		cv::moveWindow("corrs", 200, 100);
		cv::imshow("corrs", debug_corr);
		cv::waitKey(0);
#endif // DEBUG_SILCORR
	}
}



// #############################################################################
// #############################################################################
// EnergyTerm Silhouette Constraint
// #############################################################################
// #############################################################################

EnergyTerm_SilFit::EnergyTerm_SilFit( const std::string& DISPLAYenv,
                                double alphaF,
								const SklRgedPatchedCloud::Patching&     patches,
								const SklRgedPatchedCloud::PatchedCloud& cloud,
								const std::vector<IndexedMesh3D::Triangle>& triangles_r)
{
	m_alpha = alphaF;

	//##########################################################
	//glut thing by Paul, replacing "OpenGLContextMS"
	char *myargv [1];
	int myargc=1;
	myargv [0]=strdup ("Myappname");
	glutInit(&myargc, myargv);
    
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);	//glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);   
    glutInitWindowPosition(5, 5);  
	glutInitWindowSize(10, 10);
	glutCreateWindow("hello world"); 

	GLenum GlewInitResult;
	GlewInitResult = glewInit();
	if (GLEW_OK != GlewInitResult) {
		fprintf(stderr, "ERROR: %s\n", glewGetErrorString(GlewInitResult));
		exit(EXIT_FAILURE);
	}
	//##########################################################
	
	m_glRenderer = new OffScreenRendering::Renderer(64, 64,
	                                        OffScreenRendering::Renderer::GRAY8,
	                                     OffScreenRendering::Renderer::DEPTH32);
	
	// image Buffers
	m_sil_rendered = new unsigned char[2];
	m_sil_border = new unsigned char[2];
	m_depth = new float[2];

	// interpolation buffers 
	int numVertices = patches.numVertices();
	m_numPatches = patches.numPatches();
	m_pv_bounds = patches.pv_boundaries();
	m_vp_reidxed = patches.vp();
	m_dX0 = new float3[numVertices];

	// GL Buffers
	
	std::vector<IndexedMesh3D::Triangle> pv_triangles = triangles_r;
	
	int numTriangles = pv_triangles.size();

	
	m_vBufferSize = 3*numVertices;
	m_iBufferSize = 3 * numTriangles;
	m_vBuffer = new GLfloat[m_vBufferSize];
	m_iBuffer = new GLuint[m_iBufferSize];


	// Copy the triangle data
	for(int ti=0;ti<numTriangles;++ti) {
		m_iBuffer[3 * ti + 0] = pv_triangles[ti].v0;
		m_iBuffer[3 * ti + 1] = pv_triangles[ti].v1;
		m_iBuffer[3 * ti + 2] = pv_triangles[ti].v2;
	}
		
	// copy the vertex data 
	for (int pi = 0; pi<m_numPatches; ++pi) {
		float3 ci = cloud.RT0()[pi].m_t;
		float3* dX0_itr = m_dX0 + m_pv_bounds[pi];
		float3* dX0_end = m_dX0 + m_pv_bounds[pi + 1];
		
		std::vector<float3>::const_iterator X0_itr = cloud.X0().begin()	+ m_pv_bounds[pi];
		std::vector<int>::iterator vp_itr = m_vp_reidxed.begin() + m_pv_bounds[pi];
		//float* v_itr = m_vBuffer + 3 * pv_bounds[pi];

		while( dX0_itr != dX0_end ) {
			/**v_itr = X0_itr->x;		v_itr++;
			*v_itr = X0_itr->y;		v_itr++;
			*v_itr = X0_itr->z;		v_itr++;*/
			//*dX0_itr++ = *X0_itr - ci; 			// Paul: this is bug!
			*dX0_itr = *X0_itr - ci; 
			dX0_itr++;
			X0_itr++;
			*vp_itr++ = pi;
		}
	}		
}


EnergyTerm_SilFit::~EnergyTerm_SilFit()
{
	delete[] m_sil_rendered;
	delete[] m_sil_border;
	delete[] m_depth;

	delete[] m_dX0;
	delete[] m_iBuffer;
	delete[] m_vBuffer;

	delete m_glRenderer;
}


/**
 * Loading camera matrices
 */
void EnergyTerm_SilFit::setCameras( const std::list<int>& ids,
                                 const std::string& camBasename,
                                 const std::string& silBasename,
                                 double znear, double zfar)
{
	//clear existing data
	m_PhotoData.clear();
	delete[] m_sil_rendered;
	delete[] m_sil_border;
	delete[] m_depth;
	delete m_glRenderer;

	int maxW = std::numeric_limits<int>::min();
	int maxH = std::numeric_limits<int>::min();

	// load the cameras
	std::map<int, MVLib::Camera> cameras;
	//std::cout<<"get dimension before"<<std::endl;
	MVLib::loadCameras(camBasename.c_str(),
					   ids,
					   cameras);
	//std::cout<<"get dimension before"<<std::endl;
	
	// create the GL Matrices
	for(std::map<int, MVLib::Camera>::iterator c_itr = cameras.begin(); 
	                                            c_itr != cameras.end(); ++c_itr)
	{
		int id                   = c_itr->first;
		const MVLib::Camera& cam = c_itr->second;
		int w,h;		
		//std::cout<<"get dimension before"<<std::endl;
		cv::Mat img_tmp = cv::imread(MVLib::buildFilename(silBasename.c_str(), id), CV_LOAD_IMAGE_GRAYSCALE);
		w = img_tmp.cols;		h = img_tmp.rows;
		//std::cout<<"get dimension done"<<std::endl;
		CamDataPtr camptr = CamDataPtr(new SilhouetteData(cam, w, h, znear, zfar) );
		m_PhotoData[id] = camptr;

		maxW = std::max( maxW, w);
		maxH = std::max( maxH, h);
	}

	m_sil_rendered = new unsigned char[maxW*maxH];
	m_sil_border = new unsigned char[maxW*maxH];
	m_depth = new float[maxW*maxH];

	m_glRenderer = new OffScreenRendering::Renderer(maxW, maxH,
	                                        OffScreenRendering::Renderer::GRAY8,
	                                     OffScreenRendering::Renderer::DEPTH32);
}

/**
 * Loading image data
 */
void EnergyTerm_SilFit::setImages( const std::string& silBasename )
{
	for( std::map<int, CamDataPtr>::iterator c_itr = m_PhotoData.begin(); 
	                                         c_itr != m_PhotoData.end(); ++c_itr)
	{
		int id                 = c_itr->first;		
		SilhouetteData& camData = *(c_itr->second);
		std::string imageName  = MVLib::buildFilename( silBasename.c_str(), id);				
		camData.sil() = cv::imread(imageName, CV_LOAD_IMAGE_GRAYSCALE);
	}
}
