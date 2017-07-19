/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef ENERGYTERM_SILFIT_H_DEFINED
#define ENERGYTERM_SILFIT_H_DEFINED
// MVLib
#include <MVLib.h>
//MetaImg
#include <MetaImg.h>
// libGMorpher
#include <IEnergyTerm_ceres.h>
// libPatchedMesh
#include <PatchedCloud.h>
// libRendererOffScreen
#include <Renderer.h>
// std
#include <list>
#include <map>
#include <boost/shared_ptr.hpp>
// OpenCV
#include <cv.h>
#include <opencv2/highgui/highgui.hpp>

class EnergyTerm_SilFit : public GMorpher::IEnergyTermCeres
{
	public :

		EnergyTerm_SilFit(const std::string& DISPLAYenv,
	               double alphaF,
				   const SklRgedPatchedCloud::Patching&     patches,
				   const SklRgedPatchedCloud::PatchedCloud& cloud,
				   const std::vector<IndexedMesh3D::Triangle>& pv_triangles );
		~EnergyTerm_SilFit();

	// public stuff for the application
	void setCameras( const std::list<int>& ids,
	                 const std::string& camBasename,
	                 const std::string& silBasename,
	                 double znear, double zfar);

	void setImages( const std::string& silBasename );


	// Interfaces for the solver
	virtual void addResidualBlock(ceres::Problem& problem, double* const RT, double* const JX);


	protected :

	typedef boost::shared_ptr<SilhouetteData> CamDataPtr;


	// data
	double							m_alpha;
	std::map<int, CamDataPtr>       m_PhotoData;
	// Rendering tools	
	OffScreenRendering::Renderer*	m_glRenderer;
	// image buffers
	unsigned char*					m_sil_rendered;
	unsigned char*					m_sil_border;
	float*							m_depth;
	// interpolation buffers 
	int								m_numPatches;
	std::vector<int>				m_pv_bounds;
	std::vector<int>				m_vp_reidxed;
	CC3D::float3*					m_dX0;
	// GL buffers 
	GLfloat*						m_vBuffer;
	int								m_vBufferSize;
	GLuint*							m_iBuffer;
	int								m_iBufferSize;	
};












#endif
