/* *************************************************
 * Copyright (2015) : Paul Huang
 * *************************************************/
#ifndef JOINTANGLEPRIOR_H_DEFINED
#define JOINTANGLEPRIOR_H_DEFINED

#pragma comment(lib, "libmat.lib")
#pragma comment(lib, "libmx.lib")

#include "mat.h"
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <cstdio>
#include <algorithm>
#include <tuple>
#include <stdexcept>
#include <iostream>
#include "CC3D\float3.h"

namespace GMorpher 
{

	using namespace CC3D;

class JointAnglePrior{
	public:
		JointAnglePrior(const char* const fileName);
		~JointAnglePrior();

		void JointAnglePrior::global2local(const std::vector<double>& S, std::vector<float3>& dSl);
		void JointAnglePrior::local2global(const std::vector<float3>& dSl, std::vector<float3>& dS);
		
		void JointAnglePrior::isValid(const std::vector<double>& S, std::vector<bool>& flag);
		void JointAnglePrior::isValid(const std::vector<double>& S, std::vector<bool>& flag, std::vector<float3>& S2);
		

		inline const std::vector<std::vector<double>>&				arrayData()	const { return m_arrayData; }
		inline const std::vector<std::vector<std::vector<double>>>& cellData()	const { return m_cellData; }
		inline const std::vector<std::vector<std::vector<bool>>>&	BcellData()	const { return m_BcellData; }
		inline const std::map<int, int>&								MappingIjaz2Paul() const { return m_jointMappingIjaz2Paul; }

	private:
		std::vector<std::string>						m_VarName;
		std::map<int, std::string>						m_BPName;
		std::map<std::string, int>						m_arrMap;
		std::map<std::string, int>						m_cellMap;
		std::map<std::string, int>						m_BcellMap;

		std::vector<std::vector<double>> 				m_arrayData;			// 6 variables: edges, a, di, jmp, chlds, prnts
		std::vector<std::vector<std::vector<double>>>	m_cellData;				// 4 variables: sepPlane, bounds, E2, boundries
		std::vector<std::vector<std::vector<bool>>>		m_BcellData;			// 1 variables: angleSprd

		std::vector<int>								m_chldsT;


		
		std::map<int, int>								m_jointMappingIjaz2Paul;		// this mapping maps the index of Ijaz's joints to mine, but it's subject to different dataset.

};




} // end namespace GMorpher













#endif
