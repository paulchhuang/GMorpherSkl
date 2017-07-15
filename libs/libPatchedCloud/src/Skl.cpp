/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#include <Skl.h>
#include <fstream>
#include <cmath>

#define NUM_JOINTS_SUB 16
using std::string;


inline KBNode::NodeType StringToKBNodeType( const string& s ) {
	if     ( s == "KBJOINT_REVOLUTE" )  return KBNode::KBNODE_REVOLUTE;
	else if( s == "KBJOINT_BALL" )      return KBNode::KBNODE_BALL;
	else if( s == "KBJOINT_TRANSLATE" ) return KBNode::KBNODE_TRANSLATE;
	else if( s == "KBJOINT_END" )       return KBNode::KBNODE_END;
	throw std::runtime_error( "not a valid node type" );
}

inline bool isRoot( const KBNode& node ) {
	return node.m_parentName == "NULL";
}

inline int findString( const string& s, const std::vector<string>& vs) {
	std::vector<string>::const_iterator f_itr = std::find(vs.begin(), vs.end(), s );
	if( f_itr == vs.end() ) return -1;
	else return f_itr - vs.begin();
}

inline int findNamedNode( const std::string& name,
                          const std::vector<KBNode>& KBNodes ) {
	for( std::vector<KBNode>::const_iterator n_itr = KBNodes.begin(); 
	                                       n_itr != KBNodes.end(); ++ n_itr ) {
		if (n_itr->m_name == name ) return n_itr - KBNodes.begin();
	}
	return -1;
}


void KBParse( const char*               filename,
              std::vector<KBNode>&      KBNodes)
{
	KBNodes.clear();

	// open file
	std::ifstream fin( filename );
	
	if (!fin.is_open()){ throw std::runtime_error(filename); }
	
	// read header
	string buff;
	int numNodes, numSubs;
	fin >> buff >> numNodes;
	if( buff != "KBSKEL" ) throw std::runtime_error( "KBParse: expected KBSKEL" );

	numSubs = numNodes/NUM_JOINTS_SUB;	
	
	// read joints
	for (int si = 0; si < numSubs; si++){
		for(int ji=0;ji<NUM_JOINTS_SUB;++ji){
			string name, type, parentName;
			fin >> name >> type >> parentName;
			KBNode node;
			node.m_name       = name;
			node.m_parentName = parentName;
			node.m_type       = StringToKBNodeType( type );
			// make sure that the parent is in the map
			if( !isRoot(node) && (findNamedNode(parentName, KBNodes) < 0) ) 
				throw std::runtime_error( "KBParse: the skel is not hierarchically organized");
			// push in the vectors
			KBNodes.push_back(node);
		}
	}
	
	// read parameters 
	for (int si = 0; si < numSubs; si++){
		for(int ji=0;ji < NUM_JOINTS_SUB;++ji) {
			string name;
			fin >> name;
			int nodeId = findNamedNode( name, KBNodes);
			if( nodeId < 0 ) throw std::runtime_error( "KBParse: node with unknown name in params");
			KBNode& node = KBNodes[nodeId + si*NUM_JOINTS_SUB];
			switch( node.m_type ) {
				case KBNode::KBNODE_REVOLUTE:  
					fin >> node.m_x >> node.m_y >> node.m_z 
						>> node.m_ax >>node.m_ay >> node.m_az; break;
				default:
					fin >> node.m_x >> node.m_y >> node.m_z; break;
			}
			/*node.m_x /= 1000;
			node.m_y /= 1000;
			node.m_z /= 1000;*/
		}
	}	
}


void KBGenJoints( const std::vector<KBNode>&      KBNodes,
                  std::vector<std::string>&       KBJointNames,
                  std::vector<KBJoint>&           KBJoints,
                  std::vector<double>&            KBParamsVec,
                  std::vector<double>&            KBStateVec )
{
	KBJointNames.clear();
	KBJoints.clear();
	KBParamsVec.clear();
	KBStateVec.clear();


	int numNodes = KBNodes.size();
	int numSubs	 = numNodes/NUM_JOINTS_SUB;	
	for (int si = 0; si < numSubs; si++){
		for(int ni=0;ni < NUM_JOINTS_SUB;++ni){
			const KBNode& node = KBNodes[ni + si*NUM_JOINTS_SUB];
			switch( node.m_type ) 
			{
				case KBNode::KBNODE_REVOLUTE : {
					int parentId = findString( node.m_parentName, KBJointNames );
					if( !isRoot(node) && ( parentId < 0 ) ) throw std::runtime_error("badly-formed tree");
					KBJoint joint;
					joint.m_parentId = parentId + si*NUM_JOINTS_SUB;
					joint.m_type = KBJoint::KBJOINT_REVOLUTE;
					KBJoints.push_back(joint);
					KBJointNames.push_back(node.m_name);
					KBParamsVec.push_back(node.m_x);  KBParamsVec.push_back(node.m_y); KBParamsVec.push_back(node.m_z);
					//KBParamsVec.push_back(node.m_ax); KBParamsVec.push_back(node.m_ay); KBParamsVec.push_back(node.m_az);
					KBStateVec.push_back(0); // 0 angle
					break;
				}
				case KBNode::KBNODE_BALL : {
					int parentId = findString( node.m_parentName, KBJointNames );
					if( !isRoot(node) && ( parentId < 0 ) ) throw std::runtime_error("badly-formed tree");
					KBJoint joint;
					joint.m_parentId = parentId + si*NUM_JOINTS_SUB;
					joint.m_type = KBJoint::KBJOINT_BALL;
					KBJoints.push_back(joint);
					KBJointNames.push_back(node.m_name);
					KBParamsVec.push_back(node.m_x);  KBParamsVec.push_back(node.m_y); KBParamsVec.push_back(node.m_z);
					KBStateVec.push_back(1); KBStateVec.push_back(0); KBStateVec.push_back(0); KBStateVec.push_back(0); // id quaternion
					break;
				}
				case KBNode::KBNODE_TRANSLATE : {
					int parentId = findString( node.m_parentName, KBJointNames );
					if( !isRoot(node) && ( parentId < 0 ) ) throw std::runtime_error("badly-formed tree");
					KBJoint joint;
					joint.m_parentId = parentId;
					joint.m_type = KBJoint::KBJOINT_TRANSLATE;
					KBJoints.push_back(joint);
					KBJointNames.push_back(node.m_name);
					KBParamsVec.push_back(node.m_x);  KBParamsVec.push_back(node.m_y); KBParamsVec.push_back(node.m_z);
					KBStateVec.push_back(0); KBStateVec.push_back(0);KBStateVec.push_back(0); // 0 translation
					break;
				}
				case KBNode::KBNODE_END :
					// in Kineben before we do nothing here, but in this project we have to consider leaf node.
					int parentId = findString( node.m_parentName, KBJointNames );
					if( !isRoot(node) && ( parentId < 0 ) ) throw std::runtime_error("badly-formed tree");
					KBJoint joint;
					joint.m_parentId = parentId + si*NUM_JOINTS_SUB;
					joint.m_type = KBJoint::KBJOINT_END;
					KBJoints.push_back(joint);
					KBJointNames.push_back(node.m_name);
					KBParamsVec.push_back(node.m_x);  KBParamsVec.push_back(node.m_y); KBParamsVec.push_back(node.m_z);
					break;
			}
		}

	}
		
}


