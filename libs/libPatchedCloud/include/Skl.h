/* *************************************************
 * Copyright (2017) : Cedric Cagniart
 * *************************************************/
#ifndef SKL_H_DEFINED
#define SKL_H_DEFINED

#include <vector>
#include <string>

typedef struct JointNode{
	enum {
		KBJOINT_REVOLUTE,
		KBJOINT_BALL,
		KBJOINT_TRANSLATE,
		KBJOINT_END
	};
	int m_type;
	int m_parentId; /* is -1 if none */
} KBJoint;


typedef struct {
	enum NodeType { KBNODE_REVOLUTE, KBNODE_BALL, KBNODE_TRANSLATE, KBNODE_END };
	std::string m_name;
	std::string m_parentName;
	NodeType    m_type;
	float       m_x, m_y, m_z;
	float       m_ax, m_ay, m_az;
}KBNode ;



void KBParse( const char*               filename,
              std::vector<KBNode>&      KBNodes);

void KBGenJoints( const std::vector<KBNode>&      KBNodes,
                  std::vector<std::string>&       KBJointNames,
                  std::vector<KBJoint>&           KBJoints,
                  std::vector<double>&            KBParamsVec,
                  std::vector<double>&            KBStateVec );

//void KBGenBones( const std::vector<KBNode>&           KBNodes,
//                 const std::vector<std::string>&      KBJointNames,
//                 std::vector<bone>&                   bones );
//
//void KBGenAxisBones( const std::vector<KBNode>&           KBNodes, 
//                     const std::vector<std::string>&      KBJointNames,
//                     std::vector<bone>&                   bones );

class KBBVHWriter 
{
	public :
	KBBVHWriter( const char* filename, 
	             const std::vector<KBNode>&        KBNodes,
	             const std::vector<std::string>&   KBJointNames,
				 const std::vector<KBJoint>&     KBJoints);

	~KBBVHWriter();

	void writeFrame( const std::vector<double>& KBStateVec );

	struct Priv;
	protected :
	Priv* priv;
};

#endif
