/***
 * millipede: CubeFace.h
 * Copyright Stuart Golodetz, 2010. All rights reserved.
 ***/

#ifndef H_MILLIPEDE_CUBEFACE
#define H_MILLIPEDE_CUBEFACE

#include <common/util/EnumUtil.h>

namespace mp {

/**
@brief	A CubeFace maintains a map from local face node indices to global node indices.

There are five possible local positions for a node on the cube face, corresponding
to the midpoints of the four edges and the face centre. Not every node will be used
on each cube face.
*/
class CubeFace
{
	//#################### CONSTANTS ####################
public:
	/**
	@brief	An enum representing the possible face node locations.
	*/
	enum NodeDesignator
	{
		TOP_NODE,
		LEFT_NODE,
		MIDDLE_NODE,
		RIGHT_NODE,
		BOTTOM_NODE,
		POTENTIAL_NODE_COUNT,	// = 5, since there are five potential nodes on a face: top, left, middle, right and bottom (numbered here in that order)
	};

	/**
	@brief	An enum representing the locations of the face's vertices.
	*/
	enum VertexDesignator
	{
		TOP_LEFT_VERTEX,
		TOP_RIGHT_VERTEX,
		BOTTOM_LEFT_VERTEX,
		BOTTOM_RIGHT_VERTEX,
	};

private:
	static const int UNUSED = -1;
	static const int USED = -2;

	//#################### NESTED CLASSES ####################
public:
	/**
	@brief	A struct representing cube face edges (i.e. edges between nodes on the cube face).
	*/
	struct Edge
	{
		NodeDesignator u;
		NodeDesignator v;

		Edge(NodeDesignator u_, NodeDesignator v_)
		:	u(u_), v(v_)
		{}
	};

	//#################### PRIVATE VARIABLES ####################
private:
	int m_localToGlobalNodeMap[POTENTIAL_NODE_COUNT];

	//#################### CONSTRUCTORS ####################
public:
	/**
	@brief	Constructs a cube face in which all the local nodes are initially unused.
	*/
	CubeFace();

	//#################### PUBLIC METHODS ####################
public:
	/**
	@brief	Returns the global node index of the specified local node.

	@param[in]	n	The node designator of the local node
	@return	The global node index of the local node, if any, or an unspecified value < 0 otherwise
	*/
	int global_node_index(NodeDesignator n) const;

	/**
	@brief	Returns whether or not the specified local node is in use.

	@param[in]	n	The node designator of the local node
	@return	true, if it is in use, or false otherwise
	*/
	bool is_used(NodeDesignator n) const;

	/**
	@brief	Sets the global node index of the specified local node.

	@param[in]	n		The node designator of the local node
	@param[in]	index	The global node index to assign to it
	@post
		-	global_node_index(n) == index
	*/
	void set_global_node_index(NodeDesignator n, int index);

	/**
	@brief	Marks the specified local node as used, without setting its global node index yet.

	@param[in]	n	The designator of the local node
	@post
		-	is_used(n) == true
	*/
	void set_used(NodeDesignator n);
};

//#################### ENUMERATION HELPERS ####################
template <> CubeFace::NodeDesignator enum_begin();
template <> CubeFace::NodeDesignator enum_end();
CubeFace::NodeDesignator& operator++(CubeFace::NodeDesignator& e);

}

#endif
