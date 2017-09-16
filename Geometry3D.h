/***************************************************************************
*                                                                          *
*   Copyright (c) 2017                                                     *
*   FastFieldSolvers S.R.L.  http://www.fastfieldsolvers.com               *
*                                                                          *
*   This program is free software; you can redistribute it and/or modify   *
*   it under the terms of the GNU Lesser General Public License (LGPL)     *
*   as published by the Free Software Foundation; either version 2 of      *
*   the License, or (at your option) any later version.                    *
*   for detail see the LICENCE text file.                                  *
*                                                                          *
*   This program is distributed in the hope that it will be useful,        *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
*   GNU Library General Public License for more details.                   *
*                                                                          *
*   You should have received a copy of the GNU Library General Public      *
*   License along with this program; if not, write to the Free Software    *
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307   *
*   USA                                                                    *
*                                                                          *
***************************************************************************/


// Geometry3D.h : definition of geometry-related classes
// Enrico Di Lorenzo 2003/03/12

#ifndef C3D_GEOMETRY_DEFS
#define C3D_GEOMETRY_DEFS

#include <math.h>

//#include "Global.h"
#include "List.h"
#include "GeoGlobal.h"
#include "Vector3D.h"
#include "Vector2D.h"

/////////////////////////////////////////////////////
// 3D objects classes
// these are the classes which describes 3D objects
//


// C3D object types
//
#define C3D_UNINIT_ID		0x00000000
#define C3D_LISTHEAD_ID		0x00000001
#define C3D_POINTER_ID		0x00000002
#define C3D_VERTEX_ID		0x00010000
#define C3D_VERTEXUSE_ID	0x00020000
#define C3D_EDGE_ID			0x00040000
#define C3D_EDGEUSE_ID		0x00080000
#define C3D_LOOP_ID			0x00100000
#define C3D_FACEGEO_ID		0x00110000
#define C3D_FACE_ID			0x00120000
#define C3D_SHELL_ID		0x00140000
// edgeuse direction w.r.t. edge
#define C3D_D_SAME		0
#define C3D_D_OPPOSITE	1

// transform operation codes
#define C3D_TRANS_TRANSLATE			0
#define C3D_TRANS_SCALE				1
#define C3D_TRANS_ROTATEX			2
#define C3D_TRANS_ROTATEY			3
#define C3D_TRANS_ROTATEZ			4

// Binary Space Partition tree definitions
//
#define C3D_BSP_NONE		0
#define C3D_BSP_LEFT		1
#define C3D_BSP_RIGHT		2
#define C3D_BSP_ONPLANE		3

// anticipate class definition
class C3DOperation;

// declared classes
class C3DPlane;
class C3DBBox;
class C3DIntersection;
class C3DEntity;
class C3DVertex;
class C3DVertexUse;
class C3DEdge;
class C3DEdgeUse;
class C3DLoop;
class C3DFace;
class C3DShell;
class C3DRegion;
class C3DModel;


// A plane in 3D is represented with the equation
// ax + by + cz + d = 0
// where a,b,c are the components of the plane normal vector
// and -d is the distance of the plane from the origin along normal vector
// ( a(x-x0) + b(y-y0) + c(z-z0) = 0 is the equation, where (x0,y0,z0)
//   is a point on the plane. Then, -d = (a*x0+b*y0+c*z0) = dot_prod(normal,(x0,y0,z0))
//   is the distance of the plane from the origin along normal vector)
// Here, 'm_vecNormal' is [a b c] and 'm_dD' is -d
class C3DPlane
{
public:
	C3DVector	m_vecNormal;
	double		m_dD;

	// default constructor
	C3DPlane()
	{
		m_vecNormal.x = 0.0;
		m_vecNormal.y = 0.0;
		m_vecNormal.z = 1.0;
		m_dD = 0.0;
	}

	// constructor from plane normal and distance from origin
	C3DPlane(C3DVector normal, double dist) : m_vecNormal(normal), m_dD(dist)
	{
	}

	void InvertNormal()
	{
		m_vecNormal.Invert();
		m_dD = -m_dD;
	}

	// Find a point on the plane, given only two coordinates
	// To find the point x0, solve
	// Dot_prod(m_vecNormal, x0) = m_dD
	C3DVector Point3Dfrom2D(C2DVector &vec2d, char c1, char c2, char missingc)
	{
		C3DVector vec3d;

		vec3d[c1] = vec2d.x;
		vec3d[c2] = vec2d.y;

		// avoid divide by zero!
		if( m_vecNormal[missingc] > C3D_TOL || m_vecNormal[missingc] < -C3D_TOL)
			vec3d[missingc] = (m_dD - m_vecNormal[c1]*vec2d.x - m_vecNormal[c2]*vec2d.y) / m_vecNormal[missingc];
		else
			vec3d[missingc] = m_dD;

		return vec3d;
	}
};

class C3DBBox
{
public:
	bool			empty;
	C3DVector		min_point;
	C3DVector		max_point;

	C3DBBox() : empty(true) {}

	void Clear()
	{
		min_point = C3DVector(0.0, 0.0, 0.0);
		max_point = C3DVector(0.0, 0.0, 0.0);
		empty = true;
	}

	unsigned char MaxSide()
	{
		double sidex, sidey, sidez;

		sidex = max_point.x - min_point.x;
		sidey = max_point.y - min_point.y;
		sidez = max_point.z - min_point.z;

		if( (sidex >= sidey) && (sidex >= sidez) ) {
			return 0;
		}
		else if( (sidey >= sidex) && (sidey >= sidez) ) {
			return 1;
		}
		else {
			return 2;
		}
	}

	double MaxSideLen()
	{
		double sidex, sidey, sidez;

		sidex = max_point.x - min_point.x;
		sidey = max_point.y - min_point.y;
		sidez = max_point.z - min_point.z;

		if( (sidex >= sidey) && (sidex >= sidez) ) {
			return sidex;
		}
		else if( (sidey >= sidex) && (sidey >= sidez) ) {
			return sidey;
		}
		else {
			return sidez;
		}
	}

	inline void operator=(const C3DBBox &bbox)
	{
		if (bbox.empty == false) {
			min_point = bbox.min_point;
			max_point = bbox.max_point;
			empty = false;
		}
		else
			empty = true;
	}

	inline void operator=(C3DVector &point)
	{
		min_point = point;
		max_point = point;
		empty = false;
	}

	inline void operator+=(C3DBBox &bbox)
	{
		if (bbox.empty == false) {
			if (empty == false) {
				min_point.x = bbox.min_point.x < min_point.x ? bbox.min_point.x : min_point.x;
				min_point.y = bbox.min_point.y < min_point.y ? bbox.min_point.y : min_point.y;
				min_point.z = bbox.min_point.z < min_point.z ? bbox.min_point.z : min_point.z;
				max_point.x = bbox.max_point.x > max_point.x ? bbox.max_point.x : max_point.x;
				max_point.y = bbox.max_point.y > max_point.y ? bbox.max_point.y : max_point.y;
				max_point.z = bbox.max_point.z > max_point.z ? bbox.max_point.z : max_point.z;
			}
			else {
				min_point = bbox.min_point;
				max_point = bbox.max_point;
				empty = false;
			}
		}
	}

	inline void operator+=(C3DVector &point)
	{
		if (empty == false) {
			min_point.x = point.x < min_point.x ? point.x : min_point.x;
			min_point.y = point.y < min_point.y ? point.y : min_point.y;
			min_point.z = point.z < min_point.z ? point.z : min_point.z;
			max_point.x = point.x > max_point.x ? point.x : max_point.x;
			max_point.y = point.y > max_point.y ? point.y : max_point.y;
			max_point.z = point.z > max_point.z ? point.z : max_point.z;
		}
		else {
			min_point = point;
			max_point = point;
			empty = false;
		}
	}

	inline void operator+=(C3DVector_float &point)
	{
		if (empty == false) {
			min_point.x = point.x < min_point.x ? point.x : min_point.x;
			min_point.y = point.y < min_point.y ? point.y : min_point.y;
			min_point.z = point.z < min_point.z ? point.z : min_point.z;
			max_point.x = point.x > max_point.x ? point.x : max_point.x;
			max_point.y = point.y > max_point.y ? point.y : max_point.y;
			max_point.z = point.z > max_point.z ? point.z : max_point.z;
		}
		else {
			min_point = point;
			max_point = point;
			empty = false;
		}
	}

	inline void operator+=(C2DVector &point)
	{
		if (empty == false) {
			min_point.x = point.x < min_point.x ? point.x : min_point.x;
			min_point.y = point.y < min_point.y ? point.y : min_point.y;
			max_point.x = point.x > max_point.x ? point.x : max_point.x;
			max_point.y = point.y > max_point.y ? point.y : max_point.y;
		}
		else {
			min_point.pos(point.x, point.y, 0.0);
			max_point.pos(point.x, point.y, 0.0);
			empty = false;
		}
	}

	inline void operator+=(C2DVector_float &point)
	{
		if (empty == false) {
			min_point.x = point.x < min_point.x ? point.x : min_point.x;
			min_point.y = point.y < min_point.y ? point.y : min_point.y;
			max_point.x = point.x > max_point.x ? point.x : max_point.x;
			max_point.y = point.y > max_point.y ? point.y : max_point.y;
		}
		else {
			min_point.pos(point.x, point.y, 0.0);
			max_point.pos(point.x, point.y, 0.0);
			empty = false;
		}
	}

	inline C3DBBox operator+(C3DBBox &bbox)
	{
		C3DBBox newBbox;

		if (bbox.empty == false) {
			if( empty == false) {
				newBbox.min_point.x = bbox.min_point.x < min_point.x ? bbox.min_point.x : min_point.x;
				newBbox.min_point.y = bbox.min_point.y < min_point.y ? bbox.min_point.y : min_point.y;
				newBbox.min_point.z = bbox.min_point.z < min_point.z ? bbox.min_point.z : min_point.z;
				newBbox.max_point.x = bbox.max_point.x > max_point.x ? bbox.max_point.x : max_point.x;
				newBbox.max_point.y = bbox.max_point.y > max_point.y ? bbox.max_point.y : max_point.y;
				newBbox.max_point.z = bbox.max_point.z > max_point.z ? bbox.max_point.z : max_point.z;
				newBbox.empty = false;
			}
			else {
				newBbox.min_point = bbox.min_point;
				newBbox.max_point = bbox.max_point;
				newBbox.empty = false;
			}
		}
		else {
			if (empty == false) {
				newBbox.min_point = min_point;
				newBbox.max_point = max_point;
				newBbox.empty = false;
			}
			else {
				newBbox.empty = true;
			}
		}

		return newBbox;
	}

	inline C3DBBox operator+(C3DVector &point)
	{
		C3DBBox newBbox;

		if (empty == false) {
			newBbox.min_point.x = point.x < min_point.x ? point.x : min_point.x;
			newBbox.min_point.y = point.y < min_point.y ? point.y : min_point.y;
			newBbox.min_point.z = point.z < min_point.z ? point.z : min_point.z;
			newBbox.max_point.x = point.x > max_point.x ? point.x : max_point.x;
			newBbox.max_point.y = point.y > max_point.y ? point.y : max_point.y;
			newBbox.max_point.z = point.z > max_point.z ? point.z : max_point.z;
			newBbox.empty = false;
		}
		else {
			newBbox.min_point = point;
			newBbox.max_point = point;
			newBbox.empty = false;
		}

		return newBbox;
	}

};

/*
class C3DIntersection
{
public:
	// Intersection can be either a single point,
	// defined by two coincident endpoints,
	// or a segment, defined by two different endpoints.
	C3DVector p1, p2;
};

class C3DEntity
{
public:
	long			id;
	long			index;
	char			status;

	C3DEntity() : id(C3D_UNINIT_ID), index(0), status(0) {}
	C3DEntity(long newid) : id(newid), index(0), status(0) {}
};

class C3DVertex : public C3DEntity
{
public:
	// vertex spatial position
	C3DVector		pos;
	// radial vertexes list
	C3DList<C3DVertexUse*> vertexUseList;

	// constructors & destructors
	C3DVertex(C3DVector newpos) : C3DEntity(C3D_VERTEX_ID), pos(newpos) {}
	~C3DVertex() {}

	void SetPos(C3DVector newpos)
	{
		pos = newpos;
	}
};

class C3DVertexUse : public C3DEntity
{
public:
	// pointer to vertex used
	C3DVertex		*pVertex;
	// vertexuse father: can be shell, loop or edge
	union {
		C3DEntity		*pFather;
		C3DShell		*pFatherShell;
		C3DLoop			*pFatherLoop;
		C3DEdge			*pFatherEdge;
	};

	// constructors & destructors
	C3DVertexUse(C3DVector *newPos) : C3DEntity(C3D_VERTEXUSE_ID), pFatherEdge(NULL)
	{
		C3DVertex newVert(*newPos);
		CBSPTree::iterator it;

		// if there's no vertex alredy in given position
		if( vertexTree.find(it, &newVert) == false ) {
			// create new vertex
			pVertex = new C3DVertex(*newPos);
			// and insert it in tree
			vertexTree.insert(it, pVertex);
		}
		else
			// otherwise reference to existing one
			pVertex = *it;

		// add new vertexuse to vertex vertexuse list
		pVertex->vertexUseList.push_front(this);

	}

	~C3DVertexUse()
	{
		C3DList<C3DVertexUse*>::iterator it;

		// remove vertexuse from vertex vertexUseList

		it = pVertex->vertexUseList.find(this);
		// the vertexuse must be in vertexuse list
		ASSERT(*it != NULL);
		pVertex->vertexUseList.pop(it);

		// remove vertex from BSP tree and delete it,
		// if this was its last vertexuse
		if( pVertex->vertexUseList.is_empty() == true ) {
			vertexTree.erase(pVertex);
			delete pVertex;
		}
	}

	void SetFather(C3DEdge *pFather)
	{
		pFatherEdge = pFather;
	}
	void SetFather(C3DLoop *pFather)
	{
		pFatherLoop = pFather;
	}
	void SetFather(C3DShell *pFather)
	{
		pFatherShell = pFather;
	}


protected:
	class CBSPTree;

	// The BSP tree must be in common with all C3DVertexUse,
	// therefore is declared static
	static CBSPTree vertexTree;

	// Cartesian Binary Space Partition tree,
	// for maintaining unique vertex positions in space
	class CBSPTree
	{
	public:
		// tree node
		class Node
		{
		public:
			char divDir;
			union {
				double divPlane;
				C3DVertex *pObj;
			};
			bool child;
			Node *up;
			Node *left;
			Node *right;

			Node()
			{
			}
			Node(Node *upnode, C3DVertex *ppoint) :
				up(upnode), pObj(ppoint), left(NULL), right(NULL), child(true)
			{
			}
		};

		// CBSPTree iterator
		class iterator;
		friend class iterator;

		class iterator
		{
			friend class CBSPTree;
			Node *currentNode;
		public:
			iterator() : currentNode(NULL)
			{
			}
			iterator(CBSPTree& starttree) : currentNode(starttree.root)
			{
			}
			iterator(iterator &it) : currentNode(it.currentNode)
			{
			}
			iterator(Node *node) : currentNode(node)
			{
			}
			C3DVertex *operator*() {
				if( currentNode == NULL )
					return NULL;

				if(currentNode->child == true)
					return currentNode->pObj;
				else
					return NULL;
			}
		};

	protected:
		// tree variables
		Node *root;
		double m_dTol;

	public:
		CBSPTree() : root(NULL), m_dTol(C3D_TOL) {}
		CBSPTree(C3DVertex *pObj) : m_dTol(C3D_TOL)
		{
			init(pObj);
		}
		void init(C3DVertex *pObj)
		{
			root = new Node;
			root->divDir = C3D_DIR_NONE;
			root->child = true;
			root->up = NULL;
			root->left = NULL;
			root->right = NULL;
			root->pObj = pObj;
		}
		void insert(C3DVertex *ppoint)
		{
			// if first vertex in tree, must init
			if(root == NULL)
				init(ppoint);
			else
				insert(root, ppoint);
		}
		void insert(iterator it, C3DVertex *ppoint)
		{
			// if first vertex in tree, must init
			if(root == NULL)
				init(ppoint);
			// otherwise, *it must be != NULL
			else {
				ASSERT(*it);
				insert(it.currentNode, ppoint);
			}
		}
		bool find(iterator &it, C3DVertex *ppoint)
		{
			Node *node;
			bool res;

			res = find(root, ppoint, &node);

			it.currentNode = node;

			return res;
		}
		bool erase(iterator &it)
		{
			Node *father, *thischild, *otherchild;

			// check iterator; this validates not only
			// the object but also the pointer and the
			// childness status
			if(*it == NULL)
				return false;

			// remove point from tree

			// the old father will become the new node
			thischild = it.currentNode;
			father = thischild->up;

			// if father is NULL, this should be the root
			if(father == NULL) {
				ASSERT(thischild == root);
				// delete the only node in the tree
				delete root;
				// and set root to NULL, so everybody knows
				// that the tree is empty
				root = NULL;

				return true;
			}

			// find out the other node between the two children
			if(father->left == thischild)
				otherchild = father->right;
			else
				otherchild = father->left;

			ASSERT(thischild);
			ASSERT(otherchild);

			// Copy from other child to father, according to otherchild's
			// childness status.
			// Copy everything but 'up', which must stay the same
			if( otherchild->child == true ) {
				father->pObj = otherchild->pObj;
				father->child = true;
				father->left = father->right = NULL;
			}
			else {
				father->divDir = otherchild->divDir;
				father->divPlane = otherchild->divPlane;
				father->child = false;
				father->left = otherchild->left;
				father->right = otherchild->right;
				father->left->up = father;
				father->right->up = father;
			}

			// delete children
			delete thischild;
			delete otherchild;

			return true;

		}
		void erase(C3DVertex *ppoint)
		{
			iterator it;

			if(find(it, ppoint))
				erase(it);
		}

	protected:
		char planeSide(C3DVertex *ppoint, Node *node);
		bool insert(Node *node, C3DVertex *ppoint);
		bool find(Node *node, C3DVertex *ppoint, Node **lastnode);


	};

};

// TBC warning: remember that when father is set, bbox must propagate!


class C3DEdge : public C3DEntity
{
public:
	// pointers to edge vertex uses
	C3DVertexUse	*pVertexUse1;
	C3DVertexUse	*pVertexUse2;
	// radial edges list
	C3DList<C3DEdgeUse*> edgeUseList;

	// constructors & destructors
	C3DEdge(C3DVertexUse *vu1, C3DVertexUse *vu2) : C3DEntity(C3D_EDGE_ID), pVertexUse1(vu1),
													pVertexUse2(vu2)
	{
		// set vertexuses father
		pVertexUse1->pFatherEdge = pVertexUse2->pFatherEdge = this;
	}

	~C3DEdge()
	{
		// must exist!
		ASSERT(pVertexUse1);
		ASSERT(pVertexUse2);

		delete pVertexUse1;
		delete pVertexUse2;
	}
};

class C3DEdgeUse : public C3DEntity
{
public:
	// pointer to edge used
	C3DEdge		*pEdge;
	// direction w.r.t. pEdge
	char		direction;
	// edgeuse father: can be shell, loop or edge
	union {
		C3DEntity		*pFather;
		C3DLoop			*pFatherLoop;
		C3DShell		*pFatherShell;
	};

	// constructors & destructors

	// using this version of the constructor, vu1 and vu2 must be
	// two already-existing vertexuses; they should not have been
	// created ad-hoc, since it's the edge which references the
	// vertexuses, not the edgeuse!
	C3DEdgeUse(C3DVertexUse *vu1, C3DVertexUse *vu2, C3DLoop *father) : C3DEntity(C3D_EDGEUSE_ID),
															pFatherLoop(father), pEdge(NULL)
	{
		CreateEdgeUse(vu1, vu2);
	}

	C3DEdgeUse(C3DVertexUse *vu1, C3DVertexUse *vu2, C3DShell *father) : C3DEntity(C3D_EDGEUSE_ID),
															pFatherShell(father), pEdge(NULL)
	{
		CreateEdgeUse(vu1, vu2);
	}

	C3DEdgeUse(C3DVector *v1, C3DVector *v2, C3DLoop *father) : C3DEntity(C3D_EDGEUSE_ID),
															pFatherLoop(father), pEdge(NULL)
	{
		C3DVertexUse *vu1, *vu2;

		// Must create vertexuses to know if there already is
		// an edge connecting the two given points v1 and v2.
		// If not, can create a new edge referencing the two
		// edgeuses; but if exists, it does already
		// reference two other vertexuses, so must delete
		// these ones (since the edgeuse references the edge,
		// NOT the vertexuses!)

		vu1 = new C3DVertexUse(v1);
		vu2 = new C3DVertexUse(v2);

		if(CreateEdgeUse(vu1, vu2) == true) {
			delete vu1;
			delete vu2;
		}
	}

	C3DEdgeUse(C3DVector *v1, C3DVector *v2, C3DShell *father) : C3DEntity(C3D_EDGEUSE_ID),
															pFatherShell(father), pEdge(NULL)
	{
		C3DVertexUse *vu1, *vu2;

		// Must create vertexuses to know if there already is
		// an edge connecting the two given points v1 and v2.
		// If not, can create a new edge referencing the two
		// edgeuses; but if exists, it does already
		// reference two other vertexuses, so must delete
		// these ones (since the edgeuse references the edge,
		// NOT the vertexuses!)

		vu1 = new C3DVertexUse(v1);
		vu2 = new C3DVertexUse(v2);

		if(CreateEdgeUse(vu1, vu2) == true) {
			delete vu1;
			delete vu2;
		}
	}

	~C3DEdgeUse()
	{
		C3DList<C3DEdgeUse*>::iterator it;

		// remove edgeuse from edge's edgeUseList

		it = pEdge->edgeUseList.find(this);
		// the edgeuse must be in vertexuse list
		ASSERT(*it != NULL);
		pEdge->edgeUseList.pop(it);

		// delete edge, if this was its last edgeuse
		if( pEdge->edgeUseList.is_empty() == true ) {
			delete pEdge;
		}
	}

	C3DVertexUse *GetFirstVertexUse()
	{
		if(direction == C3D_D_SAME)
			return pEdge->pVertexUse1;
		else
			return pEdge->pVertexUse2;
	}

	C3DVertexUse *GetSecondVertexUse()
	{
		if(direction == C3D_D_SAME)
			return pEdge->pVertexUse2;
		else
			return pEdge->pVertexUse1;
	}

	void InvertDirection()
	{
		if( direction == C3D_D_SAME)
			direction = C3D_D_OPPOSITE;
		else
			direction = C3D_D_SAME;
	}

private:
	// if edge between given vertexuses already exists, return true
	bool CreateEdgeUse(C3DVertexUse *vu1, C3DVertexUse *vu2)
	{
		C3DVertex *v1, *v2;
		C3DEdge *edge;
		C3DList<C3DVertexUse*>::iterator it;
		bool reverseDir, ret;

		// set return value to 'edge not existing'
		ret = false;
		// do not reverse edgeuse direction w.r.t edge
		reverseDir = false;

		// get vertexuse's vertexes
		v1 = vu1->pVertex;
		v2 = vu2->pVertex;

		// v1 is the vertex with the smaller vertexuse list
		if( v2->vertexUseList.size() < v1->vertexUseList.size()) {
			v1 = vu2->pVertex;
			v2 = vu1->pVertex;
			// reverse direction, because going to scan edges going
			// from vu2 to vu1
			reverseDir = !reverseDir;
		}

		// scan every vertexues in vertex list
		for(it = v1->vertexUseList.begin(); it != v1->vertexUseList.end(); it++) {
			// copy reference to father in local variable
			edge = (*it)->pFatherEdge;
			// if vertexuse father exists
			if( edge != NULL ) {
				// if father is really an edge
				if( edge->id == C3D_EDGE_ID ) {
					// if an edge which connects the two given vertexuse vertexes
					// already exists, this will be another edgeuse of its
					if( edge->pVertexUse1->pVertex == v2 ) {
						// reference the existing edge
						pEdge = edge;
						// edge exists
						ret = true;
						// reverse direction, because edge goes from v2 to v1
						reverseDir = !reverseDir;
						break;
					}
					if( edge->pVertexUse2->pVertex == v2 ) {
						// reference the existing edge
						pEdge = edge;
						// edge exists
						ret = true;
						break;
					}
				}
			}
		}

		// no edge already exists, must create one
		if(pEdge == NULL) {
			pEdge = new C3DEdge(vu1, vu2);
			// in this case direction is of course the same
			direction = C3D_D_SAME;
		}
		else
			// reflect direction w.r.t. edge
			reverseDir == true ? direction = C3D_D_OPPOSITE : direction = C3D_D_SAME;

		// add new edgeuse to edge's edgeuse list
		pEdge->edgeUseList.push_front(this);

		return ret;
	}

};

class C3DLoop : public C3DEntity
{
public:
	// Both edgeuses and vertex can be child. If list of edgeuses
	// in loop is empty, then vertexuse is the child and vice-versa.
	C3DList<C3DEdgeUse*>	edgeUseList;
	C3DVertexUse	    *pVertexUse;
	// loop bounding box
	C3DBBox			bbox;
	// loop father: can be face or shell
	union {
		C3DEntity		*pFather;
		C3DFace			*pFatherFace;
		C3DShell		*pFatherShell;
	};

	// constructors & destructors
	C3DLoop() : C3DEntity(C3D_LOOP_ID), pVertexUse(NULL), pFatherFace(NULL) {}
	C3DLoop(C3DFace *father) : C3DEntity(C3D_LOOP_ID), pVertexUse(NULL), pFatherFace(father) {}
	C3DLoop(C3DShell *father) : C3DEntity(C3D_LOOP_ID), pVertexUse(NULL), pFatherShell(father) {}
	C3DLoop(C3DVertexUse *vertexUse) : C3DEntity(C3D_LOOP_ID), pVertexUse(vertexUse), pFatherFace(NULL) {}
	C3DLoop(C3DFace *father, C3DList<C3DVector*> &points);
	~C3DLoop();

	C3DLoop *CopyLoop(C3DOperation &op, C3DEntity *father);
};

class C3DFace : public C3DEntity
{
	friend class C3DShell;
	friend class C3DBoolean;

public:
	// plane on which face lies
	C3DPlane		plane;
	// face bounding box
	C3DBBox			bbox;
	// list of loops in face
	C3DList<C3DLoop*>	loopList;
	// face father: can only be shell
	C3DShell			*pShell;
	// members used for point in shell classification routines
	// in C3DOperation (see C3DFaceBSPTree class)
	C3DFace *testedFace;
	char inOut;
	// members used in face intersection routines
	C3DList<C3DIntersection> intersectList;

	// constructors & destructors
	C3DFace() : C3DEntity(C3D_FACE_ID), pShell(NULL) {}
	C3DFace(C3DShell *father) : C3DEntity(C3D_FACE_ID), pShell(father) {}
	C3DFace(C3DShell *father, C3DList<C3DVector*> &points);
	C3DFace(C3DShell *father, C3DList<C3DVector*> &points, C3DPlane &facePlane);
	~C3DFace();

protected:
	void Triangulate(double minAngle);
	void FlipNormal();

public:
	void EnlargeBBox(C3DBBox enlBbox);
	C3DFace *CopyFace(C3DOperation &op, C3DShell *parentShell);
};

class C3DShell : public C3DEntity
{
	friend class C3DFace;
	friend class C3DBoolean;

public:
	// shell father: can only be another shell
	C3DShell			*pShell;
	// shell bounding box
	C3DBBox				bbox;
	// list of shells in shell
	C3DList<C3DShell*>	shellList;
	// list of faces in shell
	C3DList<C3DFace*>	faceList;
	// list of loops in shell
	C3DList<C3DLoop*>	loopList;
	// list of edges in shell
	C3DList<C3DEdgeUse*>	edgeUseList;
	// vertexuse in shell (if any)
	C3DVertexUse		*pVertexUse;

	// constructors & destructors
	C3DShell() : C3DEntity(C3D_SHELL_ID), pShell(NULL), pVertexUse(NULL) {}
	C3DShell(C3DShell *pFather) : C3DEntity(C3D_SHELL_ID), pShell(pFather), pVertexUse(NULL)
	{
		ASSERT(pShell != NULL);
		// add new subshell to parent shell's shell list
		pShell->shellList.push_front(this);
	}
	~C3DShell()
	{
		Empty();
	}

	void MergeFaceList(C3DList<C3DFace*> &listToBeMerged)
	{
		C3DFace *tmpData;

		while( listToBeMerged.is_empty() == false ) {
			// get element
			tmpData = listToBeMerged.pop_front();
			// change link to parent
			tmpData->pShell = this;
			// and store in local list
			faceList.push_back(tmpData);
		}
	}

	void MergeLoopList(C3DList<C3DLoop*> &listToBeMerged)
	{
		C3DLoop *tmpData;

		while( listToBeMerged.is_empty() == false ) {
			// get element
			tmpData = listToBeMerged.pop_front();
			// change link to parent
			tmpData->pFatherShell = this;
			// and store in local list
			loopList.push_back(tmpData);
		}
	}

	void MergeEdgeUseList(C3DList<C3DEdgeUse*> &listToBeMerged)
	{
		C3DEdgeUse *tmpData;

		while( listToBeMerged.is_empty() == false ) {
			// get element
			tmpData = listToBeMerged.pop_front();
			// change link to parent
			tmpData->pFatherShell = this;
			// and store in local list
			edgeUseList.push_back(tmpData);
		}
	}

	void SetFather(C3DShell *father)
	{
		ASSERT(father != NULL);
		pShell = father;
		father->shellList.push_back(this);
	}

protected:
	// Sphere generation related variables
	C3DShell *m_clsSphereSubshell;
	C3DVector m_clsSphereCenter;
	double m_dSphereRadius;

public:
	void Clear();
	C3DFace *MakeFace(C3DList<C3DVector> &points);
	C3DFace *MakeFace(C3DList<C3DVector*> &points);
	C3DFace *MakeTriFace(C3DVector v1,  C3DVector v2, C3DVector v3);
	C3DFace *MakeFace(C3DVector points[], int numPoints);
	void EnlargeBBox(C3DBBox enlBbox);
	void UpdateBBox();
	void Translate(C3DVector shift);
	void Scale(C3DVector scale);
	void RotateX(double angle);
	void RotateY(double angle);
	void RotateZ(double angle);
	C3DShell *CopyTranslate(C3DVector shift);
	C3DShell *CopyScale(C3DVector scale);
	C3DShell *CopyRotateX(double angle);
	C3DShell *CopyRotateY(double angle);
	C3DShell *CopyRotateZ(double angle);
	void Triangulate(double minAngle);
	void Flatten();
	bool DeleteSubShell(C3DShell *subshell);
	C3DShell *C3DCreateCube(C3DVector &v1, C3DVector &v2);
	C3DShell *C3DCreateCylinder(C3DVector &basePoint, double radius, double height, unsigned int numSides);
	C3DShell *C3DCreateSphere(C3DVector &center, double radius, unsigned int depth);

protected:
	void Empty();
	void Transform(C3DOperation &op);
	C3DShell *CopyTransform(C3DOperation &op, C3DShell *parentShell);
	C3DFace *MakeTriFace(C3DVector v1,  C3DVector v2, C3DVector v3, C3DPlane &facePlane);
	C3DShell *C3DCreateNonTriCube(C3DVector &v1, C3DVector &v2);
	C3DShell *C3DCreateNonTriCylinder(C3DVector &basePoint, double radius, double height, unsigned int numSides);
	void SphereSubdivide(C3DVector &v1, C3DVector &v2, C3DVector &v3, long depth);
	bool RemoveSubShell(C3DShell *subshell);
};
*/

#endif //C3D_GEOMETRY_DEFS

