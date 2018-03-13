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


// Geometry2D.h : definition of geometry-related classes
// Enrico Di Lorenzo 2003/03/12

#ifndef C2D_GEOMETRY_DEFS
#define C2D_GEOMETRY_DEFS

#include <cstddef>
#include <math.h>
#include <vector>

//#include "../Global.h"
#include "GeoGlobal.h"
#include "Vector2D.h"

/////////////////////////////////////////////////////
// 2D objects classes
// these are the classes which describes 2D objects
//


// Binary Space Partition tree definitions
//
#define C2D_BSP_NONE		0
#define C2D_BSP_LEFT		1
#define C2D_BSP_RIGHT		2
#define C2D_BSP_ONPLANE		3

#define C2D_BSP_NEWPOINT	-1

using namespace std;

// declared classes
class C2DBSPTree;


// Cartesian Binary Space Partition tree,
// for maintaining unique vertex positions in space
class C2DBSPTree
{
public:
	// tree node
	class Node
	{
	public:
		char divDir;
		union {
			double divPlane;
			long obj;
		};
		Node *up;
		Node *left;
		Node *right;
		bool child;

		Node()
		{
		}
		Node(Node *upnode, long index) :
		obj(index), up(upnode), left(NULL), right(NULL), child(true)
		{
		}
	};

	// C2DBSPTree iterator
	class iterator;
	friend class iterator;

	class iterator
	{
		friend class C2DBSPTree;
		Node *currentNode;
	public:
		iterator() : currentNode(NULL)
		{
		}
		iterator(C2DBSPTree& starttree) : currentNode(starttree.root)
		{
		}
		iterator(iterator &it) : currentNode(it.currentNode)
		{
		}
		iterator(Node *node) : currentNode(node)
		{
		}
		long operator*() {
			if( currentNode == NULL )
				return 0;

			if(currentNode->child == true)
				return currentNode->obj;
			else
				return 0;
		}
	};

	protected:
		// tree variables
		Node *root;
		double m_dTol;
		vector<C2DVector> *vertexes;

	public:
		C2DBSPTree(vector<C2DVector> *vertexArray) : root(NULL), m_dTol(C2D_TOL), vertexes(vertexArray)
		{
		}

		~C2DBSPTree()
		{
			clear();
		}

		void init(long obj)
		{
			root = new Node;
			root->divDir = C2D_BSP_NONE;
			root->child = true;
			root->up = NULL;
			root->left = NULL;
			root->right = NULL;
			root->obj = obj;
		}

		// Return C2D_BSP_NEWPOINT if a the vertex with the same coordinates
		// of the vertex pointed to by 'index' was not already in the tree
		// (within tolerance), otherwise return the index of the vertex already
		// in the tree
		long insert(long index)
		{
			// if first vertex in tree, must init
			if(root == NULL) {
				init(index);
				// first point is for sure a new point
				return C2D_BSP_NEWPOINT;
			}
			else
				return insert(root, index);
		}

		void clear()
		{
			// already empty
			if( root == NULL)
				return;

			// delete left branch
			erase(root->left);
			// delete right branch
			erase(root->right);

			delete root->left;
			delete root->right;

			// delete root
			delete root;
			root = NULL;
		}

	protected:
		char planeSide(long index, Node *node);
		long insert(Node *node, long index);
		void erase(Node *node);
};



class C2DBBox
{
public:
	bool			empty;
	C2DVector		min_point;
	C2DVector		max_point;

	C2DBBox() : empty(true) {}

	void Clear()
	{
		min_point = C2DVector(0.0, 0.0);
		max_point = C2DVector(0.0, 0.0);
		empty = true;
	}

	unsigned char MaxSide()
	{
		double sidex, sidey;

		sidex = max_point.x - min_point.x;
		sidey = max_point.y - min_point.y;

		if( sidex >= sidey ) {
			return 0;
		}
		else {
			return 1;
		}
	}

	inline void operator=(const C2DBBox &bbox)
	{
		if (bbox.empty == false) {
			min_point = bbox.min_point;
			max_point = bbox.max_point;
			empty = false;
		}
		else
			empty = true;
	}

	inline void operator=(C2DVector &point)
	{
		min_point = point;
		max_point = point;
		empty = false;
	}

	inline void operator+=(C2DBBox &bbox)
	{
		if (bbox.empty == false) {
			if (empty == false) {
				min_point.x = bbox.min_point.x < min_point.x ? bbox.min_point.x : min_point.x;
				min_point.y = bbox.min_point.y < min_point.y ? bbox.min_point.y : min_point.y;
				max_point.x = bbox.max_point.x > max_point.x ? bbox.max_point.x : max_point.x;
				max_point.y = bbox.max_point.y > max_point.y ? bbox.max_point.y : max_point.y;
			}
			else {
				min_point = bbox.min_point;
				max_point = bbox.max_point;
				empty = false;
			}
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
			min_point = point;
			max_point = point;
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
			min_point = point;
			max_point = point;
			empty = false;
		}
	}

	inline C2DBBox operator+(const C2DBBox &bbox)
	{
		C2DBBox newBbox;

		if (bbox.empty == false) {
			if( empty == false) {
				newBbox.min_point.x = bbox.min_point.x < min_point.x ? bbox.min_point.x : min_point.x;
				newBbox.min_point.y = bbox.min_point.y < min_point.y ? bbox.min_point.y : min_point.y;
				newBbox.max_point.x = bbox.max_point.x > max_point.x ? bbox.max_point.x : max_point.x;
				newBbox.max_point.y = bbox.max_point.y > max_point.y ? bbox.max_point.y : max_point.y;
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

	inline C2DBBox operator+(C2DVector &point)
	{
		C2DBBox newBbox;

		if (empty == false) {
			newBbox.min_point.x = point.x < min_point.x ? point.x : min_point.x;
			newBbox.min_point.y = point.y < min_point.y ? point.y : min_point.y;
			newBbox.max_point.x = point.x > max_point.x ? point.x : max_point.x;
			newBbox.max_point.y = point.y > max_point.y ? point.y : max_point.y;
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


#endif //C2D_GEOMETRY_DEFS

