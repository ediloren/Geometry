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


// Geometry2D.cpp : implementation of the geometry-related classes
// Enrico Di Lorenzo 2003/07/25

#include "stdafx.h"

#include "Geometry2D.h"

// find on which side of plane in 'node' the vertex pointed to by 'index' lies
char C2DBSPTree::planeSide(long index, Node *node)
{
	double var, dist;

	if( node->divDir == C2D_DIR_X )
		var = (*vertexes)[index].x;
	else if( node->divDir == C2D_DIR_Y )
		var = (*vertexes)[index].y;

	dist = var - node->divPlane;

	if( fabs(dist) * 2.0 < m_dTol )
		return C2D_BSP_ONPLANE;

	if( dist > 0 )
		return C2D_BSP_RIGHT;
	else
		return C2D_BSP_LEFT;
}

// Insert vertex pointed to by 'index' in BSP tree
long C2DBSPTree::insert(Node *node, long index)
{
	char side;
	long present;
	double fdx, fdy, dx, dy;
	double 	newVertexX, newVertexY, localVertexX, localVertexY;
	Node *node1, *node2;

	// if insertion point not null
	if(node != NULL) {
		// if not a child
		if(node->child == false) {
			// decide on which side the point lies
			side = planeSide(index, node);
			if( side == C2D_BSP_LEFT) {
				return insert(node->left, index);
			}
			else if(side == C2D_BSP_RIGHT) {
				return insert(node->right, index);
			}
			else {
				// point is on the plane, insert in both left and right
				// branches. But if point is already in tree (seen in
				// left branch) don't bother following also right branch
				present = insert(node->left, index);
				if(present == C2D_BSP_NEWPOINT)
					return insert(node->right, index);
				else
					return present;
			}
		}
		// if node is a child
		else {
			newVertexX = (*vertexes)[index].x;
			newVertexY = (*vertexes)[index].y;
			localVertexX = (*vertexes)[node->obj].x;
			localVertexY = (*vertexes)[node->obj].y;

			dx = newVertexX - localVertexX;
			dy = newVertexY - localVertexY;

			fdx = fabs(dx);
			fdy = fabs(dy);

			// if point already exist
			if(fdx < m_dTol && fdy < m_dTol)
				return node->obj;

			// otherwise insert it

			node1 = new Node(node, node->obj);
			node2 = new Node(node, index);

			node->child = false;

			if(fdx > fdy) {
				node->divDir = C3D_DIR_X;
				node->divPlane = (newVertexX + localVertexX) / 2.0;
				if( dx > 0 ) {
					node->left = node1;
					node->right = node2;
				}
				else {
					node->right = node1;
					node->left = node2;
				}
			}
			else {
				node->divDir = C3D_DIR_Y;
				node->divPlane = (newVertexY + localVertexY) / 2.0;
				if( dy > 0 ) {
					node->left = node1;
					node->right = node2;
				}
				else {
					node->right = node1;
					node->left = node2;
				}
			}
		}
	}

	return C2D_BSP_NEWPOINT;
}

void C2DBSPTree::erase(Node *node)
{
	if(node->child == true)
		return;

	erase(node->left);
	erase(node->right);

	delete node->left;
	delete node->right;
}
