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


// Operation3D.cpp : implementation of the Non-Manifold Geometry basic operations
// Enrico Di Lorenzo 2002/12/30

#include "stdafx.h"

#include <math.h>
#include <vector>
#include <algorithm>
#include <functional>

//#include "../MessageHandler.h"
#include "Operation3D.h"
//#include "../Global.h"

// DEBUG WARNING
//#include "Debug.h"


#define C3D_MAX_FACES_IN_BSP_LEAF	30
#define C3D_MAX_LOOP				1000000

// Static data members must be initialized at file scope, even
// if private.
//
// static initializations:
long C3DOperation::m_lMark = 1;

// Compute plane equation, given three non-collinear points.
// Plane normal direction is point2-point1 x point3-point1
//
//         3
//         |\
//    N    | \
//     \   |  \
//      \  |   \
//       \ |    \
//        \|     \
//         1------2
//
// If plane equation can be found, return C3D_OK and plane 'plane',
// otherwise return error condition

int C3DOperation::PlaneFromThreePoints(C3DPlane &plane, C3DVector &point1,
									   C3DVector &point2, C3DVector &point3)
{
	C3DVector normal;
	double mod;

	normal = CrossProd(point2-point1, point3-point1);
	mod = Mod(normal);

	// if the result of the cross product is too small, 
	// can only means that the points are almost collinear
	if( mod < m_dTol ) {
		return C3D_POINTS_ARE_COLLINEAR;
	}

	// unify normal vector
	normal = normal / mod;

	plane.m_vecNormal = normal;
	plane.m_dD = DotProd(normal, point1);

	return C3D_OK;
}

// Compute plane equation, given non-parallel vectors and a point, 
// all lying on the plane.
// Plane normal direction is vector1 x vector2
//
//             3
//    ^        |\
//    |   N    | \
//    |    \   |  \
//    |     \  |   \
//    |      \ |    \
//    |       \|     \
//    |        1------2
// Vector2
//             -------> Vector1
//
//
// If plane equation can be found, return C3D_OK and plane 'plane',
// otherwise return error condition

int C3DOperation::PlaneFromTwoVectorsAndPoint(C3DPlane &plane, C3DVector &vector1,
									   C3DVector &vector2, C3DVector &point)
{
	C3DVector normal;
	double mod;

	normal = CrossProd(vector1, vector2);
	mod = Mod(normal);

	// if the result of the cross product is too small, 
	// can only means that the vectors are almost parallel
	if( mod < m_dTol ) {
		return C3D_VECTORS_ARE_PARALLEL;
	}

	// unify normal vector
	normal = normal / mod;

	plane.m_vecNormal = normal;
	plane.m_dD = DotProd(normal, point);

	return C3D_OK;
}

// Finds the closest points on two skewed lines
//
// computes point1 on line1 closest to line2
// computes point2 on line2 closest to line1
//
// Note that if lines cross, point1 and point2 can be coincident

int C3DOperation::LineToLineClosestPoints(C3DVector &point1, C3DVector &point2,
										  C3DLine &line1, C3DLine &line2)
{
	C3DVector connline;
	C3DPlane plane1, plane2;
	int ret1, ret2;
	double dist, dist2;

    // connecting line is perpendicular to both
    connline = CrossProd(line1, line2);

	// lines are nearly parallel -> all points are equally close
	if( Mod(connline) < m_dTol ) {
		point1 = line1.m_vecPoint;
		point2 = line2.m_vecPoint;

		return C3D_LINES_ARE_PARALLEL;
	}

    // form plane1 containing line1, parallel to connline
    ret1 = PlaneFromTwoVectorsAndPoint(plane1, connline, line1.m_vecDir, line1.m_vecPoint);

    // form plane2 containing line2, parallel to connline
    ret2 = PlaneFromTwoVectorsAndPoint(plane2, connline, line2.m_vecDir, line2.m_vecPoint);

	// as connline was already checked, this condition should never happen
	if( ret1 != C3D_OK || ret2 != C3D_OK ) {
		return C3D_INTERNAL_ERROR;
	}

    // closest point on line1 is obtained intersecting line1 with plane2
	ret1 = IntersLinePlane(point1, line1, plane2, dist, dist2);

     // closest point on line2 is obtained intersecting line2 with plane1
	ret2 = IntersLinePlane(point2, line2, plane1, dist, dist2);

	// as connline was already checked, this condition should never happen
	if( ret1 != C3D_OK || ret2 != C3D_OK ) {
		return C3D_INTERNAL_ERROR;
	}

	return C3D_OK;
}

char C3DOperation::IsPointInFace(C3DVector &point, C3DFace *face)
{
	char i, j;
	C3DList<C3DLoop*>::iterator itl;
	C3DList<C3DEdgeUse*>::iterator iteu;
	C3DVector v1, v2;
	double x1, y1, x2, y2;
	char quadrant1, quadrant2, angle;
	int totalAngle;

	// find the axis plane on which the face projects the biggest area 
	FindProjPlane(&i, &j, face->plane.m_vecNormal);

	totalAngle = 0;

	// for every loop in face
	for(itl = face->loopList.begin(); itl != face->loopList.end(); itl++) {
		
		// get first vertex of first edge
		iteu = (*itl)->edgeUseList.begin();
		v1 = (*iteu)->GetFirstVertexUse()->pVertex->pos;
		x1 = v1[i] - point[i];
		y1 = v1[j] - point[j];

		// for every edgeuse in loop
		for(iteu = (*itl)->edgeUseList.begin(); iteu != (*itl)->edgeUseList.end(); iteu++) {
			
			// get new edge vertex
			v2 = (*iteu)->GetSecondVertexUse()->pVertex->pos;
			x2 = v2[i] - point[i];
			y2 = v2[j] - point[j];
		
			// find vertex quadrant position
			quadrant1 = FindQuadrant(x1, y1);
			quadrant2 = FindQuadrant(x2, y2);

			// compute angle
			angle = ComputeAngle(quadrant1, quadrant2, x1, y1, x2, y2);

			// if point is on edge or on vertex, end here
			if( angle == C3D_POINT_ON_EDGE )
				return C3D_POINT_ON_FACE_EDGE;
			if( angle == C3D_POINT_ON_VERTEX)
				return C3D_POINT_ON_FACE_VERTEX;

			// sum up angles
			totalAngle += angle;

			// and move to next point
			x1 = x2;
			y1 = y2;
		}
	}

	if( totalAngle == 0 )
		return C3D_POINT_OUT_FACE;
	else
		return C3D_POINT_IN_FACE;
}

// Remark: it is assumed that the shell is 'flat' (no subshells)
void C3DFaceBSPTree::InitPointInShell(C3DShell &shell, double tolerance)
{
	C3DList<C3DFace*>::iterator itf;

	// check that the shell is 'flat'
	ASSERT(shell.shellList.begin() == shell.shellList.end());

	// set tolerance
	m_dTol = tolerance;

	// delete the old tree, if any
	clear();

	for(itf = shell.faceList.begin(); itf != shell.faceList.end(); itf++) {
		// clear face references and inOut statuses, as required 
		// by C3DFaceBSPTree's insert() routine
		(*itf)->testedFace = NULL;
		(*itf)->inOut = C3D_FACEBSP_NONE;

		// insert face in BSP tree
		insert(*itf);
	}
}

// Insert faces recursively into the tree
// Warinig: no algorithm to balance tree is used
void C3DFaceBSPTree::insert(Node *node, C3DFace *pFace)
{
	char inOut;
	C3DFace *nodeFace;
	C3DList<C3DLoop*>::iterator itl;
	C3DList<C3DEdgeUse*>::iterator iteu;
	C3DVector point;
	C3DPlane nodePlane;
	double dist;
	C3DOperation op;
				
	nodeFace = node->face;
	nodePlane = nodeFace->plane;
				
	// clear inOut status
	inOut = C3D_FACEBSP_ON;
				
	// if current face 'pFace' have already been tested
	// against the face stored in current node, don't
	// waste time performing the test again
	// (can happen when inserting face both inside and outside)
	if(nodeFace->testedFace == pFace)
		inOut = nodeFace->inOut;
	else {
		// scan all current face 'pFace' vertexes 
		// and test them against the face stored in current node
		for(itl = pFace->loopList.begin(); itl != pFace->loopList.end(); itl++) {
			for(iteu = (*itl)->edgeUseList.begin(); iteu != (*itl)->edgeUseList.end(); iteu++) {
				
				// get edge first end point
				point = (*iteu)->GetFirstVertexUse()->pVertex->pos;
							

					// DEBUG WARNING
//				int figlio = 0;
//					if( Mod(point - C3DVector(-0.618, 1.618, -2.2204)) < 0.1 )
//						figlio++;

					
					
				// calculate the distance from the end point to the face plane
				dist = op.PointPlaneDist(point, nodePlane);
							
				// if distance is negative (within tolerance), point is inside;
				// if distance is positive (within tolerance), point is outside;
				// otherwise is on. 
				// Note that the following combinations are possible, giving
				// rise to the following cases: 
				// on + inside = inside, on + outside = outside, 
				// inside + on + outside = inout,  on = on
				if( dist < -m_dTol )
					inOut |= C3D_FACEBSP_INSIDE;
				else if( dist > m_dTol )
					inOut |= C3D_FACEBSP_OUTSIDE;
				else
					inOut |= C3D_FACEBSP_ON;

			}
		}
					
		// store pFace and inOut status in the face 
		// stored in current node
		nodeFace->testedFace = pFace;
		nodeFace->inOut = inOut;
	}
				
	// recurse through the tree
	// (remark: if faces are coplanar, go inside; this is arbitrary, but find()
	// must work accordingly)
	if( (inOut & C3D_FACEBSP_INSIDE) != 0 || inOut == C3D_FACEBSP_ON ) {
		// increment counter of nodes in 'inside' branch
		node->num_inside++;
		// if no 'inside' node, create new one, otherwise go
		// deeper in the tree
		if(node->inside == NULL)
			node->inside = new Node(pFace);
		else
			insert(node->inside, pFace);
	}
	if( (inOut & C3D_FACEBSP_OUTSIDE) != 0 ) {
		// increment counter of nodes in 'inside' branch
		node->num_outside++;
		// if no 'inside' node, create new one, otherwise go
		// deeper in the tree
		if(node->outside == NULL)
			node->outside = new Node(pFace);
		else
			insert(node->outside, pFace);
	}	
}

char C3DFaceBSPTree::find(Node *node, C3DVector &point, FILE *fp)
{
	double dist;
	char inface, res;
	C3DOperation op;
	
	// calculate the distance from the point to the face plane
	// for face in current node
	dist = op.PointPlaneDist(point, node->face->plane);

	if(fp != NULL) {
		fprintf(fp, "* dist = %f\r\n", dist);
		OutputFace(node->face, fp);
	}

	// if distance is negative (within tolerance)
	if( dist < -m_dTol ) {
		// if no more 'inside', the point itself is inside the shell
		if( node->inside == NULL )
			return C3D_POINT_IN_SHELL;
		else {
			res = find(node->inside, point, fp);

			// DEBUG WARNING
//			ScanFace(node->face);


			if( res == C3D_POINT_UNDECIDABLE )
				return C3D_POINT_IN_SHELL;
			else
				return res;
		}

	}
	// if distance is positive (within tolerance)
	else if( dist > m_dTol ) {
		// if no more 'outside', the point itself is outside the shell
		if( node->outside == NULL ) {
			// DEBUG WARNING
//			ScanFace(node->face);

			return C3D_POINT_OUT_SHELL;
		}
		else {
			res = find(node->outside, point, fp);

			// DEBUG WARNING
//			ScanFace(node->face);

			
			if( res == C3D_POINT_UNDECIDABLE )
				return C3D_POINT_OUT_SHELL;
			else
				return res;
		}
	}
	// in this case the point lies on the face plane (within tolerance),
	// so must check if it is in the face, on its perimeter or ouside
	else {
		inface = op.IsPointInFace(point, node->face);

		if(inface == C3D_POINT_IN_FACE)
			return C3D_POINT_ON_SHELL_FACE;

		if(inface == C3D_POINT_ON_FACE_EDGE)
			return C3D_POINT_ON_SHELL_EDGE;

		if(inface == C3D_POINT_ON_FACE_VERTEX)
			return C3D_POINT_ON_SHELL_VERTEX;

		// If none of the above, point is outside given face
		ASSERT(inface == C3D_POINT_OUT_FACE);

		// In this case, we don't know if the point is inside or
		// outside shell. So further testing is needed: 

		// Go inside to test, among the other things, for coplanar faces;
		// this is because the BSP tree face insertion algorithm  
		// insert a new face only inside, if two faces are coplanar
		// (Remark: the opposite choice could have been made, that is,
		// insert coplanar faces both in and out, then choose the 
		// search direction for the point find operation according to
		// the branch with less nodes.)
		if(node->inside != NULL) {
			res = find(node->inside, point, fp);
		
			// if found an answer, return it
			if( res != C3D_POINT_UNDECIDABLE )
				return res;
		}

		// otherwise, return whatever found in outside branch (last chance)
		if(node->outside != NULL)
			return (find(node->outside, point, fp));

		// If no more branches, or still undecidable, the situation
		// is undecidable at this point and must be solved upward in the tree,
		// so simply report a dead branch
		return C3D_POINT_UNDECIDABLE;

	}
}

// Remark: it is assumed that the shell is 'flat' (no subshells)
void C3DCartBSPTree::InitBSPTree(C3DShell &shell, double tolerance)
{
	C3DList<C3DFace*>::iterator itf;

	// check that the shell is 'flat'
	ASSERT(shell.shellList.begin() == shell.shellList.end());

	// set tolerance
	m_dTol = tolerance;

	// delete the old tree, if any
	clear();

	// create new root node
	root = new Node;

	// copy the list of faces into the root
	root->faceList.Copy(shell.faceList);

	// start recursive distribution of faces into the BSP tree
	insert(root, shell.bbox.max_point - shell.bbox.min_point);

	// store shell bbox
	globalBbox = shell.bbox;
}

// Insert faces recursively into the tree
// Uses algorithm to balance tree
void C3DCartBSPTree::insert(Node *pNode, C3DVector bboxext)
{
	long numNodes, i;
	C3DList<C3DFace*>::iterator it;
	std::pair<double, double> minMax;
	double min, max;
	C3DVector bboxIn(bboxext), bboxOut(bboxext);

	numNodes = pNode->faceList.size();

	// if the node already contains less than C3D_MAX_FACES_IN_BSP_LEAF faces,
	// it is a leaf node, therefore return
	if(numNodes <= C3D_MAX_FACES_IN_BSP_LEAF)
		return;

	// declare vectors for splitting plane calculation
	std::vector<double> minv, maxv;
	minv.reserve(numNodes);
	maxv.reserve(numNodes);

	// define split direction
	pNode->divDir = m_clsOp.FindLargestProj(bboxext);

	//
	// calculate the split value
	//

	// get min and max extent of faces along chosen split direction
	for(it = pNode->faceList.begin(), i=0; it != pNode->faceList.end(); it++, i++) {
		minMax = getMinMax(*it, pNode->divDir);
		minv.push_back(minMax.first);
		maxv.push_back(minMax.second);
	}

	// new vectors for sorting (to preserve the old ones)
	std::vector<double> minvSort(minv), maxvSort(maxv);

	// sort min[] in ascending order, sort max in descending order
	std::sort(minvSort.begin(), minvSort.end(), std::less<double>());
	std::sort(maxvSort.begin(), maxvSort.end(), std::greater<double>());

	// initialize split point, in case split point is not found 
	// (e.g. only one primitive in list, min[0] is never >= max[0],
	// or two primitives with equal mins and maxes)
	// To this goal, let's take as split plane the minimum maximum
	// (could also choose the maximum minimum) so at least one
	// face will NOT be split by the plane, therefore the agorithm
	// is guaranteed to terminate in the end.
	pNode->divPlane = maxvSort[numNodes-1];
	// search for the split point
	for(i=0; i<numNodes; i++) {
		if(minvSort[i] >= maxvSort[i]) {
			if(i != 0) {
				max = max(maxvSort[i], minvSort[i-1]);
				min = min(minvSort[i], maxvSort[i-1]);
			}
			else {
				max = maxvSort[0];
				min = minvSort[0];
			}

			pNode->divPlane = (max + min) / 2.0;
			break;
		}
	}

	// create child nodes
	pNode->inside = new Node;
	pNode->outside = new Node;

	// and distribute the faces in the children according to the split plane
	for(it = pNode->faceList.begin(), i=0; it != pNode->faceList.end(); it++, i++) {
		// inside
		// Remark: test using tolerance is made only when *searching*
		// the tree, assuming that when in proximity of a splitting
		// plane, the face we are searching for could have been inserted
		// also in the other subtree
		if(maxv[i] <= pNode->divPlane) {
			pNode->inside->faceList.push_back(*it);
		}
		// outside
		else if(minv[i] >= pNode->divPlane) {
			pNode->outside->faceList.push_back(*it);
		}
		// crossing, insert in both
		else {
			pNode->inside->faceList.push_back(*it);
			pNode->outside->faceList.push_back(*it);
		}
	}

	// clear face list of this node
	pNode->faceList.clear();

	//
	// recursively build tree
	//

	max = maxvSort[0];
	min = minvSort[0];
	bboxIn[pNode->divDir] = (pNode->divPlane - min);
	bboxOut[pNode->divDir] = (max - pNode->divPlane);

	// recurse into subtree only if the number of faces
	// in the subtree is smaller than the number of faces
	// in the current node, to avoid infinite recursion
	// because of faces split by the split plane
	// (the bbox mechanism guarantees that we tried 
	// the best splitting direction before giving up partition)
	if( pNode->inside->faceList.size() < numNodes )
		insert(pNode->inside, bboxIn);
	if( pNode->outside->faceList.size() < numNodes )
		insert(pNode->outside, bboxOut);
}

// Containment test, using ray casting approach: a ray is fired
// in a random direction and the number of intersections with 
// the faces is counted. If this number is odd, point is inside shell.
// If the ray hits a vertex or an edge, the ray is considered invalid
// and another ray is fired (to avoid handling the special cases).
// Of course the starting point is tested also against being ON the shell.
char C3DCartBSPTree::IsPointInShell(C3DVector &point)
{
	C3DVector dforw, dback, midpoint;
	std::pair<char, long> res;
	long i;
	double maxdist;

	// verify if point is inside BSP global bbox
	if(m_clsOp.IsPointInBBox(globalBbox, point) == false)
		return C3D_POINT_OUT_SHELL;

	m_clsRay.m_vecPoint = point;

	res.first = C3D_POINT_UNINITIALIZED;
	i = 0;
	do {
		// choose a ray direction 
		
		// init random seed
		srand((unsigned)time(NULL)); 
		
		// generate random number (either positive or negative)
		m_clsRay.m_vecDir.x = (double)(rand() - RAND_MAX / 2.0);
		m_clsRay.m_vecDir.y = (double)(rand() - RAND_MAX / 2.0);
		m_clsRay.m_vecDir.z = (double)(rand() - RAND_MAX / 2.0);
		
		m_clsRay.m_vecDir.Normalize();


		// DEBUG
/*		FILE *fp = fopen("CubeFaces.txt", "w");
		fprintf(fp, "0 NMG test file '%s'\r\n", "CubeFaces.txt");
		fprintf(fp, "T  Face  %f %f %f  %f %f %f  %f %f %f\r\n", 
					m_clsRay.m_vecPoint.x, m_clsRay.m_vecPoint.y, m_clsRay.m_vecPoint.z,
					m_clsRay.m_vecPoint.x, m_clsRay.m_vecPoint.y, m_clsRay.m_vecPoint.z,
					m_clsRay.m_vecPoint.x + m_clsRay.m_vecDir.x * 4.0,
					m_clsRay.m_vecPoint.y + m_clsRay.m_vecDir.y * 4.0,
					m_clsRay.m_vecPoint.z + m_clsRay.m_vecDir.z * 4.0);
		fclose(fp);
*/



		// discard ray if parallel to any cartesian plane (that is the
		// planes making up the cartesian BSP)
		if( m_clsRay.m_vecDir.x > -m_dTol && m_clsRay.m_vecDir.x < m_dTol )
			continue;
		if( m_clsRay.m_vecDir.y > -m_dTol && m_clsRay.m_vecDir.y < m_dTol )
			continue;
		if( m_clsRay.m_vecDir.z > -m_dTol && m_clsRay.m_vecDir.z < m_dTol )
			continue;

		// get a unique marker to mark faces, since the same face
		// could appear more than one time if was intersected
		// by the split plane (and so inserted in both inside
		// and outside tree branches)
		m_lMark = m_clsOp.GetMarker();
		
		// fire the ray and recursively find face intersections
		// Remark: the maximum search distance for the ray is:
		// - if ray origin is outside bbox of the shell, the distance
		//   of the ray origin to the center of the bbox, 
		//   multiplied by two to mirror it 
		// - if the ray origin is inside the bbox of the shell,
		//   the diagonal of the bbox
		midpoint = (globalBbox.max_point + globalBbox.min_point) / 2.0;
		maxdist = max(Mod(m_clsRay.m_vecPoint - midpoint) * 2.0, Mod(globalBbox.max_point - globalBbox.min_point));
		res = fireRay(root, 0.0, maxdist + m_dTol);
//		DEBUG
//		res = visitTree(root, 0.0, maxdist, 0);
		i++;



		// TBC warning DEBUG
//		if(res.first == C3D_RAY_INVALID)
//			i++;




	} while (res.first == C3D_RAY_INVALID && i < C3D_MAX_LOOP);

	if( i >= C3D_MAX_LOOP )
		ASSERT(FALSE);
	
	// if point was not ON the shell boundary,
	// count the number of face intersections
	if(res.first == C3D_POINT_UNINITIALIZED) {
		// if even number of intersections
		if( (res.second % 2) == 0 )
			return C3D_POINT_OUT_SHELL;
		else
			return C3D_POINT_IN_SHELL;
	}
	else
		return res.first;
}

std::pair<char, long> C3DCartBSPTree::fireRay(Node *pNode, double min_dist, double max_dist)
{
	Node *farNode, *nearNode;
	std::pair<char, long> res1, res2;
	int status;
	C3DList<C3DFace*>::iterator it;
	char inface;
	C3DVector intersection, planenormal;
	double dist, orthdist;

	// if child node, test every face in list for intersection
	// and return either the number of intersections or 
	// (in case of ray origin on a face) the status
	// (on vertex, on edge or simply on face)
	if( pNode->divDir == C3D_DIR_NONE ) {

		res1 = std::make_pair((char)C3D_POINT_UNINITIALIZED, (long)0);



		// DEBUG
//		FILE *fp = fopen("CubeFaces.txt", "a+");
//		fprintf(fp, "*  Leaf node\r\n");



		for(it = pNode->faceList.begin(); it != pNode->faceList.end(); it++) {
			

			// DEBUG WARNING
			//ScanFace(*it);


			// if face already visited, skip
			if( (*it)->index == m_lMark )
				continue;

			// mark face as visited
			(*it)->index = m_lMark;




			// DEBUG
//			OutputFace(*it, fp);





			// intersect ray and face support plane
			status = m_clsOp.IntersLinePlane(intersection, m_clsRay, (*it)->plane, dist, orthdist);
			
			// if ray is parallel but not on the face plane, go on
			if( status == C3D_LINE_DONT_INTERSECT )
				continue;
			// if ray lies on the face plane, must check for no intersection,
			// otherwise ray is considered invalid, to avoid special case handling 
			if( status == C3D_LINE_IS_ON_PLANE ) {
				// TBC warning: to be completed
				ASSERT(false);
				res1.first = C3D_RAY_INVALID;
				break;
			}

			// test if the intersection point is coincident with the ray origin
			if( orthdist < m_dTol && orthdist > -m_dTol) {
				
				// test *ray origin* containment: if 'intersection' were tested
				// for containment, we could fail in case the ray is almost 
				// parallel to the plane: the intersection can lie far from 
				// the origin, even if 'orthdist' is very small
				inface = m_clsOp.IsPointInFace(m_clsRay.m_vecPoint, *it);

				if(inface == C3D_POINT_IN_FACE) {
					res1.first = C3D_POINT_ON_SHELL_FACE;
					break;
				}
				
				if(inface == C3D_POINT_ON_FACE_EDGE) {
					res1.first = C3D_POINT_ON_SHELL_EDGE;
					break;
				}
				
				if(inface == C3D_POINT_ON_FACE_VERTEX) {
					res1.first = C3D_POINT_ON_SHELL_VERTEX;
					break;
				}
				
				// otherwise is out of face (so not on shell boundary),
				// therefore continue with ray intersection
				continue;
			}

			// check if point of intersection is valid, that is, on the positive
			// direction along ray; if not, go on
			// (it has been already verified that ray is not parallel to the face plane)
			if( dist < -m_dTol )
				continue;

			// ray is not parallel to the face plane, the intersection point
			// is not the ray origin and the intersection is valid
			//

			// test point containment
			inface = m_clsOp.IsPointInFace(intersection, *it);
			
			// the intersection point is not the ray origin 
			if(inface == C3D_POINT_IN_FACE) {
				// increment the intersection count, then go on
				// scanning the faces in the node (there may be 
				// more than one intersection)
				res1.second++;
			}
			// therefore if it happens to be on a face edge or vertex, 
			// must discard ray to avoid handing of special cases 
			else if(inface == C3D_POINT_ON_FACE_EDGE || inface == C3D_POINT_ON_FACE_VERTEX) {
				res1.first = C3D_RAY_INVALID;
				break;
			}
			
			// otherwise is out of face (so not on shell boundary),
			// therefore continue with ray intersection
		}



			// DEBUG
//			fclose(fp);


		return res1;
	}
	// in this case the node is not a leaf, so recurse in the tree
	else {

		// compute division plane
		planenormal[pNode->divDir] = 1.0;
		C3DPlane divplane(planenormal, pNode->divPlane);

		// intersect ray and division plane to find distance of ray origin to split plane
		// along ray direction
		status = m_clsOp.IntersLinePlane(intersection, m_clsRay, divplane, dist, orthdist);
			
		if( status != C3D_OK )
			ASSERT(false);

		// determine near and far children, according to ray origin halfplane position
		if( orthdist > 0.0) {
			nearNode = pNode->inside;
			farNode = pNode->outside;
		}
		else {
			farNode = pNode->inside;
			nearNode = pNode->outside;
		}
		

		// DEBUG
/*		FILE *fp = fopen("CubeFaces.txt", "a+");
		fprintf(fp, "* Tree node\r\n");
		fprintf(fp, "* divDir %d, divPlane %f, dist %f, orthdist %f, min_dist %f, max_dist %f\r\n",
					pNode->divDir, (float)pNode->divPlane, (float)dist, (float)orthdist, (float)min_dist, (float)max_dist);
		if( pNode->divDir == 0 ) {
			fprintf(fp, "*Q  Face  %f %f %f  %f %f %f  %f %f %f %f %f %f\r\n", 
					pNode->divPlane, -4.0, -4.0,
					pNode->divPlane, -4.0, +4.0,
					pNode->divPlane, +4.0, +4.0,
					pNode->divPlane, +4.0, -4.0);
		}
		else if( pNode->divDir == 1 ) {
			fprintf(fp, "*Q  Face  %f %f %f  %f %f %f  %f %f %f %f %f %f\r\n", 
					-4.0, pNode->divPlane, -4.0,
					-4.0, pNode->divPlane, +4.0,
					+4.0, pNode->divPlane, +4.0,
					+4.0, pNode->divPlane, -4.0);
		}
		else {
			fprintf(fp, "*Q  Face  %f %f %f  %f %f %f  %f %f %f %f %f %f\r\n", 
					-4.0, -4.0, pNode->divPlane,
					-4.0, +4.0, pNode->divPlane,
					+4.0, +4.0, pNode->divPlane,
					+4.0, -4.0, pNode->divPlane);
		}
		fclose(fp);
*/


		// if ray origin is far enough from the split plane, 'dist' can
		// be considered valid, so can be used in comparisons
		if( orthdist > m_dTol || orthdist < -m_dTol) {	
			// if 'dist' is negative or greater than max_dist, then whole interval is on near side
			// (that is, the half plane containing the ray origin)
			if( dist < -m_dTol || dist > max_dist + m_dTol)
				return fireRay(nearNode, min_dist, max_dist);
			
			// if 'dist' is smaller than 'min_dist', then whole interval is on far side
			// (that is, the half plane not containing the ray origin)
			if( dist < min_dist - m_dTol)
				return fireRay(farNode, min_dist, max_dist);

			// otherwise ray crosses both children
			res1 = fireRay(nearNode, min_dist, dist + m_dTol);
			// if ray origin is not 'on' or ray is not invalid
			if( res1.first == C3D_POINT_UNINITIALIZED ) {
				res2 = fireRay(farNode, max(dist - m_dTol, 0.0), max_dist);
				// sum up the number of intersections from children
				res1.second += res2.second;
				// and report 'on' case, if has happened
				res1.first = res2.first;
			}
		}
		// otherwise ray must be regarded as crossing both children;
		// also, cannot safely split min_dist-max_dist segment at 'dist' 
		else {
			res1 = fireRay(nearNode, min_dist, max_dist);
			// if ray origin is not 'on' or ray is not invalid
			if( res1.first == C3D_POINT_UNINITIALIZED ) {
				res2 = fireRay(farNode, min_dist, max_dist);
				// sum up the number of intersections from children
				res1.second += res2.second;
				// and report 'on' case, if has happened
				res1.first = res2.first;
			}
		}

		return res1;
	}
}


// DEBUG
std::pair<char, long> C3DCartBSPTree::visitTree(Node *pNode, double min_dist, double max_dist, int level)
{
	Node *farNode, *nearNode;
	std::pair<char, long> res1, res2;
	int status;
	C3DList<C3DFace*>::iterator it;
	char inface;
	C3DVector intersection, planenormal;
	double dist, orthdist;

	// if child node, test every face in list for intersection
	// and return either the number of intersections or 
	// (in case of ray origin on a face) the status
	// (on vertex, on edge or simply on face)
	if( pNode->divDir == C3D_DIR_NONE ) {

		res1 = std::make_pair((char)C3D_POINT_UNINITIALIZED, (long)0);



		// DEBUG
//		FILE *fp = fopen("CubeFaces.txt", "a+");
//		fprintf(fp, "*  Leaf node\r\n");



		for(it = pNode->faceList.begin(); it != pNode->faceList.end(); it++) {
			

			// DEBUG WARNING
			//ScanFace(*it);

			// DEBUG
			// if face already visited, comment
//			if( (*it)->index == m_lMark )
//				fprintf(fp, "**");
//			OutputFace(*it, fp);

			// if face already visited, skip
			if( (*it)->index == m_lMark )
				continue;

			// mark face as visited
			(*it)->index = m_lMark;

			// intersect ray and face support plane
			status = m_clsOp.IntersLinePlane(intersection, m_clsRay, (*it)->plane, dist, orthdist);
			
			// if ray is parallel but not on the face plane, go on
			if( status == C3D_LINE_DONT_INTERSECT )
				continue;
			// if ray lies on the face plane, must check for no intersection,
			// otherwise ray is considered invalid, to avoid special case handling 
			if( status == C3D_LINE_IS_ON_PLANE ) {
				// TBC warning: to be completed
				ASSERT(false);
				res1.first = C3D_RAY_INVALID;
				break;
			}

			// test if the intersection point is coincident with the ray origin
			if( orthdist < m_dTol && orthdist > -m_dTol) {
				
				// test *ray origin* containment: if 'intersection' were tested
				// for containment, we could fail in case the ray is almost 
				// parallel to the plane: the intersection can lie far from 
				// the origin, even if 'orthdist' is very small
				inface = m_clsOp.IsPointInFace(m_clsRay.m_vecPoint, *it);

				if(inface == C3D_POINT_IN_FACE) {
					res1.first = C3D_POINT_ON_SHELL_FACE;
					break;
				}
				
				if(inface == C3D_POINT_ON_FACE_EDGE) {
					res1.first = C3D_POINT_ON_SHELL_EDGE;
					break;
				}
				
				if(inface == C3D_POINT_ON_FACE_VERTEX) {
					res1.first = C3D_POINT_ON_SHELL_VERTEX;
					break;
				}
				
				// otherwise is out of face (so not on shell boundary),
				// therefore continue with ray intersection
				continue;
			}

			// check if point of intersection is valid, that is, on the positive
			// direction along ray; if not, go on
			// (it has been already verified that ray is not parallel to the face plane)
			if( dist < -m_dTol )
				continue;

			// ray is not parallel to the face plane, the intersection point
			// is not the ray origin and the intersection is valid
			//

			// test point containment
			inface = m_clsOp.IsPointInFace(intersection, *it);
			
			// the intersection point is not the ray origin 
			if(inface == C3D_POINT_IN_FACE) {
				// increment the intersection count, then go on
				// scanning the faces in the node (there may be 
				// more than one intersection)
				res1.second++;
			}
			// therefore if it happens to be on a face edge or vertex, 
			// must discard ray to avoid handing of special cases 
			else if(inface == C3D_POINT_ON_FACE_EDGE || inface == C3D_POINT_ON_FACE_VERTEX) {
				res1.first = C3D_RAY_INVALID;
				break;
			}
			
			// otherwise is out of face (so not on shell boundary),
			// therefore continue with ray intersection
		}



			// DEBUG
//			fclose(fp);


		return res1;
	}
	// in this case the node is not a leaf, so recurse in the tree
	else {

		// compute division plane
		planenormal[pNode->divDir] = 1.0;
		C3DPlane divplane(planenormal, pNode->divPlane);

		// intersect ray and division plane to find distance of ray origin to split plane
		// along ray direction
		status = m_clsOp.IntersLinePlane(intersection, m_clsRay, divplane, dist, orthdist);
			
		if( status != C3D_OK )
			ASSERT(false);

		// determine near and far children, according to ray origin halfplane position
		if( orthdist > 0.0) {
			nearNode = pNode->inside;
			farNode = pNode->outside;
		}
		else {
			farNode = pNode->inside;
			nearNode = pNode->outside;
		}
/*		

		// DEBUG
		FILE *fp = fopen("CubeFaces.txt", "a+");
		fprintf(fp, "* Tree node\r\n");
		fprintf(fp, "* divDir %d, divPlane %f, dist %f, orthdist %f, level %d\r\n",
					pNode->divDir, (float)pNode->divPlane, (float)dist, (float)orthdist, level);
		if( pNode->divDir == 0 ) {
			fprintf(fp, "*Q  Face  %f %f %f  %f %f %f  %f %f %f %f %f %f\r\n", 
					pNode->divPlane, -4.0, -4.0,
					pNode->divPlane, -4.0, +4.0,
					pNode->divPlane, +4.0, +4.0,
					pNode->divPlane, +4.0, -4.0);
		}
		else if( pNode->divDir == 1 ) {
			fprintf(fp, "*Q  Face  %f %f %f  %f %f %f  %f %f %f %f %f %f\r\n", 
					-4.0, pNode->divPlane, -4.0,
					-4.0, pNode->divPlane, +4.0,
					+4.0, pNode->divPlane, +4.0,
					+4.0, pNode->divPlane, -4.0);
		}
		else {
			fprintf(fp, "*Q  Face  %f %f %f  %f %f %f  %f %f %f %f %f %f\r\n", 
					-4.0, -4.0, pNode->divPlane,
					-4.0, +4.0, pNode->divPlane,
					+4.0, +4.0, pNode->divPlane,
					+4.0, -4.0, pNode->divPlane);
		}
		fclose(fp);
*/

		level++;

		res1 = visitTree(nearNode, min_dist, max_dist, level);
		res2 = visitTree(farNode, min_dist, max_dist, level);
					
		return res1;
	}
}

