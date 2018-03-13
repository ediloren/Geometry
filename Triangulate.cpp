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


// Triangulate.cpp : implementation of Delaunay triangulation routines
// Enrico Di Lorenzo 2003/05/05

#include "stdafx.h"

// for binary_function
#include <functional>

#include "Triangulate.h"
#include "MathDefs.h"
#include "Geometry2D.h"

// max # of iterations in loop
#define C2D_MAX_LOOP	1000000

// For sorting vertexes must use this operator function;
// as a matter of fact, the 'sortVertexes' contains only
// references to C2DVector array of vertex points. Therefore,
// must order the array according to the C2DVectors pointed to;
// this means that this operator must access the 'vertexes' array.
// I did not find how to access directly the 'vertexes' array
// in class C2DTriangulate (a C2DTriangulate::LessThanV() function
// did not work), so this struct must be initialized with a
// reference to the 'vertexes' array
//
// Remark 1: the tests using C3D_TOL have been added to prevent wrong
//  sorting of points within tolerance (i.e. if points are coincident
//  within tolerance but slightly offset, if we use the strict '<' operator
//  sorting is ok, but subsequent triangulation, which uses tolerance, may fail).
//
//  Warning: this sort routine could fail to give the required
//  vertex order because of tolerances in input points
//  (e.g. if there are three points, say 1, 2, 3 and 1,2 are within C3D_TOL,
//  2,3 too but 1,3 are more distant (though < 2*C3D_TOL) )
//
// Remark 2: no C3D_TOL test is needed on 'y' coordinate, because
//  coincident points within tolerance have already been removed
//  from the input
//
// Remark 3: in case this routine is changed, also LessThanVdirty() must be
//  changed acordingly.
struct LessThanV  : public binary_function<C2DVertex &, C2DVertex &, bool>
{
	LessThanV(vector<C2DVector> &varray)
	{
		vertexes = &varray;
	}

	// compare vertex positions according to 'pointed to' vertex
	bool operator()(const C2DVertex &p1, const C2DVertex &p2)
	{
		if( (*vertexes)[p1.sortedVtoV].x + C3D_TOL < (*vertexes)[p2.sortedVtoV].x ) {
			return true;
		}
		else if( fabs((*vertexes)[p1.sortedVtoV].x - (*vertexes)[p2.sortedVtoV].x) <= C3D_TOL ) {
			if( (*vertexes)[p1.sortedVtoV].y < (*vertexes)[p2.sortedVtoV].y ) {
				return true;
			}
		}

		return false;
	}

	vector<C2DVector> *vertexes;
};

// Dirty trick to find a vertex in the sorted array, given the vertex position.
// 'vIndex' must be an index which cannot be found in the sortVertex array
// (i.e. > sortVertex.size() or < 0); in this case, if p1 or p2 reference to
// this 'vIndex', the given vector is compared instead of the 'vertexes' element
//
// Remark1: see comments about LessThanV structure
//
// Remark 2: here C3D_TOL test is needed also on 'y' coordinate, because
//  we are searching for a point in the array, that is, coincident points!
//  And though the strict '<' should work, since we are searching exactly
//  for a point that IS in the array, the routine may still fail;
//  for istance, if in the input array we have two coincident points
//  within tolerance, only the reference to the first one is kept;
//  then, if we search for the *second* one, the 'y' coordinate
//  could be NOT equal!
struct LessThanVdirty  : public binary_function<C2DVertex &, C2DVertex &, bool>
{
	LessThanVdirty(vector<C2DVector> &varray, C2DVector vec, long vIndex)
	{
		vertexes = &varray;
		vVector = vec;
		vectIndex = vIndex;
	}

	// compare vertex positions according to 'pointed to' vertex
	bool operator()(const C2DVertex &p1, const C2DVertex &p2)
	{
		if(p1.sortedVtoV != vectIndex && p2.sortedVtoV != vectIndex) {
			if( (*vertexes)[p1.sortedVtoV].x + C3D_TOL < (*vertexes)[p2.sortedVtoV].x ) {
				return true;
			}
			else if( fabs((*vertexes)[p1.sortedVtoV].x - (*vertexes)[p2.sortedVtoV].x) <= C3D_TOL ) {
				if( (*vertexes)[p1.sortedVtoV].y  + C3D_TOL < (*vertexes)[p2.sortedVtoV].y ) {
					return true;
				}
			}
		}
		else if(p1.sortedVtoV == vectIndex && p2.sortedVtoV != vectIndex) {
			if( vVector.x + C3D_TOL < (*vertexes)[p2.sortedVtoV].x ) {
				return true;
			}
			else if( fabs(vVector.x - (*vertexes)[p2.sortedVtoV].x) <= C3D_TOL ) {
				if( vVector.y  + C3D_TOL < (*vertexes)[p2.sortedVtoV].y ) {
					return true;
				}
			}
		}
		else if(p1.sortedVtoV != vectIndex && p2.sortedVtoV == vectIndex) {
			if( (*vertexes)[p1.sortedVtoV].x + C3D_TOL < vVector.x ) {
				return true;
			}
			else if( fabs((*vertexes)[p1.sortedVtoV].x - vVector.x) <= C3D_TOL) {
				if( (*vertexes)[p1.sortedVtoV].y  + C3D_TOL < vVector.y ) {
					return true;
				}
			}
		}

		return false;
	}

	vector<C2DVector> *vertexes;
	C2DVector vVector;
	long vectIndex;
};

// The 'sortVertexes' array is used to sort input vertexes
// without modifying input 'vertexes' array. To do so,
// it contains member index references to the 'vertexes' array,
// called 'sortedVtoV'. In this way, we can find an element
// in 'vertexes' starting from 'sortVertexes'. However, we must
// be able also to do the opposite: from an element in 'vertexes',
// find an element in 'sortVertexes'. To do so, after the sort,
// the references vector 'vToSortedV' is initialized with
// the correct sequence of indexes.
void C2DTriangulate::SortVertexes()
{
	unsigned long i, vertexIndex;
	LessThanV sortComp(vertexes);

	sort(sortVertexes.begin(), sortVertexes.end(), sortComp);

	for(i=0; i < sortVertexes.size(); i++) {
		vertexIndex = sortVertexes[i].sortedVtoV;
		vToSortedV[vertexIndex] = i;
	}
}

// Initialize the 'sortVertexes' array, assiging values
// to member index references to the 'vertexes' array;
// moreover, remove duplicate vertexes (within tolerance)
void C2DTriangulate::InitSortVertexes()
{
	C2DBSPTree vertexTree(&vertexes);
	unsigned long i, j;
	long res;

	// assure that the 'vToSortedV' and 'sortVertexes' arrays
	// have exactly the same length of the 'vertexes' array
	// (if coincident points will be found, 'sortVertexes'
	// array will be cut short)
	vToSortedV.resize(vertexes.size());
	sortVertexes.resize(vertexes.size());

	for(i=0, j=0; i < vToSortedV.size(); i++) {
		// insert vertex in 2D BSP tree to test for uniqueness
		res = vertexTree.insert(i);
		// if not unique, mark it in the 'sortedVtoV' array
		// with a reference to the equivalent point
		// (reference is negated index number, -1 based;
		// it is -1 based because 0 cannot be negated)
		if( res != C2D_BSP_NEWPOINT) {
			vToSortedV[i] = -res-1;
		}
		// otherwise, insert reference to vertex in 'sortVertexes' array
		else {
			vToSortedV[i] = 0;
			sortVertexes[j].sortedVtoV = i;
			j++;
		}
	}

	// Following coincident points detection, 'sortVertexes'
	// may be shorter: cut it down
	sortVertexes.resize(j);

}

// find the 'sortVertexes' index corresponding to given C2DVector
long C2DTriangulate::FindVertex(C2DVector point)
{
	vector<C2DVertex>::iterator pos;
	LessThanVdirty sortComp(vertexes, point, C2D_DUMMY_INDEX);
	// dummy vertex (sortComp operator() will reference to
	// 'point' when finding index 'C2D_DUMMY_INDEX'
	C2DVertex dummyV(C2D_DUMMY_INDEX, NULL);

	// search for first occurence of point in sorted vector array;
	pos = lower_bound(sortVertexes.begin(), sortVertexes.end(), dummyV, sortComp) ;

    // if no matching element found
    if (pos == sortVertexes.end())
		return C2D_VERTEX_NOT_FOUND;

	return (pos - sortVertexes.begin());
}


// Delaunay triangulation of a set of points
// Remark: coincient points (also within tolerance!) are allowed
void C2DTriangulate::Triangulate()
{
	C2DEdge *right_cw, *left_ccw;

	// verify that the list contains at least 3 points
	ASSERT(vertexes.size() > 2);

	// clear the current mesh
	ClearMesh();

	// initialize the 'sortVertexes' array, which will contain
	// the sorted vertexes; also, duplicate vertexes are removed
	InitSortVertexes();

	// sort points first according to X, then to Y
	SortVertexes();

	// divide and conquer Delaunay triangulation
	DivideAndMerge(0, sortVertexes.size()-1, &right_cw, &left_ccw);

	// constrain the boundary edges of the convex hull
	// (needed in case a constrained triangulation will follow)
	ConstrainBoundary(left_ccw);

	// signal that mesh is present (memory has been allocated)
	meshIsPresent = true;
}

C2DTriangulate::~C2DTriangulate()
{
	ClearMesh();
}

void C2DTriangulate::ClearMesh()
{
	unsigned long i;
	C2DEdge *firstedge, *edge, *oldedge;

	if(meshIsPresent == true) {
		// scan all vertexes
		for(i=0; i<sortVertexes.size(); i++) {

			firstedge = edge = GetVertexEdge(i);

			// scan every edge sharing the same vertex, in CCW order
			do {
				// store reference to edge to be deleted
				oldedge = edge;
				// get next edge
				edge = edge->Onext();
				// delete old edge
				delete oldedge;

			} while(edge != firstedge);
		}

		meshIsPresent = false;
	}

	// empty vectors
	sortVertexes.clear();
	triangles.clear();
	smalltri.clear();
}

// Remark: all generated triangles are CCW
void C2DTriangulate::GenerateTriangles()
{
	unsigned long i, trinum;
	long mark;
	C2DEdge *firstedge, *edge, *edge2, *edge3;

	// if no mesh has been created, cannot generate triangles
	if(meshIsPresent == false)
		return;

	trinum = 0;
	mark = NewIndex();

	// reserve enough space for the triangulation
	// (max number of triangles is 2*N-2-K, where K is the
	// number of points on the convex hull)
	triangles.resize(sortVertexes.size()*2 - 2);

	// scan all vertexes
	for(i=0; i<sortVertexes.size(); i++) {

		firstedge = edge = GetVertexEdge(i);

		// scan every edge sharing the same vertex, in CCW order
		do {
			// if edge has not already been visited, check for triangles
			if(edge->GetId() != mark)
			{

				// if edge has a CCW triangle at its left,
				// output it
				if(HasEdgeLeftTri(edge, &edge2, &edge3)) {
					// store triangle
					triangles[trinum].SetEdges(edge, edge2, edge3);
					// mark edges as visited
					edge2->SetId(mark);
					edge3->SetId(mark);
					trinum++;
				}

				// mark edges as visited
				edge->SetId(mark);
			}

			// get next edge
			edge = edge->Onext();

		} while(edge != firstedge);
	}

	ASSERT(trinum <= triangles.size());

	// resize triangle list to actual length
	triangles.resize(trinum);
}

void C2DTriangulate::ConstrainBoundary(C2DEdge *boundaryedge)
{
	long i;
	C2DEdge *walkingedge;

	walkingedge = boundaryedge;

	// this should be an infinite loop; to avoid
	// a lock, a max loop count has been set
	i = 0;
	do {
		walkingedge->BoundaryConstrain();
		walkingedge = walkingedge->Mate()->Onext();
		i++;
	} while (walkingedge != boundaryedge && i < C2D_MAX_LOOP);

	// if max loop limit has been broken, assert
	if( i >= C2D_MAX_LOOP) {
		ASSERT(FALSE);
	}
}

void C2DTriangulate::DivideAndMerge(long left, long right, C2DEdge **right_cw, C2DEdge **left_ccw)
{
	long n, split;
	C2DEdge *edgea, *edgeb, *edgec;
	C2DEdge *right_cw_left, *left_ccw_left, *right_cw_right, *left_ccw_right;
	C2DEdge *lower_tangent;
	char joinside;

	// calculate the number of points left to triangulate,
	// extremes included
	n = right - left + 1;

	// bottom of the recursion, must build an edge
	if( n == 2 ) {
		// construct edge
		edgea = CreateEdge(left, right);
		// and set the initial edges used to find the
		// lowest tangent to the convex hull when merging
		// two triangulations; remember that points are
		// ordered: left always come before right
		*right_cw = edgea->Mate();
		*left_ccw = edgea;
	}
	// bottom of the recursion, must build a triangle
	else if( n == 3 ) {
		// construct first two edges
		edgea = CreateEdge(left, left+1);
		edgeb = CreateEdge(left+1, right);
		// and link them (on left+1 point);
		// note that edges can be oriented /\ or \/,
		// but always edgea->Mate will be onext to edgeb
		edgeb->LinkCCW(edgea->Mate());

		// To close the triangle with the third
		// edge, order matters: /\ or \/ is, in general, different.
		// In this case, however, it is the same, since there is
		// no other triangle present (so Onext() and Oprev() of edgea
		// is always edgea, the same for edgeb); if other neighbouring
		// triangles are present, this is not true.
		edgec = JoinEdges(edgea, edgeb, &joinside);

		if( joinside == C2D_LEFT_OF ) {
			*left_ccw = edgec;
			*right_cw = edgec->Mate();
		}
		else if( joinside == C2D_RIGHT_OF ) {
			*left_ccw = edgea;
			*right_cw = edgeb->Mate();
		}
		else if( joinside == C2D_ON_SEGMENT ) {
			*left_ccw = edgea;
			*right_cw = edgeb->Mate();
		}
		else {
			// this should never happen
			ASSERT(FALSE);
		}
	}
	// more than 3 points, divide again, triangulating both
	// left and right halves, then merge the two triangulations
	else {
		// compute the split point
		split = (left + right) / 2;

		// tringulate left and right halves
		DivideAndMerge(left, split, &right_cw_left, &left_ccw_left);
		DivideAndMerge(split+1, right, &right_cw_right, &left_ccw_right);

		// merge
		lower_tangent = Merge(right_cw_left, left_ccw_right);

		// the lower tangent added by merge operation could
		// have substituted either left_ccw_left or right_cw_right;
		// check them and update, if necessary

		// if lower_tangent touches the old left half convex hull
		// at the leftmost point, then the boundary edge is now lower_tangent;
		// note that the old left_ccw_left could have been completly
		// invalidated, since not only it's not on the boundary any more,
		// but could have been deleted during the merge step.
		if( lower_tangent->Vertex() == left )
			left_ccw_left = lower_tangent;

		// if lower_tangent->Mate() touches the old right half convex hull
		// at the rightmost point, then the boundary edge is now lower_tangent;
		// note that the old right_cw_right could have been completly
		// invalidated, since not only it's not on the boundary any more,
		// but could have been deleted during the merge step.
		if( lower_tangent->MateVertex() == right )
			right_cw_right = lower_tangent->Mate();

		// update references to outermost boundary edges
		// to be passed back
		*left_ccw = left_ccw_left;
		*right_cw = right_cw_right;
	}
}

void C2DTriangulate::FindLowestTangent(C2DEdge *right_cw_left, C2DEdge *left_ccw_right,
									   C2DEdge **lowest_left, C2DEdge **lowest_right)
{
	C2DEdge *r_cw_l, *l_ccw_r;
	C2DVector *origin_left, *dest_left, *origin_right, *dest_right;
	long i;

	// copy starting segments
	r_cw_l = right_cw_left;
	l_ccw_r = left_ccw_right;

	// this should be an infinite loop; to avoid
	// a lock, a max loop count has been set
	for( i = 0; i < C2D_MAX_LOOP; i++) {
		// get all involved points
		origin_left = Vertex(r_cw_l);
		dest_left = MateVertex(r_cw_l);
		origin_right = Vertex(l_ccw_r);
		dest_right = MateVertex(l_ccw_r);

		// if the origin_left is still left of the (supposed) lowest tangent,
		// the segment cannot be the lowest tangent; so move on around
		// left contour.
		// Note that if origin_left is on the lowest tangent, it's ok to
		// stop here, since this is the nearest point (we are circling cw)
		if(origin_left->IsLeftOrRight(dest_left, origin_right) == C2D_LEFT_OF) {
			r_cw_l = r_cw_l->Mate()->Oprev();
		}
		// if the origin_right is still right of the (supposed) lowest tangent,
		// the segment cannot be the lowest tangent; so move on around
		// right contour.
		// Note that if origin_right is on the lowest tangent, it's ok to
		// stop here, since this is the nearest point (we are circling ccw)
		else if(origin_right->IsLeftOrRight(dest_right, origin_left) == C2D_RIGHT_OF) {
			l_ccw_r = l_ccw_r->Mate()->Onext();
		}
		else
			// lowest tangent found
			break;
	}

	// if max loop limit has been broken, assert
	if( i >= C2D_MAX_LOOP) {
		ASSERT(FALSE);
	}

	*lowest_left = r_cw_l;
	*lowest_right = l_ccw_r;
}

// merge two Delaunay triangulations;
// returns the two group's lowest tangent, half edge oriented from left to right
C2DEdge *C2DTriangulate::Merge(C2DEdge *right_cw_left, C2DEdge *left_ccw_right)
{
	C2DEdge *lowest_left, *lowest_right, *lowest_tangent;
	C2DEdge *base_edge, *left_candidate, *right_candidate, *next_left_cand, *next_right_cand;
	long i, j;
	char incircle_left;

	// find lowest tangent of the two Delaunay triangulations
	FindLowestTangent(right_cw_left, left_ccw_right, &lowest_left, &lowest_right);

	// create new edge (lowest tangent)
	lowest_tangent = JoinEdgesRight(lowest_left, lowest_right);

	// initailize base edge
	base_edge = lowest_tangent;

	// merge the two triangulations

	// this should be an infinite loop; to avoid
	// a lock, a max loop count has been set
	for( j = 0; j < C2D_MAX_LOOP; j++) {

		// initialize left and right candidates
		left_candidate = base_edge->Onext();
		right_candidate = base_edge->Mate()->Oprev();

		//
		// find left candidate
		//

		// this should be an infinite loop; to avoid
		// a lock, a max loop count has been set
		for( i = 0; i < C2D_MAX_LOOP; i++) {

			// verify that the angle between the base_edge and the left_candidate is less than 180 degrees
			if( MateVertex(left_candidate)->IsLeftOrRight(Vertex(base_edge), MateVertex(base_edge)) != C2D_LEFT_OF) {
				// if not, no left candidate is submitted
				left_candidate = NULL;
				break;
			}

			// verify that the circle passing through the base_edge and left_candidate vertexes
			// does not contain the next left candidate end point

			next_left_cand = left_candidate->Onext();

			// verify that the angle between the base_edge and the next_left_cand is less than 180 degrees
			if( MateVertex(next_left_cand)->IsLeftOrRight(Vertex(base_edge), MateVertex(base_edge)) != C2D_LEFT_OF) {
				// if not, there is no next left candidate, so left candidate is ok, can submit it
				break;
			}

			// in-circle test
			if( MateVertex(next_left_cand)->InCircle(MateVertex(left_candidate), Vertex(base_edge), MateVertex(base_edge)) != C2D_IN_CIRCLE) {
				// if not, left candidate is ok, can submit it
				break;
			}

			// otherwise, must delete old left_candidate; next_left_cand will become the
			// new left_candidate

			DeleteEdge(left_candidate);

			left_candidate = next_left_cand;

		}

		// if max loop limit has been broken, assert
		if( i >= C2D_MAX_LOOP) {
			ASSERT(FALSE);
		}

		//
		// find right candidate
		//

		// this should be an infinite loop; to avoid
		// a lock, a max loop count has been set
		for( i = 0; i < C2D_MAX_LOOP; i++) {

			// verify that the angle between the base_edge and the right_candidate is less than 180 degrees
			if( MateVertex(right_candidate)->IsLeftOrRight(MateVertex(base_edge), Vertex(base_edge)) != C2D_RIGHT_OF) {
				// if not, no right candidate is submitted
				right_candidate = NULL;
				break;
			}

			// verify that the circle passing through the base_edge and right_candidate vertexes
			// does not contain the next right candidate end point

			next_right_cand = right_candidate->Oprev();

			// verify that the angle between the base_edge and the next_right_cand is less than 180 degrees
			if( MateVertex(next_right_cand)->IsLeftOrRight(MateVertex(base_edge), Vertex(base_edge)) != C2D_RIGHT_OF) {
				// if not, there is no next right candidate, so right candidate is ok, can submit it
				break;
			}

			// in-circle test
			if( MateVertex(next_right_cand)->InCircle(MateVertex(right_candidate), Vertex(base_edge), MateVertex(base_edge)) != C2D_IN_CIRCLE) {
				// if not, right candidate is ok, can submit it
				break;
			}

			// otherwise, must delete old right_candidate; next_right_cand will become the
			// new right_candidate

			DeleteEdge(right_candidate);

			right_candidate = next_right_cand;

		}

		// if max loop limit has been broken, assert
		if( i >= C2D_MAX_LOOP) {
			ASSERT(FALSE);
		}

		// now must choose between the two candidates

		// if no candidate, merge is completed
		if( left_candidate == NULL && right_candidate == NULL ) {
			break;
		}

		// if no left_candidate, connect right candidate
		if( left_candidate == NULL ) {
			base_edge = JoinEdgesRight(base_edge, right_candidate->Mate());
		}
		// if no right_candidate, connect left candidate
		else if( right_candidate == NULL ) {
			base_edge = JoinEdgesRight(left_candidate->Mate(), base_edge->Mate());
		}
		// if both candidates are present, choose the one which does
		// not 'InCircle' the other
		else {
			// should test both candidates, one against the other,  for incircle property;
			// however, for the existence and uniqueness of the Delaunay triangulation,
			// one and only one of the two tests must fail. So one suffices.
			incircle_left = MateVertex(left_candidate)->InCircle(MateVertex(right_candidate), Vertex(base_edge), MateVertex(base_edge));

			// if right candidate contains left candidate, cannot be Delaunay;
			// so the correct candidate is the left one
			if(incircle_left == C2D_IN_CIRCLE) {
				// connect left candidate
				base_edge = JoinEdgesRight(left_candidate->Mate(), base_edge->Mate());
			}
			else {
				// connect right candidate
				base_edge = JoinEdgesRight(base_edge, right_candidate->Mate());
			}

		}
	}

	return lowest_tangent;
}

// insert a constrained edge between points p1 and p2
// (of course, they must already be in the mesh)
bool C2DTriangulate::InsertConstrEdge(C2DVector p1, C2DVector p2)
{
	long a, b;
	bool res;

	// Must be sure that vertexes are sorted, since crossing
	// constrained edges may have indroduced new points
	// on the array tail; these non-ordered points may
	// cause the array template binary_search() function in FindVertex()
	// to fail (depending on how it is implemented)
	SortVertexes();

	// find index corresponding to given vertex
	a = FindVertex(p1);
	b = FindVertex(p2);

	// if either vertex is not in the mesh, return
	if( a == C2D_VERTEX_NOT_FOUND || b == C2D_VERTEX_NOT_FOUND )
		return false;

	// convert index from 'sortVertexes' to 'vertexes'
	// (InsertConstrEdge() refers to indexes of 'vertexes'
	// array, which is the user input)
	a = GetVIndexFromSortV(a);
	b = GetVIndexFromSortV(b);

	res = InsertConstrEdge(a,b);
	// a,b must always exist!
	ASSERT(res == true);

	return true;
}

// insert a constrained edge between points with 'vertexes' index a and b
// (of course, they must already be in the mesh)
// Remark: vertexes array is shuffled by sort operation
// when performing divide-and-conquer Delaunay triangulation.
// So vertex order in vertexes array is not any more the
// initial one; it follows that a and b must refer to
// new position of the points in the array
bool C2DTriangulate::InsertConstrEdge(long a, long b)
{
	C2DEdge *edgea, *edgeb, *walkingedge, *walkinglnext, *newedge;
	char halfplane, res;
	bool isAligned;
	long i, j, size;
	C2DLine lineab;
	C2DVector point;

	// test indexes a and b for validity
	size = vertexes.size();
	if( a < 0 || a >= size || b < 0 || b >= size )
		return false;

	// convert index from 'vertexes' to 'sortVertexes'
	// (InsertConstrEdge() refers to indexes of 'vertexes'
	// array, which is the user input)
	a = GetSortVIndexFromV(a);
	b = GetSortVIndexFromV(b);

	// since coincident input points are allowed (and mapped to the same indexes),
	// a and b could be coincident. In this case, no constrained edge
	// is needed (only one point, already in the mesh)
	if( a == b )
		return true;

	// must connect a with b using a sequence of one or more edges;
	// (more than one edge could be used when some segments or some
	// vertexes are already on the path of the edge)

	// line through a,b
	lineab = C2DLine(*Vertex(a), *Vertex(b));

	// find two edges with origin a and b respectively
	edgea = GetVertexEdge(a);
	edgeb = GetVertexEdge(b);

	// loop until the origin of edgea is not coincident with the last
	// point, b, of the constrained edge
	// (note: this loop could be avoided using recursion, splitting
	// the operation of connecting a to b in two subtasks consisting of
	// a connection operation of a to the first vertex already on (a,b)
	// and then calling again InsertConstrEdge() to connect this
	// vertex to b. However, infinite loop cases are more easily
	// detected in a loop than in recursion)

	// this should be an infinite loop; to avoid
	// a lock, a max loop count has been set
	for( i = 0; i < C2D_MAX_LOOP && edgea->Vertex() != b; i++) {

		// Search the the first edge to the right of edge (a,b) or
		// aligned with (a,b) and store it in edgea
		// To perform the search, move ccw from current edgea

		isAligned = false;

		// this should be an infinite loop; to avoid
		// a lock, a max loop count has been set
		for( j = 0; j < C2D_MAX_LOOP; j++) {

			// if the destination of edgea is on the (a,b) line
			if( lineab.IsPointOnLine(*MateVertex(edgea)) ) {

				// and in the positive direction from a (that is,
				// not before a)

				halfplane = lineab.HalfPlane(*MateVertex(edgea));
				// if coincident with line origin, this is an internal error
				ASSERT(halfplane != C2D_ON_ORIGIN);
				if( halfplane == C2D_POS_HALFPLANE) {

					// if so, we have already found a piece of the
					// constrained edge
					isAligned = true;
					break;
				}
			}

			// if the destination of edgea->Onext() is on the (a,b) line
			if( lineab.IsPointOnLine(*MateVertex(edgea->Onext())) ) {

				// and in the positive direction from a (that is,
				// not before a)

				halfplane = lineab.HalfPlane(*MateVertex(edgea->Onext()));
				// if coincident with line origin, this is an internal error
				ASSERT(halfplane != C2D_ON_ORIGIN);
				if( halfplane == C2D_POS_HALFPLANE) {

					// if so, we have already found a piece of the
					// constrained edge
					edgea = edgea->Onext();
					isAligned = true;
					break;
				}
			}

			// if edgea is on the right of line (a,b)
			// Remark: no need to test within tolerance: the edges
			// have already been tested NOT to be aligned with (a,b)
			// within tolerance
			if( ::CrossProd(C2DVector(*Vertex(a), *MateVertex(edgea)), lineab) > 0 ) {
				// and edgea->Onext() is on the left of line (a,b)
				if( ::CrossProd(C2DVector(*Vertex(a), *MateVertex(edgea->Onext())), lineab) < 0 ) {
					// edgea is the first edge to the right of (a.b)
					break;
				}
			}

			// otherwise continue the search
			edgea = edgea->Onext();

		}

		// if max loop limit has been broken, assert
		if( j >= C2D_MAX_LOOP) {
			ASSERT(FALSE);
		}

		// if edgea is aligned with segment (a,b), simply constrain
		// edgea and go on inserting a new constrained edge
		// staring from edgea->MateVertex()
		if( isAligned == true ) {
			edgea->Constrain();
			// new starting point 'a'
			a = edgea->MateVertex();
			// and new line (a,b), but only if new a is not the end of the line
			if( a != b )
				lineab = C2DLine(*Vertex(a), *Vertex(b));
			// new edgea, as near as possible to the rightmost edge w.t.r (a.b)
			edgea = edgea->Mate()->Onext();
			continue;
		}

		// initialize edge which will 'walk' along the simplex
		// contour around the constrained edge to be built;
		// in this way, all edges crossing the path of the
		// constrained edge will be identified and deleted
		walkingedge = edgea;

		// this should be an infinite loop; to avoid
		// a lock, a max loop count has been set
		for( j = 0; j < C2D_MAX_LOOP; j++) {
			// get next edge on the left, since it will be tested intensively
			walkinglnext = walkingedge->Lnext();
			// if there is a vertex along the path of the edge (a,b),
			// must insert a segment of the total constrained edge
			if( lineab.IsPointOnLine(*MateVertex(walkinglnext)) ) {
				// create new edge
				newedge = JoinEdgesRight(edgea, walkinglnext->Mate());
				newedge->Constrain();
				// new edgea
				edgea = walkinglnext->Mate();
				// and re-triangulate the regions on the right and on the left
				// of the new edge, which could be empty (the edges crossing
				// the new constrained one have been deleted)
				Retriangulate(newedge->Oprev());
				Retriangulate(newedge->Lnext());
				// start again from beginning
				break;
			}
			// if that edge is on the left of line (a,b), this means that
			// we have crossed from right to left the path of the constrained edge
			// Remark: no need to test within tolerance: the edge
			// have already been tested NOT to be aligned with (a,b)
			// within tolerance
			if( ::CrossProd(C2DVector(*Vertex(a), *MateVertex(walkinglnext)), lineab) < 0 ) {
				// therefore, if the crossing edge is not constrained
				if(walkinglnext->IsConstrained() == false) {
					DeleteEdge(walkinglnext);
				}
				// otherwise, must split crossing constrained edge
				else {
					// intersect edges
					res = IntersSegmentLine(point, *Vertex(walkinglnext), *MateVertex(walkinglnext), lineab);
					ASSERT(res != C2D_SEGMENT_DONT_INTERSECT);

					// split crossing constrained edge, inserting a new vertex
					// at point position
					Split(walkinglnext, point);

					// now can continue iterating; at next iteration, a vertex will
					// be found along the path of the constrained edge being inserted
				}
			}
			// if no vertex on or no edge crossing the path
			// of the constrained edge being inserted,
			// go on
			else {
				walkingedge = walkinglnext;
			}
		}

		// if max loop limit has been broken, assert
		if( j >= C2D_MAX_LOOP) {
			ASSERT(FALSE);
		}

	}

	// if max loop limit has been broken, assert
	if( i >= C2D_MAX_LOOP) {
		ASSERT(FALSE);
    }

	return true;
}

// Retriangulate a region inside a Delaunay triangulation
// made free because of edge deletion caused by constrained
// edge insertion.
// This is not a general-purpose triangulation routine, since
// the polygon to be triangulated, though not monotone, is
// close enough to have some nice properties.
// The algorithm is based on some sort of edge clipping;
// it triangulates the left face of 'first', which is assumed
// to be closed. Moreover, called 'last' the edge preceding
// 'first' ccw, it is assumed that all vertices of the face lie to
// the left of 'last'.
void C2DTriangulate::Retriangulate(C2DEdge *first)
{
	long i, j;
	C2DEdge *last, *walkingedge, *walkprev, *walknext;

	// initialize 'last' edge
	last = first->Lprev();

	// continue to loop until the left face becomes a triangle
	//
	// this should be an infinite loop; to avoid
	// a lock, a max loop count has been set
	for( i = 0; i < C2D_MAX_LOOP && first->Lnext()->Lnext() != last; i++) {

		// initialize walking edges
		// (these edges will 'walk' along the polygon boundary
		// searching for ears)
		walkingedge = first->Lnext();
		walkprev = first;

		// continue to loop until 'walkingedge' has walked the whole boundary
		// of the polygon
		//
		// this should be an infinite loop; to avoid
		// a lock, a max loop count has been set
		for( j = 0; j < C2D_MAX_LOOP && walkingedge != last; j++) {

			walknext = walkingedge->Lnext();

			// if only one triangle left (walkprev is always the edge
			// Lprev() to 'edge', walknext the edge Lnext() to 'edge')
			if(walkprev == first && walknext == last)
				break;

			// if this is an ear (that is, these two edges are
			// convex), can cut it out
			if(CrossProd(walkprev, walkingedge) > C3D_TOL) {
				// If an ear, must cut it. However, two cases
				// may happen: the 'walkprev' is also the 'first'
				// or not. In the former case, the new edge
				// will replace 'first'; in any case, will
				// also be the new 'walkprev'
				if( walkprev == first ) {
					walkprev = JoinEdgesRight(walkprev, walkingedge->Mate());
					first = walkprev;
				}
				else
					walkprev = JoinEdgesRight(walkprev, walkingedge->Mate());

				// move on with 'walkingedge'
				walkingedge = walknext;

				// perform any swap required on the two other edges
				// of the newly formed triangle for a Delaunay constrained
				// triangulation
				DelaunayEdge(walkprev->Oprev());
				DelaunayEdge(walkprev->Mate()->Onext()->Mate());
			}
			// otherwise, if not an ear
			else {
				// move on with 'walkingedge'
				walkprev = walkingedge;
				walkingedge = walknext;
			}
		}

		// if max loop limit has been broken, assert
		if( j >= C2D_MAX_LOOP) {
			ASSERT(FALSE);
		}
	}

	// if max loop limit has been broken, assert
	if( i >= C2D_MAX_LOOP) {
		ASSERT(FALSE);
	}
}

// Swap given edge, if required, for mantaining the Delaunay property;
// also recursively call the procedure on the side edges in case of
// a swap.
// Remark: the recursive call is made only on the segments on the
// right of 'edge', see Retriangulate(). Note that this is enough
// to mantain the Delaunay property (in case of a swap, the segments
// on the old left side will be tested in the next recursive calls.
void C2DTriangulate::DelaunayEdge(C2DEdge *edge)
{
	C2DEdge *oprev, *dnext;

	// if this is a constrained edge, cannot do anything
	if(edge->IsConstrained() == true)
		return;

	oprev = edge->Oprev();
	dnext = edge->Mate()->Oprev();

	// if Delaunay property is violated
	if( MateVertex(oprev)->InCircle(MateVertex(edge), MateVertex(edge->Onext()), Vertex(edge)) == C2D_IN_CIRCLE) {
		// swap edges
		SwapEdge(edge);
		// and fix the neigborhood
		DelaunayEdge(oprev);
		DelaunayEdge(dnext);
	}
}

// refine computed CDT; epsilon is minimum angle allowed, in degrees
void C2DTriangulate::Refine(double epsilon)
{
	C2DEdgeTriangle tmptri;
	C2DEdge *edges[3], *shortestEdge, *mediumEdge;
	C2DVector splitpoint;
	double len[3], cosAngle[3], t0;
	unsigned char lenOrder[3], angleOrder[3];
	bool boundaryEdges[3], isTriangChanged;
	unsigned long i;

	// if no mesh has been created, cannot refine it
	if(meshIsPresent == false)
		return;

	// save minimum angle value
	cosEpsilon = cos(epsilon * CMATH_PI / 180.0);

	//
	// find the list of the worst-angled triangles
	//

	// clear smalltri list
	smalltri.clear();

	// generate triangles
	GenerateTriangles();


	// scan all triangles, pre-processing the boundary triangles which could
	// lead to infinite recursion, in case of condition 3.
	isTriangChanged = false;
	for(i=0; i<triangles.size(); i++) {

		// get edge lengths
		GetLenOrder(triangles[i], edges, len, lenOrder);

		// find out which edges are boundary edges
		FindBoundaryEdges(edges, boundaryEdges);

		// if shortest and medium edges are boundary edges,
		// condition 3 may hold
		if( boundaryEdges[lenOrder[0]] == true && boundaryEdges[lenOrder[1]] == true) {

			// get angles
			ComputeCosAngles(triangles[i], cosAngle);
			GetDoubleOrder(cosAngle, angleOrder);

			// if smallest angle is < 30 degrees
			// and largest angle is <= 120 degrees
			// cycle condition 3 may happen
			// (note that smallest angle has largest cosine)
			if( cosAngle[angleOrder[2]] > COS30DEG && cosAngle[angleOrder[0]] >= COS120DEG) {
				// if half of medium edge is longer than shortest edge
				// condition 3 may hold
				if(len[lenOrder[1]] > len[lenOrder[0]] * 2) {

					// insert in medium edge a new point Q such that
					// the resulting two boundary edges triangle
					// is isosceles

					// relative distance of point Q
					t0 = len[lenOrder[0]] / len[lenOrder[1]];

					// now must find the order in which the shortest and
					// medium edges appear, to split medium edge
					// on the correct side w.r.t. the shortest edge
					shortestEdge = edges[lenOrder[0]];
					mediumEdge = edges[lenOrder[1]];
					if( shortestEdge->MateVertex() == mediumEdge->Vertex() )
						splitpoint = PointOnSegment(*Vertex(mediumEdge), *MateVertex(mediumEdge), t0);
					else
						splitpoint = PointOnSegment(*MateVertex(mediumEdge), *Vertex(mediumEdge), t0);

					// split boundary edge, inserting a new vertex at midpoint
					Split(mediumEdge, splitpoint);

					// insert a new edge, thus creating two new triangles
					JoinEdgesRight(mediumEdge->Onext()->Mate(), mediumEdge->Mate());

					// signal that triangulation is changed
					isTriangChanged = true;
				}
			}
		}
	}

	// generate triangles again, if previous procedure has split
	// any triangle
	if( isTriangChanged == true)
		GenerateTriangles();

	// scan all triangles, searching for any triangle whose smallest edge is < epsilon
	// and inserting it in smalltri list
	for(i=0; i<triangles.size(); i++) {

		SmallTriangleInsert(triangles[i]);
	}


	//
	// angle refinement process
	//

	// refine small-angled triangles until no more triangles
	// are present in the 'smalltri' list

	// this should be an infinite loop; to avoid
	// a lock, a max loop count has been set
	for( i = 0; i < C2D_MAX_LOOP && smalltri.size() > 0; i++) {
		// get first triangle in list
		tmptri = *(smalltri.begin());
		// and refine it
		RefineTriangle(tmptri);
	// small-angled triangle list could have changed during the process;
	// so just test that set is not empty
	}

	// if max loop limit has been broken, assert
	if( i >= C2D_MAX_LOOP) {
		ASSERT(FALSE);
	}
}

// insert given triangle in the small triangles list, if small-angled
void C2DTriangulate::SmallTriangleInsert(C2DEdgeTriangle &triangle)
{
	double cosAngle[3];
	C2DEdge *edge1, *edge2, *edge3;

	// get triangle edges
	edge1 = triangle.GetEdge(0);
	edge2 = triangle.GetEdge(1);
	edge3 = triangle.GetEdge(2);

	ComputeCosAngles(triangle, cosAngle);

	if( (cosAngle[0] > cosEpsilon && (edge1->IsConstrained() == false || edge2->IsConstrained() == false) ) ||
		(cosAngle[1] > cosEpsilon && (edge2->IsConstrained() == false || edge3->IsConstrained() == false) ) ||
		(cosAngle[2] > cosEpsilon && (edge3->IsConstrained() == false || edge1->IsConstrained() == false) ) ) {

		AddSmallTriangle(triangle);
	}
}

// refine given small-angle triangle
void C2DTriangulate::RefineTriangle(C2DEdgeTriangle &triangle)
{
	long i;

	// while 'triangle' remains an actual triangle in the mesh
	// (i.e. it has not been refined), continue to iterate

	// this should be an infinite loop; to avoid
	// a lock, a max loop count has been set
	for( i = 0; i < C2D_MAX_LOOP && IsMeshTriangle(triangle) == true; i++) {

		LEPPPointInsert(triangle);
	}

	// if max loop limit has been broken, assert
	if( i >= C2D_MAX_LOOP) {
		ASSERT(FALSE);
	}
}

// refine the triangulation following the LEPP of the given
// triangle and inserting a point on the LEPP terminal or
// boundary edge
void C2DTriangulate::LEPPPointInsert(C2DEdgeTriangle triangle)
{
	C2DEdgeTriangle leppTriangle;
	C2DEdge *edges[3], *longestEdge, *longestEdgeMate, *leppEdges[3];
	C2DEdge *mediumEdge, *mediumEdgeMate;
	double len[3], leppLen[3], cosAngle[3];
	unsigned char lenOrder[3], leppLenOrder[3], angleOrder[3];
	bool boundaryEdges[3];
	long i;

	// this should be an infinite loop; to avoid
	// a lock, a max loop count has been set
	for( i = 0; i < C2D_MAX_LOOP; i++) {

		// get edge lengths
		GetLenOrder(triangle, edges, len, lenOrder);

		//
		// check conditions for boundary edge point insertion
		//

		// if current triangle is a boundary triangle with one boundary
		// edge not equal to its smallest edge, can perform boundary
		// point insertion

		FindBoundaryEdges(edges, boundaryEdges);

		// if longest edge is a boundary edge, perform boundary
		// point insertion
		if( boundaryEdges[lenOrder[2]] == true ) {
			BoundaryPointInsert(triangle, lenOrder[2]);
			return;
		}
		// otherwise, if medium-size edge is a boundary edge, perform boundary
		// point insertion; this condition is imposed to prevent fractal
		// behaviour, due to condition 1. In this case, also condition
		// about largest angle <= 120 degrees must be checked; this relates
		// to condition 3 also and is used to avoid over-refinement
		if( boundaryEdges[lenOrder[1]] == true ) {

			ComputeCosAngles(triangle, cosAngle);
			GetDoubleOrder(cosAngle, angleOrder);

			// if largest angle is <= 120 degrees, can perform
			// boundary point insertion
			if( cosAngle[angleOrder[0]] >= COS120DEG ) {

				BoundaryPointInsert(triangle, lenOrder[1]);
				return;
			}
		}

		//
		// otherwise must check if this is a terminal triangle or not
		//

		// get next triangle in LEPP (must exist, because boundary
		// triangle conditions have already been checked)
		longestEdge = edges[lenOrder[2]];
		longestEdgeMate = longestEdge->Mate();
		leppTriangle.SetEdges(longestEdgeMate, longestEdgeMate->Lnext(), longestEdgeMate->Lprev());

		// order edges according to length
		GetLenOrder(leppTriangle, leppEdges, leppLen, leppLenOrder);

		// if this is a terminal triangle, terminal edge insert
		//
		// if the two triangles share the longest edge, this is
		// a terminal triangle
		if( longestEdge == leppEdges[leppLenOrder[2]]->Mate() )
		{
			// check the non-cycle condition; if infinite recursion could
			// happen, avoid it by inserting a point on the medium sized
			// edge of the terminal triangle instead of on the longest one
			if(NonCycleTriangle(triangle, edges, lenOrder) == false) {
				// medium size edge of 'triangle' is NOT a boundary edge;
				// this case is catched before by boundary point insertion conditions,
				// moreover the cycle condition cannot happen.
				// So another triangle sharing the medium edge must exist;
				// find it and perform terminal edge insertion on this edge.
				mediumEdge = edges[lenOrder[1]];
				mediumEdgeMate = mediumEdge->Mate();
				leppTriangle.SetEdges(mediumEdgeMate, mediumEdgeMate->Lnext(), mediumEdgeMate->Lprev());

				TerminalPointInsert(triangle, lenOrder[1], leppTriangle);
				return;
			}

			if(NonCycleTriangle(leppTriangle, leppEdges, leppLenOrder) == false) {
				// medium size edge of 'leppTriangle' is NOT a boundary edge;
				// otherwise, the cycle condition cannot happen
				// So another triangle sharing the medium edge must exist;
				// find it and perform terminal edge insertion on this edge.
				mediumEdge = leppEdges[leppLenOrder[1]];
				mediumEdgeMate = mediumEdge->Mate();
				triangle.SetEdges(mediumEdgeMate, mediumEdgeMate->Lnext(), mediumEdgeMate->Lprev());

				TerminalPointInsert(leppTriangle, leppLenOrder[1], triangle);
				return;
			}

			// otherwise no fractal behaviour can happen, so insert terminal point
			TerminalPointInsert(triangle, lenOrder[2], leppTriangle);
			return;
		}

		// if this is not a terminal triangle,
		// continue to follow the triangle's LEPP
		triangle = leppTriangle;
	}

	// if max loop limit has been broken, assert
	if( i >= C2D_MAX_LOOP) {
		ASSERT(FALSE);
	}
}

// split a boundary triangle inserting a point on the boundary edge
// and 'delaunay' edges
void C2DTriangulate::BoundaryPointInsert(C2DEdgeTriangle &triangle, char edgeindex)
{
	C2DVector midpoint;
	C2DEdge *edge, *newedge, *newedgemate;

	// remove triangle from smalltri list, if present
	RemoveSmallTriangle(triangle);

	edge = triangle.GetEdge(edgeindex);

	// find edge midpoint
	midpoint = (*MateVertex(edge) + *Vertex(edge)) / 2.0;

	// split boundary edge, inserting a new vertex at midpoint
	Split(edge, midpoint);

	// insert a new edge, thus creating two new triangles
	newedge = JoinEdgesRight(edge->Onext()->Mate(), edge->Mate());

	newedgemate = newedge->Mate();

	// if small-angled, insert in small triangles set
	SmallTriangleInsert(newedge);
	SmallTriangleInsert(newedgemate);

	// and Delaunay edge-swap the four sides
	//
	// right side
	LEPPDelaunayEdge(newedge->Oprev());
	LEPPDelaunayEdge(newedge->Mate()->Onext()->Mate());
	// left side
	LEPPDelaunayEdge(newedge->Onext());
	LEPPDelaunayEdge(newedge->Mate()->Oprev()->Mate());
}

// split a boundary triangle inserting a point on the boundary edge
void C2DTriangulate::TerminalPointInsert(C2DEdgeTriangle &triangle, char edgeindex, C2DEdgeTriangle &triangle2)
{
	C2DVector midpoint;
	C2DEdge *edge, *newedge1, *newedge2, *newedgemate1, *newedgemate2;

	// remove triangles from smalltri list, if present
	RemoveSmallTriangle(triangle);
	RemoveSmallTriangle(triangle2);

	edge = triangle.GetEdge(edgeindex);

	// find edge midpoint
	midpoint = (*MateVertex(edge) + *Vertex(edge)) / 2.0;

	// split boundary edge, inserting a new vertex at midpoint
	edge = Split(edge, midpoint);

	// insert a new edge in the first triangle, thus creating two new triangles
	newedge1 = JoinEdgesRight(edge->Onext()->Mate(), edge->Mate());
	// insert a new edge in the second triangle, thus creating other two new triangles
	newedge2 = JoinEdgesLeft(edge->Oprev()->Mate(), edge->Mate());

	newedgemate1 = newedge1->Mate();
	newedgemate2 = newedge2->Mate();

	// if small-angled, insert in small triangles set
	SmallTriangleInsert(newedge1);
	SmallTriangleInsert(newedgemate1);
	SmallTriangleInsert(newedge2);
	SmallTriangleInsert(newedgemate2);

	// and Delaunay edge-swap the six sides
	//
	// newedge1 right side
	LEPPDelaunayEdge(newedge1->Oprev());
	LEPPDelaunayEdge(newedge1->Mate()->Onext()->Mate());
	// newedge1 left side
	LEPPDelaunayEdge(newedge1->Onext());
	LEPPDelaunayEdge(newedge1->Mate()->Oprev()->Mate());
	// newedge2 left side (half left side is in common with triangle 1)
	LEPPDelaunayEdge(newedge2->Onext());
	// newedge2 right side (half left side is in common with triangle 1)
	LEPPDelaunayEdge(newedge1->Oprev());
}

// Swap given edge, if required, for mantaining the Delaunay property;
// also recursively call the procedure on the side edges in case of
// a swap.
// Remark: the recursive call is made only on the segments on the
// right of 'edge', see Retriangulate(). Note that this is enough
// to mantain the Delaunay property (in case of a swap, the segments
// on the old left side will be tested in the next recursive calls.
//
// This routine differs from DelaunayEdge() in that adds/removes
// small-angled triangles from small-angled triangles list
void C2DTriangulate::LEPPDelaunayEdge(C2DEdge *edge)
{
	C2DEdge *oprev, *dnext, *edgemate;
	C2DEdgeTriangle tri;

	// if this is a constrained edge, cannot do anything
	if(edge->IsConstrained() == true)
		return;

	oprev = edge->Oprev();
	dnext = edge->Mate()->Oprev();

	// if Delaunay property is violated
	if( MateVertex(oprev)->InCircle(MateVertex(edge), MateVertex(edge->Onext()), Vertex(edge)) == C2D_IN_CIRCLE) {

		// remove triangles from smalltri list, if present
		edgemate = edge->Mate();
		RemoveSmallTriangle(edge);
		RemoveSmallTriangle(edgemate);

		// swap edges
		edge = SwapEdge(edge);

		edgemate = edge->Mate();

		// if small-angled, insert in small triangles set
		SmallTriangleInsert(edge);
		SmallTriangleInsert(edgemate);

		// and fix the neigborhood
		LEPPDelaunayEdge(oprev);
		LEPPDelaunayEdge(dnext);
	}
}

// Delete all non-strong constrained edges
// CW of given edge (on its right)
// This routine is used to perform holes in the formed mesh
// (usually before refining it)
bool C2DTriangulate::DeleteCWEdges(C2DVector p1, C2DVector p2)
{
	long a, b;
	bool res;

	// Must be sure that vertexes are sorted, since crossing
	// constrained edges may have indroduced new points
	// on the array tail; these non-ordered points may
	// cause the array template binary_search() function in FindVertex()
	// to fail (depending on how it is implemented).
	// This works also in case Refine() has been called before DeleteCWEdges(),
	// thus adding (possibily a lot of) vertexes, though DeleteCWEdges()
	// should be called BEFORE Refine()
	SortVertexes();

	// find index corresponding to given vertex
	a = FindVertex(p1);
	b = FindVertex(p2);

	// if either vertex is not in the mesh, return
	if( a == C2D_VERTEX_NOT_FOUND || b == C2D_VERTEX_NOT_FOUND )
		return false;

	// convert index from 'sortVertexes' to 'vertexes'
	// (InsertConstrEdge() refers to indexes of 'vertexes'
	// array, which is the user input)
	a = GetVIndexFromSortV(a);
	b = GetVIndexFromSortV(b);

	res = DeleteCWEdges(a,b);
	// a,b must always exist!
	ASSERT(res == true);

	return true;
}

// Delete all non-strong constrained edges
// CW of given edge (on its right)
// This routine is used to perform holes in the formed mesh
// (usually before refining it)
bool C2DTriangulate::DeleteCWEdges(long a, long b)
{
	C2DEdge *edgea, *walkingEdge;
	long size;

	// test indexes a and b for validity
	size = vertexes.size();
	if( a < 0 || a >= size || b < 0 || b >= size )
		return false;

	// convert index from 'vertexes' to 'sortVertexes'
	// (InsertConstrEdge() refers to indexes of 'vertexes'
	// array, which is the user input)
	a = GetSortVIndexFromV(a);
	b = GetSortVIndexFromV(b);

	// find an edge with origin a
	edgea = GetVertexEdge(a);

	// find edge joining a with b
	walkingEdge = edgea->Onext();
	while( walkingEdge->MateVertex() != b && walkingEdge != edgea )
		walkingEdge = walkingEdge->Onext();

	// no such edge: return
	if(walkingEdge->MateVertex() != b)
		return true;

	m_lMark = NewIndex();

	// delete edges CW to walkingEdge, if not constrained
	if( walkingEdge->Oprev()->IsHardConstrained() == false)
		DeleteEdges(walkingEdge->Oprev());

	return true;
}

// Recursively delete all non-strong constrained edges
void C2DTriangulate::DeleteEdges(C2DEdge *edge)
{
	C2DEdge *walkingEdge, *startEdge, *tmpEdge;
	bool isHardConstr;

	// this is an error: cannot try to delete a constrained edge
	ASSERT(edge->IsHardConstrained() == false);

	// this is an error: if the edge is already marked for deletion,
	// the function should never have been called
	ASSERT(edge->GetId() != m_lMark);

	// find CCW-most constrained edge
	walkingEdge = edge->Onext();
	while(walkingEdge->IsHardConstrained() == false && walkingEdge != edge)
		walkingEdge = walkingEdge->Onext();

	// position on starting edge (next edge CW)
	startEdge = walkingEdge->Oprev();

	// this is an error: cannot have two consecutive CW constrained edges
	ASSERT(startEdge->IsHardConstrained() == false);

	// mark all edges between constrained edges (if any) for deletion,
	// otherwise all edges around CW
	walkingEdge = startEdge;
	do {
		walkingEdge->SetId(m_lMark);
		walkingEdge = walkingEdge->Oprev();
	} while(walkingEdge->IsHardConstrained() == false && walkingEdge != startEdge);


	// call recursively DeleteEdges(), if mate edge is not already marked
	// or already deleted
	walkingEdge = startEdge;
	do {
		if( walkingEdge->Mate() != NULL )
			if(	walkingEdge->Mate()->GetId() != m_lMark)
				DeleteEdges(walkingEdge->Mate());
		walkingEdge = walkingEdge->Oprev();
	} while(walkingEdge->IsHardConstrained() == false && walkingEdge != startEdge);


	// actually delete radial half-edges
	walkingEdge = startEdge;
	do {
		tmpEdge = walkingEdge->Oprev();
		// must store, otherwise could try to access a deleted edge (if last one)
		isHardConstr = tmpEdge->IsHardConstrained();
		DeleteHalfEdge(walkingEdge);
		walkingEdge = tmpEdge;
	} while(isHardConstr == false && walkingEdge != startEdge);
}

bool operator<(const C2DTriplet &t1, const C2DTriplet &t2)
{
	if( t1.vertex[0] < t2.vertex[0] ) {
		return true;
	}
	else if( t1.vertex[0] == t2.vertex[0] ) {
		if( t1.vertex[1] < t2.vertex[1] ) {
			return true;
		}
		else if( t1.vertex[1] == t2.vertex[1] ) {
			if( t1.vertex[2] < t2.vertex[2] ) {
				return true;
			}
		}
	}

	return false;
}

// return the vertex index associated to given edge #
long C2DEdgeTriangle::GetVertex(unsigned char i)
{
	ASSERT(i<=2);
	return(edges[i]->Vertex());
}




