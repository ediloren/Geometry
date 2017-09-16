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


// Operation3D.h : Non-Manifold Geometry basic operations
// Enrico Di Lorenzo 2002/12/30

#ifndef C3DOPERATION_DEFS
#define C3DOPERATION_DEFS


#include <utility>

#include "MathDefs.h"
#include "Geometry3D.h"

// routine return values
#define C3D_OK						0
#define C3D_POINTA_ON_PLANE			2
#define C3D_POINTB_ON_PLANE			3
#define C3D_BOTH_POINTS_ON_PLANE	4
#define C3D_LINE_DONT_INTERSECT		5
#define C3D_LINE_IS_ON_PLANE		6
#define C3D_POINTS_ARE_COLLINEAR	7
#define C3D_VECTORS_ARE_PARALLEL	8
#define C3D_LINES_ARE_PARALLEL		9
#define C3D_PLANES_ARE_PARALLEL		10
#define C3D_POINT_NOT_ON_LINE		11
#define C3D_INTERNAL_ERROR			1000

// FindQuadrant() return codes
#define C3D_QUADRANT_1				0x00
#define C3D_QUADRANT_2				0x01
#define C3D_QUADRANT_3				0x02
#define C3D_QUADRANT_4				0x03
#define C3D_QUADRANT_ON_X			0x10
#define C3D_QUADRANT_ON_Y			0x20
#define C3D_QUADRANT_ON_0			0x30
#define C3D_QUADRANT_SECTOR			0x0F
#define C3D_QUADRANT_ON				0xF0

// ComputeAngle() return codes
#define C3D_POINT_ON_EDGE			-126
#define C3D_POINT_ON_VERTEX			-127

// IsPointInFace() return codes
#define C3D_POINT_OUT_FACE			0
#define C3D_POINT_IN_FACE			1
#define C3D_POINT_ON_FACE_EDGE		2
#define C3D_POINT_ON_FACE_VERTEX	3

// C3DFaceBSPTree class defines
#define C3D_FACEBSP_NONE		0x00
#define C3D_FACEBSP_INSIDE		0x01
#define C3D_FACEBSP_OUTSIDE		0x02
#define C3D_FACEBSP_BOTH		0x03
#define C3D_FACEBSP_ON			0x04

// IsPointInShell() return codes
#define C3D_POINT_ON_SHELL_MASK		0x03
#define C3D_POINT_ON_SHELL_EDGE		0x01
#define C3D_POINT_ON_SHELL_VERTEX	0x02
#define C3D_POINT_ON_SHELL_FACE		0x03
#define C3D_POINT_IN_SHELL			0x04
#define C3D_POINT_OUT_SHELL			0x08
#define C3D_POINT_UNINITIALIZED		0x10
#define C3D_RAY_INVALID				0x20
// DEBUG WARNING: remove following define
// together with BSP tree routines
#define C3D_POINT_UNDECIDABLE		0x10

// IntersTrianglePlane() return codes
#define C3D_TRIANGLE_ON_PLANE		-1

// Remark: if redefining C3D_OP_TOL, take care that also CMATH_PI had
// at least the same number of decimals
#define C3D_OP_TOL			1E-12


class C3DOperation
{
public:
	C3DOperation() : m_dTol(C3D_OP_TOL), m_vecTol(C3D_OP_TOL, C3D_OP_TOL, C3D_OP_TOL) {}
	C3DOperation(double tol) : m_dTol(tol), m_vecTol(tol, tol, tol) {}

	inline char FindLargestProj(const C3DVector &vect)
	{
		double lx, ly, lz;

		lx = fabs(vect.x);
		ly = fabs(vect.y);
		lz = fabs(vect.z);

		if( lx >= ly )
			if( lx >= lz )
				return C3D_DIR_X;
			else
				return C3D_DIR_Z;
		else
			if( ly >= lz )
				return C3D_DIR_Y;
			else
				return C3D_DIR_Z;
	}

	inline char FindSmallestProj(const C3DVector &vect)
	{
		double lx, ly, lz;

		lx = fabs(vect.x);
		ly = fabs(vect.y);
		lz = fabs(vect.z);

		if( lx <= ly )
			if( lx <= lz )
				return C3D_DIR_X;
			else
				return C3D_DIR_Z;
		else
			if( ly <= lz )
				return C3D_DIR_Y;
			else
				return C3D_DIR_Z;
	}

	// Find most suitable projection plane for given vector.
	// Moreover, returns biggest vector component index
	inline char FindProjPlane(char *i, char *j, const C3DVector &vect)
	{
		double lx, ly, lz;

		lx = fabs(vect.x);
		ly = fabs(vect.y);
		lz = fabs(vect.z);

		ASSERT(lx >= C3D_TOL || ly >= C3D_TOL || lz >= C3D_TOL);

		if( lx >= ly )
			if( lx >= lz ) {
				*i = C3D_DIR_Y;
				*j = C3D_DIR_Z;
				return C3D_DIR_X;
			}
			else {
				*i = C3D_DIR_X;
				*j = C3D_DIR_Y;
				return C3D_DIR_Z;
			}
		else
			if( ly >= lz ) {
				*i = C3D_DIR_X;
				*j = C3D_DIR_Z;
				return C3D_DIR_Y;
			}
			else {
				*i = C3D_DIR_X;
				*j = C3D_DIR_Y;
				return C3D_DIR_Z;
			}
	}

	inline double Dist(C3DVector &p1, C3DVector &p2)
	{
		return Mod(p1-p2);
	}

	inline double PointPlaneDist(C3DVector &point, C3DPlane &plane)
	{
		return ( DotProd(point, plane.m_vecNormal) - plane.m_dD );
	}

	// return point at relative normalized distance t0 on segment point_a-point_b,
	// starting from point_a (t0 = 0)
	// point = (1-t0)*point_a + t0*point_b
	inline C3DVector PointOnSegment(C3DVector &point_a, C3DVector &point_b, double t0)
	{
		C3DVector point;

		ASSERT( t0>=0.0 && t0<=1.0 );

		point.x = (1-t0) * point_a.x + t0 * point_b.x;
		point.y = (1-t0) * point_a.y + t0 * point_b.y;
		point.z = (1-t0) * point_a.z + t0 * point_b.z;

		return point;
	}

	// Given a line in parametric form ([x] = [x0] + t * [dir]),
	// and a point on the line, find the corresponding t
	inline int PointAbscissaOnLine(C3DLine &line, C3DVector &point, double &t)
	{
		char bigProj;

		// find abscissa using largest 'm_vecDir' coordinate
		bigProj = FindLargestProj(line.m_vecDir);
		t = (point-line.m_vecPoint)[bigProj] / line.m_vecDir[bigProj];

		// check that point really lies on the line (within tolerance)
		if( point.x - line.m_vecPoint.x - t * line.m_vecDir.x > m_dTol ||
		    point.y - line.m_vecPoint.y - t * line.m_vecDir.y > m_dTol ||
		    point.z - line.m_vecPoint.z - t * line.m_vecDir.z > m_dTol)
			return C3D_POINT_NOT_ON_LINE;

		return C3D_OK;
	}

	// Compute segment - plane intersection point
	// If intersection exists, return C3D_OK and point 'point',
	// otherwise return status or error condition
	inline int IntersSegmentPlane(C3DVector &point, C3DVector &point_a, C3DVector &point_b, C3DPlane &plane)
	{
		double dista, distb, absdista, absdistb, len, t0;

		dista = PointPlaneDist(point_a, plane);
		distb = PointPlaneDist(point_b, plane);

		len = dista - distb;

		absdista = fabs(dista);
		absdistb = fabs(distb);

		// if both points are on the same half-space side,
		// the segment does not intersect the plane
		if( (dista > m_dTol && distb > m_dTol) || (dista < -m_dTol && distb < -m_dTol) ) {
			return C3D_LINE_DONT_INTERSECT;
		}

		// If dista == 0, point a is on the plane
		if( absdista < m_dTol ) {
			// If both are on the plane, return the first one
			// (either one would be ok) and signal it
			if( absdistb < m_dTol ) {
				point = point_a;
				return C3D_BOTH_POINTS_ON_PLANE;
			}

			point = point_a;
			return C3D_POINTA_ON_PLANE;
		}

		// If distb == 0, point b is on the plane
		if( absdistb < m_dTol ) {
			point = point_b;
			return C3D_POINTB_ON_PLANE;
		}

		// if len == 0, endpoints are coincident
		// (and cannot be both on the plane, or should have
		// been catched before)
		if( len < m_dTol && len > -m_dTol) {
			// no intersection with plane (this condition should
			// have been catched by the test on half space containing points;
			// however, due to combining tolerances (m_dTol) when calculating len,
			// the condition could have been missed
			return C3D_LINE_DONT_INTERSECT;
		}

		// compute relative intersection distance along segment
		t0 = dista / len;
		// determine intersection point at relative distance t0 from point_a
		point = PointOnSegment(point_a, point_b, t0);

		return C3D_OK;
	}

	// Compute line - plane intersection point
	// If intersection exists, return C3D_OK and point 'point',
	// otherwise return error condition
	inline int IntersLinePlane(C3DVector &point, C3DLine &line, C3DPlane &plane,
								double &len_along_line, double &dist_origin_plane)
	{
		double cosalpha;
		C3DVector point_b;

		// Calculate the angle between the line and the plane
		cosalpha = DotProd(plane.m_vecNormal, line.m_vecDir);
		// Compute distance of the plane from the line origin
		// (which is the distance of the line origin from the plane
		// with opposite sign)
		dist_origin_plane = -PointPlaneDist(line.m_vecPoint, plane);

		// if cos(angle) is almost 0, the line is parallel
		// to the plane
		if( cosalpha < m_dTol && cosalpha > -m_dTol) {
			// If line origin is on the plane, then all the line in ON the plane
			if( dist_origin_plane < m_dTol ) {
				return C3D_LINE_IS_ON_PLANE;
			}
			else
				return C3D_LINE_DONT_INTERSECT;
		}

		// compute lenght along line, starting from line origin,
		// to intersection point
		len_along_line = dist_origin_plane / cosalpha;

		// compute intersection point
		point = line.m_vecPoint + line.m_vecDir * len_along_line;

		return C3D_OK;
	}

	// Compute plane - plane intersection line
	// If intersection exists, return C3D_OK and line 'line',
	// otherwise return error condition
	inline int IntersPlanePlane(C3DLine &line, C3DPlane &plane1, C3DPlane &plane2)
	{
		char bigProj;
		double det;

		// line direction is perpendicular to both plane normals
		line.m_vecDir = CrossProd(plane1.m_vecNormal, plane2.m_vecNormal);
		line.m_vecDir /= Mod(line.m_vecDir);

		// check for parallelism of planes
		if( Mod(line.m_vecDir) < m_dTol ) {
			return C3D_PLANES_ARE_PARALLEL;
		}

		// Now must find a point on the line as origin.
		// If n1 and n2 are plane normals and -d1, -d2 their
		// distances from origin, we can find a point
		// on the intersection line by imposing to point x0
		// to lie on both planes:
		//
		// n1*x0 = -d1              (1)
		// n2*x0 = -d2
		//
		// This system is over-determined (2 equations, 3 variables),
		// since identifies a line; to find a point, simply fix one
		// coordinate. However, we must be careful in the coordinate
		// choice; e.g., if line.m_vecDir has null z-component,
		// we cannot impose z0 = 0: in this case, n1 x n2 gave
		// null z component and this means that from
		// [i j k; n1x n1y n1z; n2x n2y n2z] we have
		// k*(n1x*n2y - n1y*n2x) = k*0;
		// but [n1x n1y; n2x n2y] = 0 is the determinant of system (1)
		// when z0 = 0.
		// To avoid this problem, let's null line's largest element
		bigProj = FindLargestProj(line.m_vecDir);

		if(bigProj == C3D_DIR_Z) {
			det = plane1.m_vecNormal.x * plane2.m_vecNormal.y - plane1.m_vecNormal.y * plane2.m_vecNormal.x;
			line.m_vecPoint.x = (plane1.m_dD * plane2.m_vecNormal.y - plane2.m_dD * plane1.m_vecNormal.y) / det;
			line.m_vecPoint.y = (-plane1.m_dD * plane2.m_vecNormal.x + plane2.m_dD * plane1.m_vecNormal.x) / det;
			line.m_vecPoint.z = 0.0;
		}
		else if(bigProj == C3D_DIR_Y) {
			det = - plane1.m_vecNormal.x * plane2.m_vecNormal.z + plane1.m_vecNormal.z * plane2.m_vecNormal.x;
			line.m_vecPoint.x = (- plane1.m_dD * plane2.m_vecNormal.z + plane2.m_dD * plane1.m_vecNormal.z) / det;
			line.m_vecPoint.z = (plane1.m_dD * plane2.m_vecNormal.x - plane2.m_dD * plane1.m_vecNormal.x) / det;
			line.m_vecPoint.y = 0.0;
		}
		else {
			det = plane1.m_vecNormal.y * plane2.m_vecNormal.z - plane1.m_vecNormal.z * plane2.m_vecNormal.y;
			line.m_vecPoint.y = (plane1.m_dD * plane2.m_vecNormal.z - plane2.m_dD * plane1.m_vecNormal.z) / det;
			line.m_vecPoint.z = (- plane1.m_dD * plane2.m_vecNormal.y + plane2.m_dD * plane1.m_vecNormal.y) / det;
			line.m_vecPoint.x = 0.0;
		}

		return C3D_OK;
	}
/*
	// Compute triangle - plane intersection points on triangle boundary;
	// since there could be one, two or none, return the number of points
	// (0, 1 or 2) according to the condition.
	// Remark: if only one point, the position of this point is copied
	// in both members of the array 'points'
	inline int IntersTrianglePlane(C3DPlane &plane, C3DFace &face, C3DVector points[2])
	{
		C3DList<C3DLoop*> *looplist;
		C3DList<C3DEdgeUse*> *edgeuselist;
		C3DList<C3DEdgeUse*>::iterator iteu;
		char pointindex;
		C3DVector startpt, endpt;
		int res;

		// verify that face has only three edges; this means
		// only one loop with exactly three edges

		looplist = &(face.loopList);
		ASSERT(looplist->size() == 1);
		edgeuselist = &( (*(looplist->begin()))->edgeUseList);
		ASSERT(edgeuselist->size() == 3);

		if( plane.m_vecNormal == face.plane.m_vecNormal && fabs(plane.m_dD - face.plane.m_dD) < C3D_TOL ) {
			return C3D_TRIANGLE_ON_PLANE;
		}

		for(iteu = edgeuselist->begin(), pointindex = 0; iteu != edgeuselist->end(); iteu++) {
			// get edge vertexes
			startpt = (*iteu)->GetFirstVertexUse()->pVertex->pos;
			endpt = (*iteu)->GetSecondVertexUse()->pVertex->pos;
			// intersect edge with plane
			res = IntersSegmentPlane(points[pointindex], startpt, endpt, plane);
			// process return status (type of intersection)
			if( res == C3D_POINTA_ON_PLANE || res == C3D_OK )
				pointindex++;
			if( res == C3D_BOTH_POINTS_ON_PLANE ) {
				points[pointindex] = startpt;
				pointindex++;
			}
			// if( res == C3D_LINE_DONT_INTERSECT || res == C3D_POINTB_ON_PLANE),
			// do nothing (if point b is on plane, point a of next segment will be too,
			// so only one point will be inserted as intersection)
		}

		// verify that only zero, one or two points have been identified
		ASSERT(pointindex >= 0 && pointindex <= 2);

		// if two points have been identified,verify that they are not
		// coincident (within tolerance). This could happen if the plane
		// is intersecting the triangle near a small-angled vertex.
		//
		//        /\ 20 degrees
		//     s / |\
		//      / h| \
		//  ---o------o---
		//    /    l    \
		//   /          \
		//
		// Note that if all trangles are good quality (i.e. smallest
		// angle >= 20 degrees), the maximum distance from the vertex
		// at which this could happen is h = (l/2) / tan(20 / 2) = l * 2.84,
		// therefore if l <= m_dTol, h can be up to 2.84 times m_dTol
		// and still the intersection point does not appear to be on the
		// intersection plane (within tolerance).
		//
		// We can therefore assume that if distance l is <= m_dTol and
		// distance s from the corresponding vertex is smaller than
		// 4 * m_dTol, we have only one intersection point on the
		// vertex. Otherwise, this is an error, which could lead to
		// the formation of degenerate triangles: signal an error.
		if( pointindex == 2) {
			if( Mod(points[0] - points[1]) <= m_dTol ) {
				// find vertex
				for(iteu = edgeuselist->begin(); iteu != edgeuselist->end(); iteu++) {
					// get edge vertexes
					startpt = (*iteu)->GetFirstVertexUse()->pVertex->pos;

					if( Mod(startpt - points[0]) <= m_dTol * 4.0 &&
						Mod(startpt - points[1]) <= m_dTol * 4.0) {

						pointindex = 1;
						points[0] = startpt;

						return pointindex;
					}
				}

				GlbErrorMsg("ERROR: Numerical tolerance problem when intersecting plane and triangle");
			}
		}
		else if( pointindex == 1) {
			// if only one point, copy the same location
			// also in second point position (two coincident points)
			points[1] = points[0];
		}

		return pointindex;
	}

	class C3DIsect
	{
	public:
		C3DVector isectPt[2];
		double isectAbscissa[2];
		int isectPtsNum;

		void Order()
		{
			C3DVector tmpv;
			double tmpd;

			if( isectAbscissa[0] > isectAbscissa[1] ) {
				tmpv = isectPt[0];
				tmpd = isectAbscissa[0];
				isectPt[0] = isectPt[1];
				isectAbscissa[0] = isectAbscissa[1];
				isectPt[1] = tmpv;
				isectAbscissa[1] = tmpd;
			}
		}

	};

	inline bool Smaller(C3DIsect &isect1, C3DIsect &isect2)
	{
		return( isect1.isectAbscissa[0] + m_dTol < isect2.isectAbscissa[0] &&
			    isect1.isectAbscissa[1] + m_dTol < isect2.isectAbscissa[0] );
	}

	// Compute triangle - triangle intersection, determining a new point or a new edge
	// to be inserted in the triangles, along the line of intersection
	inline void IntersTriangleTriangle(C3DFace &tri1, C3DFace &tri2)
	{
		C3DList<C3DLoop*> *looplist1, *looplist2;
		C3DList<C3DEdgeUse*> *edgeuselist1, *edgeuselist2;
		C3DIsect isect[2];
		C3DIntersection faceIsect;
		C3DLine isectLine;
		int i;
		char res;

		// verify that each triangular face has only three edges; this means
		// only one loop with exactly three edges

		looplist1 = &(tri1.loopList);
		ASSERT(looplist1->size() == 1);
		edgeuselist1 = &( (*(looplist1->begin()))->edgeUseList);
		ASSERT(edgeuselist1->size() == 3);

		looplist2 = &(tri2.loopList);
		ASSERT(looplist2->size() == 1);
		edgeuselist2 = &( (*(looplist2->begin()))->edgeUseList);
		ASSERT(edgeuselist2->size() == 3);

		// intersect triangles with each other's plane
		isect[0].isectPtsNum = IntersTrianglePlane(tri1.plane, tri2, isect[0].isectPt);
		isect[1].isectPtsNum = IntersTrianglePlane(tri2.plane, tri1, isect[1].isectPt);

		// if either of the triangle did not intersect the plane
		// on which the other triangle lies, intersection is void
		if( isect[0].isectPtsNum == 0 || isect[1].isectPtsNum == 0 )
			return;

		// if the triangles are coplanar, return
		if( isect[0].isectPtsNum == C3D_TRIANGLE_ON_PLANE || isect[1].isectPtsNum == C3D_TRIANGLE_ON_PLANE )
			return;

		// otherwise must find if there is a common intersection
		//

		// intersect triangles support planes
		IntersPlanePlane(isectLine, tri1.plane, tri2.plane);

		// find intersection point abscissa along intersection line
		for(i=0; i<2; i++) {
			res = PointAbscissaOnLine(isectLine, isect[0].isectPt[i], isect[0].isectAbscissa[i]);
			// point must be on the intersection line; if not, error!
			ASSERT( res == C3D_OK);
			res = PointAbscissaOnLine(isectLine, isect[1].isectPt[i], isect[1].isectAbscissa[i]);
			// point must be on the intersection line; if not, error!
			ASSERT( res == C3D_OK);
		}

		// order intersection segment endpoints
		isect[0].Order();
		isect[1].Order();

		// no intersection within tolerance, do nothing
		if( Smaller(isect[0], isect[1]) || Smaller(isect[1], isect[0]) ) {
			return;
		}

		// find 1D intersection along line
		if( isect[0].isectAbscissa[0] > isect[1].isectAbscissa[0] ) {
			faceIsect.p1 = isect[0].isectPt[0];
			if( isect[0].isectAbscissa[1] < isect[1].isectAbscissa[1] )
				faceIsect.p2 = isect[0].isectPt[1];
			else
				faceIsect.p2 = isect[1].isectPt[1];
		}
		else  {
			faceIsect.p1 = isect[1].isectPt[0];
			if( isect[1].isectAbscissa[1] < isect[0].isectAbscissa[1] )
				faceIsect.p2 = isect[1].isectPt[1];
			else
				faceIsect.p2 = isect[0].isectPt[1];
		}

		// if points are coincident within tolerance, merge them
		if( Mod(faceIsect.p1 - faceIsect.p2) < m_dTol ) {
			faceIsect.p2 = faceIsect.p1;
		}

		// insert points in faces intersection lists
		tri1.intersectList.push_back(faceIsect);
		tri2.intersectList.push_back(faceIsect);
	}


	// Compute the area of a polygon, using:
	// 1. projection of the polygon on a plane (ignoring one of the coordinates, where
	//    the component of the plane normal n = (nx,ny,nz) is biggest, i.e. Max(nx,ny,nz))
	// 2. Stone's method (see Graphics Gems II, I.1) in 2D
	// 3. 3D area is recovered using an area scaling factor:
	//    the ratio of areas for the projected polygon Projc(W) and original planar polygon W
	//    with normal n = (nx,ny,nz) is  Area( Projc(W) ) / Area (W) =  Max(nx,ny,nz) / Mod(n)
	double AreaOfPolygon(C3DList<C3DEdgeUse*> &edgeuses, C3DVector normal)
	{
		C3DList<C3DEdgeUse*>::iterator it;
		double area;
		char c1, c2, c3;

		// find a suitable projection plane
	 	c3 = FindProjPlane(&c1, &c2, normal);
		// find the area of the projection
		area = AreaOfPolygonProj(edgeuses, c1, c2);
		// and recover the full area scaling per cos(alpha)
		area /= normal[c3];

		return area;
	}


	// Compute the area of a polygon projection, using Stone's method (see Graphics Gems II, I.1) in 2D
	double AreaOfPolygonProj(C3DList<C3DEdgeUse*> &edgeuses, char c1, char c2)
	{
		C3DList<C3DEdgeUse*>::iterator it;
		C3DVector point1, point2;
		double area;

		area = 0.0;
		// scan all points
		for(it = edgeuses.begin(); it != edgeuses.end(); it++) {

			// get edge vertexes coordinates;
			// note that edgeuse direction is important!
			point1 = (*it)->GetFirstVertexUse()->pVertex->pos;
			point2 = (*it)->GetSecondVertexUse()->pVertex->pos;

			// calculate 2 * the signed area of the trapezoids determined by the edges
			// (see PlaneFromPolygon() for more information)
			area += point1[c1] * point2[c2] - point2[c1] * point1[c2];
		}

		area /= 2.0;

		return area;
	}
*/

	// Compute plane of a polygon, using Newell's method (see e.g. graphics gems III, V.5)
	// Plane normal direction is so that N x point1-point3
	// is directed inside the polygon
	//
	//         3                    O
	//         |\                   O
	//    N    | \                  O
	//     \   |  \                 O
	//      \  |   \                O
	//       \ |    \               O
	//        \|     \              O
	//         1------2             O
	//
	// If plane equation can be found, return C3D_OK and plane 'plane',
	// otherwise return error condition
	int PlaneFromPolygon(C3DPlane &plane, double &mod, C3DVector *vertexes, int numvertexes)
	{
		C3DVector normal(0.0), midpoint(0.0), point1, point2;
		int pointnum, i;

		pointnum = 0;

		// scan all points
		for(i = 0; i < numvertexes; i++) {

			// get edge coordinates;
			point1 = vertexes[i];
			if(i < numvertexes-1) {
				point2 = vertexes[i+1];
			}
			else {
				point2 = vertexes[0];
			}

			// calculate components of normal vector, using the
			// fact that the area of the trapezoids projected on the
			// xy, xz, yz planes is proportional to the z, y, x
			// components of the normal
			normal.x += (point1.y - point2.y) * (point1.z + point2.z);
			normal.y += (point1.z - point2.z) * (point1.x + point2.x);
			normal.z += (point1.x - point2.x) * (point1.y + point2.y);

			// calculates the mean point of the polygon
			midpoint += point1;

			pointnum++;
		}

		mod = Mod(normal);

		// if the result of the cross product is too small,
		// can only mean that the points are almost collinear
		if( mod < m_dTol ) {
			return C3D_POINTS_ARE_COLLINEAR;
		}

		// unify normal vector
		normal = normal / mod;

		// compute mean point
		midpoint /= pointnum;

		plane.m_vecNormal = normal;
		plane.m_dD = DotProd(normal, midpoint);

		return C3D_OK;
	}

/*
	// Compute plane of a polygon, using Newell's method (see e.g. graphics gems III, V.5)
	// Plane normal direction is so that N x point1-point3
	// is directed inside the polygon
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
	int PlaneFromPolygon(C3DPlane &plane, double &mod, C3DList<C3DEdgeUse*> &edgeuses)
	{
		C3DList<C3DEdgeUse*>::iterator it;
		C3DVector normal, midpoint, point1, point2;
		int pointnum;

		pointnum = 0;

		// scan all points
		for(it = edgeuses.begin(); it != edgeuses.end(); it++) {

			// get edge vertexes coordinates;
			// note that edgeuse direction is important!
			point1 = (*it)->GetFirstVertexUse()->pVertex->pos;
			point2 = (*it)->GetSecondVertexUse()->pVertex->pos;

			// calculate components of normal vector, using the
			// fact that the area of the trapezoids projected on the
			// xy, xz, yz planes is proportional to the z, y, x
			// components of the normal
			normal.x += (point1.y - point2.y) * (point1.z + point2.z);
			normal.y += (point1.z - point2.z) * (point1.x + point2.x);
			normal.z += (point1.x - point2.x) * (point1.y + point2.y);
			// remark: the above is the original Newell's method. However, notice
			// that for next edge, the old point2 becomes the new point 1
			// therefore, developing the multiplication, we have (for instance on normal.z),
			// and considering three points 1, 2, 3:
			// first step:
			// normal.z += point1.x*point1.y - point1.x*point2.y - point2.x*point1.y - point2.x*point2.y
			// next step:
			// normal.z += point2.x*point2.y - point2.x*point3.y - point3.x*point2.y - point3.x*point3.y
			// so the third member of the first equation simplifies with the first member of the second
			// (always point2.x*point2.y) so we save two multiplications.
			// This is by the way how the Stone's method (see Graphics Gems II, I.1) in 2D works to
			// calculate the area of a polygon (remark: we should divide by two if we wanted to get the real area)
			// Remark: numerically, the above method appears more stable!! Using the below shortcut leads
			// to slight differences in the normal calculation, that in the end let the point on triangle
			// test fail!!
			//normal.x += point1.y*point2.z - point1.z*point2.y;
			//normal.y += point1.z*point2.x - point2.z*point1.x;
			//normal.z += point1.x*point2.y - point2.x*point1.y;

			// calculates the mean point of the polygon
			midpoint += point1;
			pointnum++;
		}

		mod = Mod(normal);

		// if the result of the cross product is too small,
		// can only means that the points are almost collinear
		if( mod < m_dTol ) {
			return C3D_POINTS_ARE_COLLINEAR;
		}

		// unify normal vector
		normal = normal / mod;

		// compute mean point
		midpoint /= pointnum;

		plane.m_vecNormal = normal;
		plane.m_dD = DotProd(normal, midpoint);

		return C3D_OK;
	}

	// Find the quadrant in which given point lies w.r.t. the axis origin
	// If the point is on the cartesian axes, signal it with flag
	// (in this case the point belongs to the positive halfplanes)
	inline char FindQuadrant(double x, double y)
	{
		if( x > m_dTol )
			if( y > m_dTol )
				return C3D_QUADRANT_1;
			else if( y < -m_dTol )
				return C3D_QUADRANT_4;
			else
				return C3D_QUADRANT_1 | C3D_QUADRANT_ON_X;
		else if( x < -m_dTol )
			if( y > m_dTol )
				return C3D_QUADRANT_2;
			else if( y < -m_dTol )
				return C3D_QUADRANT_3;
			else
				return C3D_QUADRANT_2 | C3D_QUADRANT_ON_X;
		else
			if( y > m_dTol )
				return C3D_QUADRANT_1 | C3D_QUADRANT_ON_Y;
			else if( y < -m_dTol )
				return C3D_QUADRANT_4 | C3D_QUADRANT_ON_Y;
			else
				return C3D_QUADRANT_ON_0;
	}

	// Calculate the 'angle' spanned by a point moving from
	// quadrant1 to quadrant2
	inline char ComputeAngle(char quadrant1, char quadrant2, double x1, double y1,
							double x2, double y2)
	{
		// array of 'angles' according to quadrant position
		static char angles[4][4] = { {0, 1, 2, -1}, {-1, 0, 1, 2}, {2, -1, 0, 1}, {1, 2, -1, 0} };
		char q1, q2, on1, on2, res;
		double dist;
		C2DVector v1_1(-x2, -y2), v1_2(x1-x2, y1-y2);
		C2DVector v2_1(-x1, -y1), v2_2(x2-x1, y2-y1);

		q1 = quadrant1 & C3D_QUADRANT_SECTOR;
		q2 = quadrant2 & C3D_QUADRANT_SECTOR;
		on1 = quadrant1 & C3D_QUADRANT_ON;
		on2 = quadrant2 & C3D_QUADRANT_ON;

		if( (on1 == C3D_QUADRANT_ON_X && on2 == C3D_QUADRANT_ON_X) ||
			(on1 == C3D_QUADRANT_ON_Y && on2 == C3D_QUADRANT_ON_Y) )
			return C3D_POINT_ON_EDGE;

		if( on1 == C3D_QUADRANT_ON_0 || on2 == C3D_QUADRANT_ON_0 )
			return C3D_POINT_ON_VERTEX;

		res = angles[q1][q2];

		// if points jump to opposite quadrant, should compute if movement
		// is positive or negative
		if( res == 2 ) {
			// Compute distance from origin; must use symmetric forumla, otherwise
			// the same segment considered from one vertex to the other or
			// vice-versa could lead to two (slightly) different results
			dist = ( CrossProd(v1_1, v1_2) / Mod(v1_2) - CrossProd(v2_1, v2_2) / Mod(v2_2)) / 2.0;

			// if distance is almost zero, point is on edge
			if( -m_dTol < dist && dist < m_dTol )
				return C3D_POINT_ON_EDGE;

			// if distance is negative, angle is spanned in the opposite
			// direction
			if( dist < 0)
				res = -2;
		}

		return res;
	}

	C3DVector TransformFunction(C3DVector &point)
	{
		C3DVector v;

		if(m_iTranOperation == C3D_TRANS_TRANSLATE)
			return TranslatePoint(point, m_vecShift);
		else if(m_iTranOperation == C3D_TRANS_SCALE)
			return ScalePoint(point, m_vecScale);
		else if(m_iTranOperation == C3D_TRANS_ROTATEX)
			return RotatePointX(point);
		else if(m_iTranOperation == C3D_TRANS_ROTATEY)
			return RotatePointY(point);
		else if(m_iTranOperation == C3D_TRANS_ROTATEZ)
			return RotatePointZ(point);
		else {
			// always error
			ASSERT(NULL);
			// return dummy vector
			return v;
		}
	}

	int GetTransformType()
	{
		return m_iTranOperation;
	}

	void SetTranslationShift(C3DVector &shift)
	{
		m_vecShift = shift;
		m_iTranOperation = C3D_TRANS_TRANSLATE;
	}

	C3DVector TranslatePoint(C3DVector &point, C3DVector &shift)
	{
		return (point + shift);
	}

	void SetScaleVector(C3DVector &scale)
	{
		m_vecScale = scale;
		m_iTranOperation = C3D_TRANS_SCALE;
	}

	C3DVector GetScaleVector()
	{
		return m_vecScale;
	}

	C3DVector ScalePoint(C3DVector &point, C3DVector &scale)
	{
		return (point * scale);
	}

	// angle is in degrees
	void SetRotationAngleAndOperation(double angle, int operation)
	{
		double alpha;

		alpha = angle * CMATH_PI / 180.0;
		m_dSinAngle = sin(alpha);
		m_dCosAngle = cos(alpha);

		m_iTranOperation = operation;
	}

	// angle is in degrees
	C3DVector RotatePointX(C3DVector &point, double angle)
	{
		SetRotationAngleAndOperation(angle, C3D_TRANS_ROTATEX);
		return RotatePointX(point);
	}

	C3DVector RotatePointX(C3DVector &point)
	{
		C3DVector res;

		res.x = point.x;
		res.y = point.y * m_dCosAngle - point.z * m_dSinAngle;
		res.z = point.y * m_dSinAngle + point.z * m_dCosAngle;

		return res;
	}

	// angle is in degrees
	C3DVector RotatePointY(C3DVector &point, double angle)
	{
		SetRotationAngleAndOperation(angle, C3D_TRANS_ROTATEY);
		return RotatePointY(point);
	}

	C3DVector RotatePointY(C3DVector &point)
	{
		C3DVector res;

		res.x = point.x * m_dCosAngle + point.z * m_dSinAngle;
		res.y = point.y;
		res.z = - point.x * m_dSinAngle + point.z * m_dCosAngle;

		return res;
	}

	// angle is in degrees
	C3DVector RotatePointZ(C3DVector &point, double angle)
	{
		SetRotationAngleAndOperation(angle, C3D_TRANS_ROTATEZ);
		return RotatePointZ(point);
	}

	// angle is in degrees
	C3DVector RotatePointZ(C3DVector &point)
	{
		C3DVector res;

		res.x = point.x * m_dCosAngle - point.y * m_dSinAngle;
		res.y = point.x * m_dSinAngle + point.y * m_dCosAngle;
		res.z = point.z;

		return res;
	}

	// This routine returns a unique marker to use to mark
	// geometry data (e.g. to signal they have already been visited)
	long GetMarker()
	{
		return m_lMark++;
	}

	inline bool AreBBoxOverlapping(C3DBBox &bbox1, C3DBBox &bbox2)
	{
		// check overlapping within tolerance; tol is a positive, small vector
		return( bbox1.min_point - m_vecTol < bbox2.max_point && bbox2.min_point - m_vecTol < bbox1.max_point );
	}

	inline bool IsPointInBBox(C3DBBox &bbox, C3DVector &point)
	{
		// check overlapping within tolerance
		return( point.x > bbox.min_point.x - m_dTol &&
				point.x < bbox.max_point.x + m_dTol &&
				point.y > bbox.min_point.y - m_dTol &&
				point.y < bbox.max_point.y + m_dTol &&
				point.z > bbox.min_point.z - m_dTol &&
				point.z < bbox.max_point.z + m_dTol );
	}

	inline bool IsCartLineInBBox(C3DBBox &bbox, C3DVector &point, char c1, char c2)
	{
		// check overlapping within tolerance
		return( point[c1] > bbox.min_point[c1] - m_dTol &&
				point[c1] < bbox.max_point[c1] + m_dTol &&
				point[c2] > bbox.min_point[c2] - m_dTol &&
				point[c2] < bbox.max_point[c2] + m_dTol );
	}

	int PlaneFromThreePoints(C3DPlane &plane, C3DVector &point1,
							C3DVector &point2, C3DVector &point3);
	int PlaneFromTwoVectorsAndPoint(C3DPlane &plane, C3DVector &vector1,
									C3DVector &vector2, C3DVector &point);
	int LineToLineClosestPoints(C3DVector &point1, C3DVector &point2,
								C3DLine &line1, C3DLine &line2);
	void IntersectTriangFaces(C3DFace &face1, C3DFace &face2);
	char IsPointInFace(C3DVector &point, C3DFace *face);
*/
protected:
	double    m_dTol;
	C3DVector m_vecTol;

	C3DVector m_vecShift, m_vecScale;
	double m_dCosAngle, m_dSinAngle;
	int m_iTranOperation;

	static long m_lMark;
};

/*
// Cartesian Binary Space Partition tree,
// for accelerating shell intersection and
// point containment classification
class C3DCartBSPTree
{
protected:

	// tree node
	class Node
	{
	public:
		C3DList<C3DFace*> faceList;
		double divPlane;
		char divDir;
		Node *inside;
		Node *outside;

		Node() : inside(NULL), outside(NULL), divDir(C3D_DIR_NONE), divPlane(0.0)
		{
		}

		~Node()
		{
			// recursively delete tree
			if(inside != NULL)
				delete inside;
			if(outside != NULL)
				delete outside;
		}
	};

	// tree variables
	Node *root;
	C3DBBox globalBbox;
	long m_lMark;
	C3DOperation m_clsOp;
	double    m_dTol;
	C3DLine m_clsRay;

public:
	C3DCartBSPTree() : root(NULL) {}

	~C3DCartBSPTree()
	{
		// delete tree
		if(root != NULL)
			delete root;
	}

	void InitBSPTree(C3DShell &shell, double tolerance);
	char IsPointInShell(C3DVector &point);

protected:
	void init(C3DFace *pFace)
	{
		root = new Node;
	}

	void clear()
	{
		if( root != NULL ) {
			delete root;
			root = NULL;
		}
	}

	void insert(Node *pNode, C3DVector bboxext);
	std::pair<char, long> fireRay(Node *pNode, double min_dist, double max_dist);
	// DEBUG
	std::pair<char, long> visitTree(Node *pNode, double min_dist, double max_dist, int level);

	// find min and max extension of face along direction 'dir'
	std::pair<double, double> getMinMax(C3DFace *pFace, char dir)
	{
		C3DList<C3DLoop*>::iterator itl;
		C3DList<C3DEdgeUse*>::iterator iteu;
		std::pair<double, double> minMax;
		C3DVector point;

		// init minMax
		//

		// face cannot be empty
		ASSERT(pFace->loopList.size() > 0);

		itl = pFace->loopList.begin();

		if((*itl)->edgeUseList.size() > 0) {
			iteu = (*itl)->edgeUseList.begin();
			point = (*iteu)->GetFirstVertexUse()->pVertex->pos;
		}
		else if( (*itl)->pVertexUse != NULL ) {
			point = (*itl)->pVertexUse->pVertex->pos;
		}
		else
			// loop cannot be empty
			ASSERT(false);

		minMax.first = point[dir];
		minMax.second = point[dir];

		// scan all current face 'pFace' vertexes
		// and test them against the face stored in current node
		for(itl = pFace->loopList.begin(); itl != pFace->loopList.end(); itl++) {
			for(iteu = (*itl)->edgeUseList.begin(); iteu != (*itl)->edgeUseList.end(); iteu++) {

				// get edge first end point
				point = (*iteu)->GetFirstVertexUse()->pVertex->pos;

				if(point[dir] < minMax.first)
					minMax.first = point[dir];

				if(point[dir] > minMax.second)
					minMax.second = point[dir];
			}
		}

		return minMax;
	}

};

// Face-induced Binary Space Partition tree,
// point vs. shell classification (in, out or on)
class C3DFaceBSPTree
{
protected:
	// tree node
	class Node
	{
	public:
		C3DFace *face;
		unsigned long num_inside;
		Node *inside;
		unsigned long num_outside;
		Node *outside;

		Node()
		{
		}

		Node(C3DFace *pFace) :
		face(pFace), num_inside(0), inside(NULL), num_outside(0), outside(NULL)
		{
		}

		~Node()
		{
			// recursively delete tree
			if(inside != NULL)
				delete inside;
			if(outside != NULL)
				delete outside;
		}
	};

	// tree variables
	Node *root;

public:
	C3DFaceBSPTree() : root(NULL) {}

	~C3DFaceBSPTree()
	{
		// delete tree
		if(root != NULL)
			delete root;
	}

	void InitPointInShell(C3DShell &shell, double tolerance);
	char IsPointInShell(C3DVector &point, CString name = "")
	{
		FILE *fp = NULL;
		char res;

		if( name.IsEmpty() == false ) {
			fp = fopen( (LPCTSTR)name, "w");

			fprintf(fp, "0 Tree hierarchy '%s'\r\n", (LPCTSTR)name);
			fprintf(fp, "* Searched point\r\n", (LPCTSTR)name);
			fprintf(fp, "Q  Point  %f %f %f  %f %f %f  %f %f %f  %f %f %f \r\n",
				point.x - 0.1, point.y - 0.1, point.z,
				point.x - 0.1, point.y + 0.1, point.z,
				point.x + 0.1, point.y + 0.1, point.z,
				point.x + 0.1, point.y - 0.1, point.z );
			fprintf(fp, "Q  Point  %f %f %f  %f %f %f  %f %f %f  %f %f %f \r\n",
				point.x, point.y - 0.1, point.z - 0.1,
				point.x, point.y - 0.1, point.z + 0.1,
				point.x, point.y + 0.1, point.z + 0.1,
				point.x, point.y + 0.1, point.z - 0.1);
			fprintf(fp, "Q  Point  %f %f %f  %f %f %f  %f %f %f  %f %f %f \r\n",
				point.x - 0.1, point.y, point.z - 0.1,
				point.x - 0.1, point.y, point.z + 0.1,
				point.x + 0.1, point.y, point.z + 0.1,
				point.x + 0.1, point.y, point.z - 0.1 );
		}


		res = find(point, fp);

		if( fp != NULL )
			fclose(fp);

		return(res);
	}

protected:
	void init(C3DFace *pFace)
	{
		root = new Node(pFace);
	}

	void clear()
	{
		if( root != NULL ) {
			delete root;
			root = NULL;
		}
	}

	// Remark: before starting inserting faces, must be sure
	// that face references and inOut statuses have been all cleared!
	void insert(C3DFace *pFace)
	{
		// if first face in tree, must init
		if(root == NULL)
			init(pFace);
		else
			insert(root, pFace);
	}

	void insert(Node *node, C3DFace *pFace);

	char find(C3DVector &point, FILE *fp)
	{
		if(root == NULL)
			return C3D_POINT_OUT_SHELL;
		return( find(root, point, fp));
	}

	char find(Node *node, C3DVector &point,  FILE *fp);

	double    m_dTol;
};
*/

#endif //C3DOPERATION_DEFS

