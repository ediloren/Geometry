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


// Vector2D.h : definition of 2D vector class
// Enrico Di Lorenzo 2003/05/05

#ifndef C2D_VECTOR2D_DEFS
#define C2D_VECTOR2D_DEFS

#include <math.h>

//#include "Global.h"
#include "GeoGlobal.h"
#include "Vector2D_float.h"

/////////////////////////////////////////////////////
// 2D vector class
//

#define C2D_ON_SEGMENT				0
#define C2D_RIGHT_OF				1
#define C2D_LEFT_OF					2

#define C2D_ON_CIRCLE				16
#define C2D_OUT_CIRCLE				17
#define C2D_IN_CIRCLE				18
#define C2D_COLINEAR_POINTS			19

#define C2D_POS_HALFPLANE			32
#define C2D_NEG_HALFPLANE			33
#define C2D_ON_ORIGIN				34

#define C2D_OK						64
#define C2D_SEGMENT_DONT_INTERSECT	65
#define C2D_BOTH_POINTS_ON_LINE		66
#define C2D_POINT_ON_LINE			67

// declared classes
class C2DVector;
class C2DLine;


// vector class specialized for 2d operations
class C2DVector
{
public:
	double			x;
	double			y;

//	C2DVector() : x(0.0), y(0.0) {}
	C2DVector() {}

	C2DVector(double xcoord, double ycoord) :
				x(xcoord), y(ycoord) {}

	C2DVector(C2DVector p1, C2DVector p2) :
				x(p2.x - p1.x), y(p2.y - p1.y) {}

	C2DVector(C2DVector_float floatVector) :
				x(floatVector.x), y(floatVector.y) {}

	friend double DotProd(const C2DVector &vect1, const C2DVector &vect2);
    friend double Mod(const C2DVector &vect);
	friend double CrossProd(const C2DVector &vect1, const C2DVector &vect2);
	friend bool FindCircle3Points(C2DVector p1, C2DVector p2, C2DVector p3,
								C2DVector *center, double *radius2 );

	inline void pos(double xcoord, double ycoord)
	{
		x = xcoord;
		y = ycoord;
	}

	inline void Invert()
	{
		x = -x;
		y = -y;
	}

	inline void Normalize()
	{
	    double mod;

	    mod = Mod();

	    if(mod != 0.0) {
            *this /= Mod();
	    }
//	    else {
//			GlbSysErrorMsg("C3DVector Normalize() called on a zero-length vector");
//	    }

	}

	inline C2DVector operator-()
	{
		C2DVector newPoint;

		newPoint.x = -x;
		newPoint.y = -y;

		return(newPoint);
	}

	inline bool operator==(C2DVector &point)
	{
		return( (*this - point).Mod() < C3D_TOL );
	}

	inline void operator=(const C2DVector &point)
	{
		x = point.x;
		y = point.y;
	}

	inline void operator=(const C2DVector_float &point)
	{
		x = point.x;
		y = point.y;
	}

	inline void operator+=(C2DVector &point)
	{
		x += point.x;
		y += point.y;
	}

	inline C2DVector operator+(const C2DVector &point)
	{
		C2DVector newPoint;

		newPoint.x = x + point.x;
		newPoint.y = y + point.y;

		return newPoint;
	}

	inline void operator-=(C2DVector &point)
	{
		x -= point.x;
		y -= point.y;
	}

	inline C2DVector operator-(const C2DVector &point)
	{
		C2DVector newPoint;

		newPoint.x = x - point.x;
		newPoint.y = y - point.y;

		return newPoint;
	}

	inline C2DVector operator*(const double &scalar)
	{
		C2DVector newPoint;

		newPoint.x = x * scalar;
		newPoint.y = y * scalar;

		return newPoint;
	}

	inline C2DVector operator*(const C2DVector &point)
	{
		C2DVector newPoint;

		newPoint.x = x * point.x;
		newPoint.y = y * point.y;

		return newPoint;
	}

	inline void operator*=(const double &scalar)
	{
		x *= scalar;
		y *= scalar;
	}

	inline C2DVector operator/(const double &scalar)
	{
		C2DVector newPoint;

		if( scalar == 0 ) {
//			GlbSysErrorMsg("C2DVector operator/ tried to divide by zero");
			return *this;
		}

		newPoint.x = x / scalar;
		newPoint.y = y / scalar;

		return newPoint;
	}

	inline C2DVector &operator/=(const double &scalar)
	{
		if( scalar == 0 ) {
//			GlbSysErrorMsg("C2DVector operator/ tried to divide by zero");
			return *this;
		}

		x /= scalar;
		y /= scalar;

		return (*this);
	}

	inline double Mod()
	{
		return sqrt(DotProd(*this, *this));
	}

	// array is 0-based (C/C++ convention)
    inline double &operator[](int i)
    {
		if( i == 0 )
			return x;
		if( i == 1 )
			return y;

//		GlbSysErrorMsg("Tried to access C2DVector element %i (out of range)", i);
		return y;
	}

	inline bool operator<(C2DVector &point)
	{
		return( x < point.x && y < point.y );
	}

	inline bool operator>(C2DVector &point)
	{
		return( x > point.x && y > point.y );
	}

	// tells if this point is left of, right of or on the vector
	// defined by p1-p2
	inline char IsLeftOrRight(C2DVector *p1, C2DVector *p2)
	{
		C2DVector v1(*p1, *this);
		C2DVector v2(*p1, *p2);
		double v2mod, normcross;

		v2mod = v2.Mod();

//		ASSERT(v2mod != 0.0);

		// compute cross product and normalize by v2mod
		// (so only sen(alpha) * v1mod remains, which is the distance
		// of point 'this' to line 'p1'-'p2')
		normcross = CrossProd(v1, v2) / v2mod;

		// side of 'this' is according to 'normcross'
		if(normcross < -C3D_TOL)
			return C2D_LEFT_OF;
		else if (normcross > C3D_TOL)
			return C2D_RIGHT_OF;
		else
			return C2D_ON_SEGMENT;
	}

	// check if given 'point' is inside, outside or on the circle
	// defined by three points p1, p2, p3
	inline char InCircle(C2DVector p1, C2DVector p2, C2DVector p3)
	{
		C2DVector center;
		double radius2, distx, disty, dist2;

		if( FindCircle3Points(p1, p2, p3, &center, &radius2) == false )
			return C2D_COLINEAR_POINTS;

		distx = (x - center.x);
		disty = (y - center.y);
		dist2 = distx*distx + disty*disty;

		if( dist2 < (radius2 - C3D_TOL * C3D_TOL / 4.0) )
			return C2D_IN_CIRCLE;
		else if( dist2 > (radius2 + C3D_TOL * C3D_TOL / 4.0) )
			return C2D_OUT_CIRCLE;
		else
			return C2D_ON_CIRCLE;
	}

	inline char InCircle(C2DVector *p1, C2DVector *p2, C2DVector *p3)
	{
		return InCircle(*p1, *p2, *p3);
	}
};

// parameters are passed by reference, but are not modified;
// therefore can be passed as const, in this way when passing
// a temporary object there are no issues with GCC compiler
// (otherwise DotProd(C2DVector(a), C2DVector(b)) would cause
// a compiler error
inline double DotProd(const C2DVector &vect1, const C2DVector &vect2)
{
	return vect1.x*vect2.x + vect1.y*vect2.y;
}

inline double DotProd(const C2DVector &vect1, const C2DVector_float &vect2)
{
	return vect1.x*vect2.x + vect1.y*vect2.y;
}

inline double DotProd(const C2DVector_float &vect1, const C2DVector &vect2)
{
	return vect1.x*vect2.x + vect1.y*vect2.y;
}

inline double Mod(const C2DVector &vect)
{
	return sqrt(DotProd(vect, vect));
}

inline double CrossProd(const C2DVector &vect1, const C2DVector &vect2)
{
	// cross-product always return a Z-oriented vector;
	// so specifying its module is enough
	return ( vect1.x*vect2.y - vect1.y*vect2.x );
}

inline C2DVector_float::C2DVector_float(C2DVector doubleVector)
{
	x = (C2DVector_float_type) doubleVector.x;
	y = (C2DVector_float_type) doubleVector.y;
}

inline void C2DVector_float::operator=(const C2DVector &point)
{
	x = (float)point.x;
	y = (float)point.y;
}

inline void C2DVector_float::operator+=(const C2DVector &point)
{
	x += (float)point.x;
	y += (float)point.y;
}


// A line in 2D can be represented as the intersection
// in parametric form, as a vector along the line direction and a point
// on the line (parameter origin on the line)
class C2DLine
{
public:
	C2DVector	m_vecDir;
	C2DVector	m_vecPoint;

	// default constructor; line starts at origin and is
	// directed towards positive Y-axis
	C2DLine()
	{
		m_vecDir.y = 1.0;
	}
	// line passing through two points, from v1 to v2
	C2DLine(C2DVector &v1, C2DVector &v2)
	{
		double mod;

		m_vecDir = v2 - v1;
		mod = Mod(m_vecDir);
		// if v1 and v2 are coincident (within tolerance) assert error
//		ASSERT(mod > C3D_TOL);
		// m_vecDir must always be normalized
		m_vecDir /= mod;

		m_vecPoint = v1;
	}

	bool IsPointOnLine(C2DVector &point)
	{
		C2DVector vec;
		double cross;

		// find vector from line origin to 'point'
		vec = point - m_vecPoint;

		// since Mod(m_vecDir) = 1.0, the result is Mod(v)*sin(alpha), that is,
		// the distance of the given point from the line
		cross = CrossProd(vec, m_vecDir);

		if(cross < C3D_TOL && cross > -C3D_TOL)
			return true;
		else
			return false;

	}

	// distance along the line of 'point' from line origin
	double DistFromOrigin(C2DVector &point)
	{
		C2DVector v;

		// find vector from line origin to 'point'
		v = point - m_vecPoint;

		return DotProd(v, m_vecDir);
	}

	// half plane in which point lies w.r.t. line origin
	// and direction
	char HalfPlane(C2DVector &point)
	{
		double dist;

		dist = DistFromOrigin(point);

		if( dist > C3D_TOL)
			return C2D_POS_HALFPLANE;
		else if ( dist < -C3D_TOL)
			return C2D_NEG_HALFPLANE;
		else
			return C2D_ON_ORIGIN;
	}
};

inline double CrossProd(const C2DLine &line1, const C2DLine &line2)
{
	return CrossProd(line1.m_vecDir, line2.m_vecDir);
}

inline double CrossProd(const C2DVector &vect, const C2DLine &line)
{
	return CrossProd(vect, line.m_vecDir);
}

// distance is positive if point is on the right half plane of line
inline double PointLineDist(C2DVector &point, C2DLine &line)
{
    C2DVector pointpoint(line.m_vecPoint, point);

	return ( CrossProd(pointpoint, line.m_vecDir) );
}

// return point at relative distance t0 on segment point_a-point_b,
// starting from point_a (t0 = 0)
// point = (1-t0)*point_a + t0*point_b
inline C2DVector PointOnSegment(C2DVector &point_a, C2DVector &point_b, double t0)
{
	C2DVector point;

//	ASSERT( t0>=0.0 && t0<=1.0 );

	point.x = (1-t0) * point_a.x + t0 * point_b.x;
	point.y = (1-t0) * point_a.y + t0 * point_b.y;

	return point;
}

// Compute segment - line intersection point
// If intersection exists, return C2D_OK and point 'point',
// otherwise return error condition
inline char IntersSegmentLine(C2DVector &point, C2DVector &point_a, C2DVector &point_b, C2DLine &line)
{
	double dista, distb, absdista, absdistb, len, t0;

	dista = PointLineDist(point_a, line);
	distb = PointLineDist(point_b, line);

	len = dista - distb;

	absdista = fabs(dista);
	absdistb = fabs(distb);

	// if both points are on the same half-plane side,
	// the segment does not intersect the line
	if( (dista > C3D_TOL && distb > C3D_TOL) || (dista < -C3D_TOL && distb < -C3D_TOL) ) {
		return C2D_SEGMENT_DONT_INTERSECT;
	}

	// If absdista == 0, point a is on the line
	if( absdista < C3D_TOL ) {
		// If both are on the line, return the first one
		// (either one would be ok) and signal it
		if( absdistb < C3D_TOL ) {
			point = point_a;
			return C2D_BOTH_POINTS_ON_LINE;
		}

		point = point_a;
		return C2D_POINT_ON_LINE;
	}

	// If distb == 0, point b is on the line
	if( absdistb < C3D_TOL ) {
		point = point_b;
		return C2D_POINT_ON_LINE;
	}

	// if len == 0, endpoints are coincident
	// (and cannot be both on the line, or should have
	// been catched before)
	if( len < C3D_TOL ) {
		// no intersection with line (this condition should
		// have been catched by the test on half plane containing points;
		// however, due to combining tolerances (C3D_TOL) when calculating len,
		// the condition could have been missed
		return C2D_SEGMENT_DONT_INTERSECT;
	}

	// compute relative intersection distance along segment
	t0 = dista / len;
	// determine intersection point at relative distance t0 from point_a
	point = PointOnSegment(point_a, point_b, t0);

	return C2D_OK;
}

inline bool FindCircle3Points(C2DVector p1, C2DVector p2, C2DVector p3, C2DVector *center, double *radius2 )
{
	double a, b, c, d, e, f, g, h, i;

	a = p2.x - p1.x;
	b = p2.y - p1.y;
	c = p3.x - p1.x;
	d = p3.y - p1.y;

	e = a * (p1.x + p2.x) + b * (p1.y + p2.y);
	f = c * (p1.x + p3.x) + d * (p1.y + p3.y);

	g = 2.0 * ( a * (p3.y - p2.y) - b * (p3.x - p2.x) );

	// if g is almost zero, points are almost collinear
	if( g < C3D_TOL )
		return false;

	// center of the circle
	center->x = (d*e - b*f) / g;
	center->y = (a*f - c*e) / g;

	h = p1.x - center->x;
	i = p1.y - center->y;
	// this is radius^2
	*radius2 = h*h + i*i;

	return true;
}

#endif //C2D_VECTOR2D_DEFS

