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


// Vector3D_float.h : definition of 3D vector class based on 'float' type
// Enrico Di Lorenzo 2010/07/07

#ifndef C3D_VECTOR3D_FLOAT_DEFS
#define C3D_VECTOR3D_FLOAT_DEFS

#include <math.h>

#include "GeoGlobal.h"
//#include "Global.h"

/////////////////////////////////////////////////////
// 3D vector class
//

/*

// optimized pre-processor cross-product macro
#define C3D_CROSSPROD(result, vect1, vect2)			\
	result.x = vect1.y*vect2.z - vect1.z*vect2.y;	\
	result.y = -vect1.x*vect2.z + vect1.z*vect2.x;	\
	result.z = vect1.x*vect2.y - vect1.y*vect2.x;

// optimized pre-processor dot-product macro
#define C3D_DOTPROD(vect1, vect2)			\
	(vect1.x*vect2.x + vect1.y*vect2.y + vect1.z*vect2.z)

// optimized pre-processor module macro
#define C3D_MOD(vect1)			\
	sqrt(vect1.x*vect1.x + vect1.y*vect1.y + vect1.z*vect1.z)


// declared classes
class C3DVector;
class C3DLine;

*/

#define C3DVector_float_type	float

class C3DVector;

// vector class specialized for 3d operations
class C3DVector_float
{

public:
	C3DVector_float_type			x;
	C3DVector_float_type			y;
	C3DVector_float_type			z;

//	C3DVector() : x(0.0), y(0.0), z(0.0) {}
	C3DVector_float() {}

	C3DVector_float(C3DVector_float_type xcoord, C3DVector_float_type ycoord, C3DVector_float_type zcoord) :
				x(xcoord), y(ycoord), z(zcoord) {}

	inline C3DVector_float(C3DVector doubleVector);

	friend double DotProd(const C3DVector_float &vect1, const C3DVector_float &vect2);
	friend double Mod(const C3DVector_float &vect);
	friend C3DVector_float Cross(C3DVector_float &vect1, C3DVector_float &vect2);

	inline void pos(C3DVector_float_type xcoord, C3DVector_float_type ycoord, C3DVector_float_type zcoord)
	{
		x = xcoord;
		y = ycoord;
		z = zcoord;
	}

	inline void Invert()
	{
		x = -x;
		y = -y;
		z = -z;
	}

	inline void Normalize()
	{
		*this /= (float)Mod();
	}

	inline bool operator==(C3DVector_float &point)
	{
		return( (*this - point).Mod() < C3D_TOL );
	}

	inline bool operator!=(C3DVector_float &point)
	{
		return( !((*this - point).Mod() < C3D_TOL) );
	}

	inline void operator=(const C3DVector_float &point)
	{
		x = point.x;
		y = point.y;
		z = point.z;
	}

	// body defined in Vector3D.h
	inline void operator=(const C3DVector &point);

	inline void operator+=(C3DVector_float &point)
	{
		x += point.x;
		y += point.y;
		z += point.z;
	}

	inline void operator+=(C3DVector &point);

	inline C3DVector_float operator+(const C3DVector_float &point)
	{
		C3DVector_float newPoint;

		newPoint.x = x + point.x;
		newPoint.y = y + point.y;
		newPoint.z = z + point.z;

		return newPoint;
	}

	inline C3DVector operator+(C3DVector &point);

	inline void operator-=(C3DVector_float &point)
	{
		x -= point.x;
		y -= point.y;
		z -= point.z;
	}

	inline C3DVector_float operator-(const C3DVector_float &point)
	{
		C3DVector_float newPoint;

		newPoint.x = x - point.x;
		newPoint.y = y - point.y;
		newPoint.z = z - point.z;

		return newPoint;
	}

	inline C3DVector operator-(C3DVector &point);

	// TBC warning: the operator-() function has changed
	// because as it was it changed 'this' directly,
	// so negated permanently any vector!
	inline C3DVector_float operator-()
	{
		C3DVector_float tmp;

		tmp.x = -x;
		tmp.y = -y;
		tmp.z = -z;

		return tmp;
	}

	inline C3DVector_float operator*(const C3DVector_float_type &scalar)
	{
		C3DVector_float newPoint;

		newPoint.x = x * scalar;
		newPoint.y = y * scalar;
		newPoint.z = z * scalar;

		return newPoint;
	}

	inline C3DVector_float operator*(const C3DVector_float &point)
	{
		C3DVector_float newPoint;

		newPoint.x = x * point.x;
		newPoint.y = y * point.y;
		newPoint.z = z * point.z;

		return newPoint;
	}

	inline C3DVector_float operator/(const C3DVector_float_type &scalar)
	{
		C3DVector_float newPoint;

		if( scalar == 0 ) {
//			GlbSysErrorMsg("C3DVector operator/ tried to divide by zero");
			return *this;
		}

		newPoint.x = x / scalar;
		newPoint.y = y / scalar;
		newPoint.z = z / scalar;

		return newPoint;
	}

	inline C3DVector_float &operator/=(const C3DVector_float_type &scalar)
	{
		if( scalar == 0 ) {
//			GlbSysErrorMsg("C3DVector operator/ tried to divide by zero");
			return *this;
		}

		x /= scalar;
		y /= scalar;
		z /= scalar;

		return (*this);
	}

	inline double Mod()
	{
		return sqrt(DotProd(*this, *this));
	}

	// array is 0-based (C/C++ convention)
    inline C3DVector_float_type &operator[](int i)
    {
		if( i == 0 )
			return x;
		if( i == 1 )
			return y;
		if( i == 2 )
			return z;

//		GlbSysErrorMsg("Tried to access C3DVector element %i (out of range)", i);
		return z;
	}

	inline bool operator<(C3DVector_float &point)
	{
		return( x < point.x && y < point.y && z < point.z );
	}

	inline bool operator>(C3DVector_float &point)
	{
		return( x > point.x && y > point.y && z > point.z );
	}

};


inline double DotProd(const C3DVector_float &vect1, const C3DVector_float &vect2)
{
	return vect1.x*vect2.x + vect1.y*vect2.y + vect1.z*vect2.z;
}

inline double Mod(const C3DVector_float &vect)
{
	return sqrt(DotProd(vect, vect));
}

inline double Mod2(const C3DVector_float &vect)
{
	return DotProd(vect, vect);
}

inline C3DVector_float CrossProd(const C3DVector_float &vect1, const C3DVector_float &vect2)
{
	C3DVector_float tmp;

	tmp.x = vect1.y*vect2.z - vect1.z*vect2.y;
	tmp.y = -vect1.x*vect2.z + vect1.z*vect2.x;
	tmp.z = vect1.x*vect2.y - vect1.y*vect2.x;

	return tmp;
}

/*
// A line in 3D can be either represented as the intersection of two planes or,
// in parametric form, as a vector along the line direction and a point
// on the line (parameter origin on the line)
class C3DLine
{
public:
	C3DVector	m_vecDir;
	C3DVector	m_vecPoint;

	// default constructor; line starts at origin and is
	// directed towards positive Z-axis
	C3DLine()
	{
		m_vecDir.z = 1.0;
	}
};

inline C3DVector CrossProd(C3DLine &line1, C3DLine &line2)
{
	return CrossProd(line1.m_vecDir, line2.m_vecDir);
}
*/

#endif //C3D_VECTOR2D_DEFS

