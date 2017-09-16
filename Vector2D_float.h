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


// Vector2D_float.h : definition of 3D vector class based on 'float' type
// Enrico Di Lorenzo 2013/02/01

#ifndef C3D_VECTOR2D_FLOAT_DEFS
#define C3D_VECTOR2D_FLOAT_DEFS

#include <math.h>

#include "GeoGlobal.h"

/////////////////////////////////////////////////////
// 2D vector class
//

#define C2DVector_float_type	float

// declared classes
class C2DVector;


// vector class specialized for 2d operations
class C2DVector_float
{
public:
	C2DVector_float_type			x;
	C2DVector_float_type			y;

	C2DVector_float() {}


	C2DVector_float(C2DVector_float_type xcoord, C2DVector_float_type ycoord) :
				x(xcoord), y(ycoord) {}

	inline C2DVector_float(C2DVector doubleVector);

	friend double DotProd(const C2DVector_float &vect1, const C2DVector_float &vect2);
    friend double Mod(const C2DVector_float &vect);
	friend double CrossProd(const C2DVector_float &vect1, const C2DVector_float &vect2);

	inline void pos(C2DVector_float_type xcoord, C2DVector_float_type ycoord)
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

	inline C2DVector_float operator-()
	{
		C2DVector_float newPoint;

		newPoint.x = -x;
		newPoint.y = -y;

		return(newPoint);
	}

	inline bool operator==(C2DVector_float &point)
	{
		return( (*this - point).Mod() < C3D_TOL );
	}

	inline void operator=(const C2DVector_float &point)
	{
		x = point.x;
		y = point.y;
	}

	// body defined in Vector2D.h
	inline void operator=(const C2DVector &point);

	// body defined in Vector2D.h
	inline void operator+=(const C2DVector &point);

	inline void operator+=(const C2DVector_float &point)
	{
		x += point.x;
		y += point.y;
	}

	inline C2DVector_float operator+(const C2DVector_float &point)
	{
		C2DVector_float newPoint;

		newPoint.x = x + point.x;
		newPoint.y = y + point.y;

		return newPoint;
	}

	inline void operator-=(C2DVector_float &point)
	{
		x -= point.x;
		y -= point.y;
	}

	inline C2DVector_float operator-(const C2DVector_float &point)
	{
		C2DVector_float newPoint;

		newPoint.x = x - point.x;
		newPoint.y = y - point.y;

		return newPoint;
	}

	inline C2DVector_float operator*(const double &scalar)
	{
		C2DVector_float newPoint;

		newPoint.x = x * scalar;
		newPoint.y = y * scalar;

		return newPoint;
	}

	inline C2DVector_float operator*(const C2DVector_float &point)
	{
		C2DVector_float newPoint;

		newPoint.x = x * point.x;
		newPoint.y = y * point.y;

		return newPoint;
	}

	inline void operator*=(const double &scalar)
	{
		x *= scalar;
		y *= scalar;
	}

	inline C2DVector_float operator/(const double &scalar)
	{
		C2DVector_float newPoint;

		if( scalar == 0 ) {
//			GlbSysErrorMsg("C2DVector operator/ tried to divide by zero");
			return *this;
		}

		newPoint.x = x / scalar;
		newPoint.y = y / scalar;

		return newPoint;
	}

	inline C2DVector_float &operator/=(const double &scalar)
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
    inline C2DVector_float_type &operator[](int i)
    {
		if( i == 0 )
			return x;
		if( i == 1 )
			return y;

//		GlbSysErrorMsg("Tried to access C2DVector element %i (out of range)", i);
		return y;
	}

	inline bool operator<(C2DVector_float &point)
	{
		return( x < point.x && y < point.y );
	}

	inline bool operator>(C2DVector_float &point)
	{
		return( x > point.x && y > point.y );
	}

};

// parameters are passed by reference, but are not modified;
// therefore can be passed as const, in this way when passing
// a temporary object there are no issues with GCC compiler
// (otherwise DotProd(C2DVector(a), C2DVector(b)) would cause
// a compiler error
inline double DotProd(const C2DVector_float &vect1, const C2DVector_float &vect2)
{
	return vect1.x*vect2.x + vect1.y*vect2.y;
}

inline double Mod(const C2DVector_float &vect)
{
	return sqrt(DotProd(vect, vect));
}

inline double CrossProd(const C2DVector_float &vect1, const C2DVector_float &vect2)
{
	// cross-product always return a Z-oriented vector;
	// so specifying its module is enough
	return ( vect1.x*vect2.y - vect1.y*vect2.x );
}

#endif //C3D_VECTOR2D_FLOAT_DEFS

