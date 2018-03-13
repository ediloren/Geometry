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


// GeoGlobal.h : definition of geometry-related global declarations
// Enrico Di Lorenzo 2003/05/05

#ifndef GEOGLOBAL_DEFS
#define GEOGLOBAL_DEFS

//#include "Global.h"
#include <new>

#define C3D_TOL		1E-12
#define C2D_TOL		1E-12

#define C3D_DIR_NONE	-1
#define C3D_DIR_X		0
#define C3D_DIR_Y		1
#define C3D_DIR_Z		2

#define C2D_DIR_X		0
#define C2D_DIR_Y		1

// minumum angle used in face quality triangulation
#define C3D_MINIMUM_ANGLE	20.0

// for multiplatform wxWidgets
#ifndef MS_VS
// for wxASSERT
#  include <wx/debug.h>
// when not using MS VisualC++
#  define ASSERT wxASSERT
#  define FALSE 0
#endif

#endif // GEOGLOBAL_DEFS

