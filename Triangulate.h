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


// Triangulate.h : Delaunay triangulation routines
// Enrico Di Lorenzo 2003/05/05

#ifndef C2D_TRIANGULATE_DEFS
#define C2D_TRIANGULATE_DEFS

// multiplatform condition for wxWidgets GCC vs. Visual Studio
#ifdef MS_VS
// disable warning C4786
// "'identifier' : identifier was truncated to 'number' characters in the debug information"
// caused by the use of the STL library
#pragma warning( disable : 4786 )
#endif

#include <algorithm>
#include <vector>
#include <list>
#include <utility>

#include "GeoGlobal.h"
#include "Vector2D.h"
#include "MathDefs.h"

#ifdef MS_VS
// for memory state and debug macros (e.g. _ASSERT), when using MS VisualC++
#  include <crtdbg.h>
#  define ASSERT _ASSERT
#endif

#define C2D_NOT_CONSTRAINED			0
#define C2D_CONSTRAINED				1
#define C2D_BOUNDARY_CONSTRAINED	2

#define C2D_NOBOUND_NOSMALLEDGE		-1

#define C2D_UNINIT_INDEX			-1
#define C2D_DUMMY_INDEX				-2

#define C2D_VERTEX_NOT_FOUND		-1


/////////////////////////////////////////////////////
//
//
//

using namespace std;

// REMARK:
// no use to declare operators as friends: due to MS acknowledged bug,
// when using 'using namespace std' directive, friends don't work;
// in case bug will be fixed, can protect class member variables
// and use friend access for the operator

// prototypes
class C2DEdge;
class C2DTriplet;
class C2DTriangle;
class C2DEdgeTriangle;
class C2DVertex;
class C2DTriangulate;


// triplet of array indexes, identifying triangle
// vertexes through vertex array
class C2DTriplet
{
	void Order()
	{
		long i, j, k;

		i = vertex[0];
		j = vertex[1];
		k = vertex[2];

		// insert first element
		vertex[0] = i;

		// insert second element
		if(j > vertex[0])
			vertex[1] = j;
		else {
			vertex[0] = j;
			vertex[1] = i;
		}

		// insert second element
		if(k > vertex[1])
			vertex[2] = k;
		else {
			if( k > vertex[0] ) {
				vertex[2] = vertex[1];
				vertex[1] = k;
			}
			else {
				vertex[2] = vertex[1];
				vertex[1] = vertex[0];
				vertex[0] = k;
			}
		}
	}

public:
	long vertex[3];

	C2DTriplet()
	{
	}

	C2DTriplet(long i, long j, long k)
	{
		Set(i, j, k);
	}

	void Set(long i, long j, long k)
	{
		vertex[0] = i;
		vertex[1] = j;
		vertex[2] = k;

		Order();
	}

	void operator=(const C2DTriplet &triplet)
	{
		vertex[0] = triplet.vertex[0];
		vertex[1] = triplet.vertex[1];
		vertex[2] = triplet.vertex[2];
	}

	// triplets are ordered by default
	bool operator==(C2DTriplet &triplet)
	{
		if( vertex[0] == triplet.vertex[0] &&
			vertex[1] == triplet.vertex[1] &&
			vertex[2] == triplet.vertex[2])
			return true;
		else
			return false;
	}

	// no use to declare it friend: due to MS acknowledged bug,
	// when using 'using namespace std' directive, friends don't work;
	// in cas bug will be fixed, can protect class member variables
	// and use friend access for the operator
	friend bool operator<(const C2DTriplet &t1, const C2DTriplet &t2);
};

class C2DTriangle
{
	// A C2DTriangle is identified by its three 'vertexes' indexes
	long vertex[3];
public:
	C2DTriangle()
	{
		vertex[0] = C2D_UNINIT_INDEX;
		vertex[1] = C2D_UNINIT_INDEX;
		vertex[2] = C2D_UNINIT_INDEX;
	}

	C2DTriangle(long v1, long v2, long v3)
	{
		vertex[0] = v1;
		vertex[1] = v2;
		vertex[2] = v3;
	}

	// return the vertex index associated to given edge #
	long GetVertex(unsigned char i)
	{
		ASSERT(i<=2);
		return(vertex[i]);
	}

	void SetVertexes(long v1, long v2, long v3)
	{
		vertex[0] = v1;
		vertex[1] = v2;
		vertex[2] = v3;
	}

	bool IsGivenTriangle(long i, long j, long k)
	{
		C2DTriplet givenTriplet(i, j, k);
		C2DTriplet localTriplet(vertex[0], vertex[1], vertex[2]);

		if( localTriplet == givenTriplet )
			return true;
		else
			return false;
	}

	bool operator==(C2DTriangle &tri)
	{
		C2DTriplet givenTriplet(tri.vertex[0], tri.vertex[1], tri.vertex[2]);
		C2DTriplet localTriplet(vertex[0], vertex[1], vertex[2]);

		if( givenTriplet == localTriplet )
			return true;
		else
			return false;
	}

};


class C2DEdgeTriangle
{
	// A C2DEdgeTriangle is identified by its three edges
	C2DEdge *edges[3];
public:

	C2DEdgeTriangle()
	{
		edges[0] = NULL;
		edges[1] = NULL;
		edges[2] = NULL;
	}

	C2DEdge *GetEdge(unsigned char i)
	{
		ASSERT(i<=2);
		return(edges[i]);
	}

	void SetEdges(C2DEdge *edge1, C2DEdge *edge2, C2DEdge *edge3)
	{
		edges[0] = edge1;
		edges[1] = edge2;
		edges[2] = edge3;
	}

	long GetVertex(unsigned char i);
};


class C2DEdge {
	// identificator (e.g. to mark visited edges)
	long id;
	// associated vertex (index in 'vertexes' array)
	long vertex;
	// constrained status
	char constrainStatus;
	// pointers to other edges
	C2DEdge *onext;
	C2DEdge *oprev;
	C2DEdge *mate;
	// pointer to triangle in small triangles list
	list<C2DEdgeTriangle>::iterator smallTriangle;
	// boolean to signal if there is small triangle associated
	bool smallTriangleLink;

	C2DEdge(C2DEdge *mymate, long v) : id(0), constrainStatus(false), smallTriangleLink(false)
	{
		onext = this;
		oprev = this;
		mate = mymate;
		vertex = v;
	}

public:

	friend char JoinEdges(C2DEdge *edgea, C2DEdge *edgeb);
	friend double CrossProd(C2DEdge *edgea, C2DEdge *edgeb);
	friend C2DEdge *JoinEdgesRight(C2DEdge *edgea, C2DEdge *edgeb);
	friend C2DEdge *JoinEdgesLeft(C2DEdge *edgea, C2DEdge *edgeb);
	friend void DeleteEdge(C2DEdge *edge);
	friend C2DEdge *CreateEdge(long v1, long v2);
	friend C2DVector *Vertex(C2DEdge *edge);
	friend C2DVector *Vertex(long vertex);
	friend C2DVector *MateVertex(C2DEdge *edge);

	C2DEdge(long v1, long v2) : id(0), constrainStatus(false), smallTriangleLink(false)
	{
		onext = this;
		oprev = this;
		vertex = v1;
		mate = new C2DEdge(this, v2);
	}

	void SetMate(C2DEdge *edge)
	{
		mate = edge;
	}

	void SetVertex(long v)
	{
		vertex = v;
	}

	// next edge moving cw around first vertex of current edge (origin next)
	C2DEdge *Oprev()
	{
		return oprev;
	}

	void SetOprev(C2DEdge *next)
	{
		oprev = next;
	}

	// next edge moving ccw around first vertex of current edge (origin prev)
	C2DEdge *Onext()
	{
		return onext;
	}

	void SetOnext(C2DEdge *next)
	{
		onext = next;
	}

	// next edge moving ccw along the contour (left next)
	C2DEdge *Lnext()
	{
		// to move ccw along the contour, must move cw around Mate() origin
		return Mate()->Oprev();
	}

	// prev edge moving ccw along the contour (left next)
	C2DEdge *Lprev()
	{
		return Onext()->Mate();
	}

	C2DEdge *Mate()
	{
		return mate;
	}

	long Vertex()
	{
		return vertex;
	}

	long MateVertex()
	{
		ASSERT(mate != NULL);

		return mate->Vertex();
	}

	long GetId()
	{
		return id;
	}

	void SetId(long newid)
	{
		id = newid;
	}

	void Constrain()
	{
		// must constrain both edge and its mate
		constrainStatus = C2D_CONSTRAINED;
		Mate()->constrainStatus = C2D_CONSTRAINED;
	}

	void BoundaryConstrain()
	{
		// must constrain both edge and its mate
		constrainStatus = C2D_BOUNDARY_CONSTRAINED;
		Mate()->constrainStatus = C2D_BOUNDARY_CONSTRAINED;
	}

	bool IsConstrained()
	{
		return(constrainStatus != C2D_NOT_CONSTRAINED);
	}

	bool IsHardConstrained()
	{
		return(constrainStatus == C2D_CONSTRAINED);
	}

	// Add an edge to the double-linked list of edges radial
	// to a vertex. Note that order matters: the new edge will
	// be the next ccw edge to the existing edge (this)
	void LinkCCW(C2DEdge *newedge)
	{
		C2DEdge *next;

		// save pointer to old next edge
		next = Onext();

		// link new edge
		newedge->SetOprev(this);
		newedge->SetOnext(next);

		// adjust this edge and this edge's onext references
		SetOnext(newedge);
		next->SetOprev(newedge);
	}

	// Add an edge to the double-linked list of edges radial
	// to a vertex. Note that order matters: the new edge will
	// be the next cw edge to the existing edge (this)
	void LinkCW(C2DEdge *newedge)
	{
		C2DEdge *prev;

		// save pointer to old prev edge
		prev = Oprev();

		// link new edge
		newedge->SetOnext(this);
		newedge->SetOprev(prev);

		// adjust this edge and this edge's onext references
		SetOprev(newedge);
		prev->SetOnext(newedge);
	}

	// add a vertex along 'this' edge
	//
	//           this
	//     --------------->
	//     <---------------
	//          Mate()
	//
	//  becomes
	//
	//       this   newedge2
	//     ------->O------>
	//     <-------O<------
	//     newedge1  old Mate()
	//
	// warning: this is only a topology operation;
	// no check is made that vertex is on the given edge!
	void Split(long vertex)
	{
		C2DEdge *newedge1, *newedge2;

		// create new edges
		newedge1 = new C2DEdge(this, vertex);
		newedge2 = new C2DEdge(Mate(), vertex);

		// make new mates
		Mate()->SetMate(newedge2);
		SetMate(newedge1);

		// and adjust link pointers
		newedge1->SetOprev(newedge2);
		newedge1->SetOnext(newedge2);
		newedge2->SetOprev(newedge1);
		newedge2->SetOnext(newedge1);

		// then inheredit constrain status information
		newedge1->constrainStatus = constrainStatus;
		newedge2->constrainStatus = constrainStatus;
	}

	void SetSmallTriangle(list<C2DEdgeTriangle>::iterator it)
	{
		smallTriangle = it;
	}

	list<C2DEdgeTriangle>::iterator GetSmallTriangle()
	{
		return smallTriangle;
	}

	bool GetSmallTriangleFlag()
	{
		return smallTriangleLink;
	}

	void SetSmallTriangleFlag(bool flag)
	{
		smallTriangleLink = flag;
	}

};


class C2DVertex {
public:
	// default constructor
	C2DVertex() : sortedVtoV(C2D_UNINIT_INDEX), edge(NULL) {}

	C2DVertex(long sVtoV, C2DEdge *e) : sortedVtoV(sVtoV), edge(e) {}

	// reference to user input vertexes list from sorted vertex
	long sortedVtoV;
	// reference to one of the radial edges
	C2DEdge *edge;
};


class C2DTriangulate
{
	// add a new vertex in the triangulation
	// Remark: this vertex must be unique!
	// No uniquneness test is done.
	long AddVertex(C2DVector &point)
	{
		long sortVertexesSize, vToSortedVsize;

		sortVertexesSize = sortVertexes.size();
		vToSortedVsize = vToSortedV.size();

		// add new vertex

		// in vertexes array
		vertexes.push_back(point);
		// in sorted vertexes array (which becomes unsorted!)
		sortVertexes.push_back(C2DVertex(vToSortedVsize, NULL));
		// in references array
		vToSortedV.push_back(sortVertexesSize);

		// return last element (array is 0-based)
		return (sortVertexesSize);
	}

public:
	C2DEdge *CreateEdge(long v1, long v2)
	{
		C2DEdge *edge;

		// create edge
		edge = new C2DEdge(v1, v2);

		// and add reference from vertex to one of the radial
		// edges, if not already assigned

		if( GetVertexEdge(v1) == NULL )
			SetVertexEdge(v1, edge);

		if( GetVertexEdge(v2) == NULL )
			SetVertexEdge(v2, edge->Mate());

		return edge;
	}

	// like CreateEdge, but uses an existing edge pointer
	void SetEdge(C2DEdge *edge, long v1, long v2)
	{
		// set edge vertexes
		edge->SetVertex(v1);
		edge->Mate()->SetVertex(v2);

		// and add reference from vertex to one of the radial
		// edges, if not already assigned

		if( GetVertexEdge(v1) == NULL )
			SetVertexEdge(v1, edge);

		if( GetVertexEdge(v2) == NULL )
			SetVertexEdge(v2, edge->Mate());
	}

	void DetachHalfEdge(C2DEdge *edge)
	{
		// detach from list of radials
		edge->Onext()->SetOprev(edge->Oprev());
		edge->Oprev()->SetOnext(edge->Onext());

		// if the vertex referenced this edge of all its radials
		if( GetVertexEdge(edge->Vertex()) == edge ) {
			// and this is not the last edge in the list
			if( edge->Onext() != edge ) {
				// change vertex reference
				SetVertexEdge(edge->Vertex(), edge->Onext());
			}
			else {
				SetVertexEdge(edge->Vertex(), NULL);
			}
		}
	}

	void DetachEdge(C2DEdge *edge)
	{
		DetachHalfEdge(edge);
		DetachHalfEdge(edge->Mate());
	}

	void DeleteEdge(C2DEdge *edge)
	{
		C2DEdge *mate;

		mate = edge->Mate();

		DetachHalfEdge(edge);
		DetachHalfEdge(mate);

		delete edge;
		delete mate;
	}

	void DeleteHalfEdge(C2DEdge *edge)
	{
		DetachHalfEdge(edge);

		// NULL pointer to mate
		if(edge->Mate() != NULL)
			edge->Mate()->SetMate(NULL);

		delete edge;
	}

	C2DVector *Vertex(C2DEdge *edge)
	{
		return &vertexes[sortVertexes[edge->Vertex()].sortedVtoV];
	}

	C2DVector *Vertex(long vertex)
	{
		return &vertexes[sortVertexes[vertex].sortedVtoV];
	}

	C2DVector *MateVertex(C2DEdge *edge)
	{
		return &vertexes[sortVertexes[edge->MateVertex()].sortedVtoV];
	}

	C2DEdge *GetVertexEdge(long vertex)
	{
		return sortVertexes[vertex].edge;
	}

	void SetVertexEdge(long vertex, C2DEdge *edge)
	{
		sortVertexes[vertex].edge = edge;
	}

	long GetVIndexFromSortV(long index)
	{
		return sortVertexes[index].sortedVtoV;
	}

	long GetSortVIndexFromV(long index)
	{
		long reference;

		reference = vToSortedV[index];
		// if reference points to a vertex which has
		// been deleted because coincident with some other
		if(reference < 0)
			return vToSortedV[-reference+1];
		else
			return reference;
	}

	// Cross-product of the vectors defined by the two edges
	//             _
	//     edgeb   /|
	//           /
    //          O------>
	//           edgea
	//
	inline double CrossProd(C2DEdge *edgea, C2DEdge *edgeb)
	{
		C2DVector v1, v2;

		v1 = *MateVertex(edgea) - *Vertex(edgea);
		v2 = *MateVertex(edgeb) - *Vertex(edgeb);

		return ::CrossProd(v1, v2);
	}

	// Cross-product of the vectors defined by the one edge and one line
	//             _
	//      line   /|
	//           /
    //          O------>
	//           edge
	//
	inline double CrossProd(C2DEdge *edge, C2DLine &line)
	{
		C2DLine l1(*Vertex(edge), *MateVertex(edge));

		return ::CrossProd(l1, line);
	}

	// join two consecutive edges
	//
	//             _
	//     edgea   /|\    edgeb
	//           /    _\|
    //          O------>O
	//           joinedge
	//
	C2DEdge *JoinEdges(C2DEdge *edgea, C2DEdge *edgeb, char *side)
	{
		C2DVector *va, *vb, *vc;

		va = Vertex(edgea);
		vb = Vertex(edgeb);
		vc = MateVertex(edgeb);

		// verify if apex point of the two edges a and b is left or right
		// of segment joining extremes ( /\ or \/ case)
		//       vb           va    vc
		//      /  \     or     \  /
		//     va   vc           vb

		*side = vb->IsLeftOrRight(va, vc);

		if(*side == C2D_LEFT_OF) {
			return JoinEdgesLeft(edgea, edgeb->Mate());
		}
		else if(*side == C2D_RIGHT_OF) {
			return JoinEdgesRight(edgea, edgeb->Mate());
		}

		// if side == C2D_ON_SEGMENT, do nothing

		return NULL;
	}

	// join two edges; 'right' refers to the side of the
	// new edge to be inserted w.r.t. the two existing edges
	//             _
	//        _    /|    _
	//       |\   /      /|
	//   edgea \ /      /  edgeb
    //          O----->O
	//          joinedge
	//
	C2DEdge *JoinEdgesRight(C2DEdge *edgea, C2DEdge *edgeb)
	{
		C2DEdge *joinedge;

		// create new joining edge
		joinedge = CreateEdge(edgea->Vertex(), edgeb->Vertex());
		// and link it radially
		edgea->LinkCCW(joinedge);
		edgeb->LinkCW(joinedge->Mate());

		return joinedge;
	}

	// join two edges; 'left' refers to the side of the
	// new edge to be inserted w.r.t. the two existing edges
	//
	//          joinedge
	//          O----->O
	//   edgea / \      \  edgeb
	//       |/   \      \|
    //        -    \|    -
	//             -
	C2DEdge *JoinEdgesLeft(C2DEdge *edgea, C2DEdge *edgeb)
	{
		C2DEdge *joinedge;

		// create new joining edge
		joinedge = CreateEdge(edgea->Vertex(), edgeb->Vertex());
		// and link it radially
		edgea->LinkCW(joinedge);
		edgeb->LinkCCW(joinedge->Mate());

		return joinedge;
	}

	// Join two edges using the given 'joinedge' edge.
	// 'left' refers to the side of the
	// new edge to be inserted w.r.t. the two existing edges
	//
	//          joinedge
	//          O----->O
	//   edgea / \      \  edgeb
	//       |/   \      \|
    //        -    \|    -
	//             -
	void JoinEdgesLeft(C2DEdge *joinedge, C2DEdge *edgea, C2DEdge *edgeb)
	{
		// create new joining edge
		SetEdge(joinedge, edgea->Vertex(), edgeb->Vertex());
		// and link it radially
		edgea->LinkCW(joinedge);
		edgeb->LinkCCW(joinedge->Mate());

	}

	// split 'edge' on vertex 'point'
	//
	//           edge
	//     --------------->
	//     <---------------
	//        edge->Mate()
	//
	//  becomes
	//
	//       edge
	//     ------->O------>
	//     <-------O<------
	//   edge->Mate()
	//
	// the function returns the new 'edge', that is, the new edge
	// starting from the same point of the old 'edge' and going
	// in the same direction
	//
	// warning: no check is made that vertex is on the given edge!
	C2DEdge * Split(C2DEdge *edge, C2DVector &point)
	{
		long index;

		// insert new vertex in triangulation
		index = AddVertex(point);

		// topologically split edge
		edge->Split(index);

		// must be NULL: verex has just been created
		ASSERT(GetVertexEdge(index) == NULL);

		// initialize vertexedge reference
		SetVertexEdge(index, edge->Mate());

		return edge;
	}

	// returns true if on the left side of the given edge
	// there is a triangle
	bool HasEdgeLeftTri(C2DEdge *edge, C2DEdge **edge2, C2DEdge **edge3)
	{
		long apex1, apex2;

		// get the two adjacent edges (for a CCW triangle)
		//
		//         *edge3
		//       O<---_O
		//       |    /|
		//  edge |   /
		//       |  / *edge2
		//       | /
		//      \//
		//       O
		//
		*edge2 = edge->Mate()->Oprev();
		*edge3 = edge->Onext()->Mate();

		// apex candidates
		apex1 = (*edge2)->Mate()->Vertex();
		apex2 = (*edge3)->Vertex();

		// if the two adjacent edges share the apex vertex,
		// the three edges form a triangle
		if( apex1 == apex2 ) {
			// note however that the triangle could also be on its right
			// in some particular cases, so must test apex point for
			// leftness
			if( Vertex(apex1)->IsLeftOrRight(Vertex(edge), MateVertex(edge)) == C2D_LEFT_OF )
				return true;
		}

		return false;
	}

	// returns true if on the right side of the given edge
	// there is a triangle
	bool HasEdgeRightTri(C2DEdge *edge, C2DEdge **edge2, C2DEdge **edge3)
	{
		return HasEdgeLeftTri(edge->Mate(), edge2, edge3);
	}

	// returns true if given triangle is present in the mesh
	bool IsMeshTriangle(C2DEdgeTriangle &triangle)
	{
		if( (triangle.GetEdge(0)->Mate()->Oprev() == triangle.GetEdge(1)) &&
			(triangle.GetEdge(0)->Onext()->Mate() == triangle.GetEdge(2)) )
			return true;
		else
			return false;
	}

	// return the number of generated triangles
	long GetTriangleNum()
	{
		return triangles.size();
	}

	// return a triangle in user coordinates from the generated triangle array
	C2DTriangle GetTriangle(long index)
	{
		C2DTriangle triangle;
		long v1, v2, v3;

		v1 = GetVIndexFromSortV(triangles[index].GetVertex(0));
		v2 = GetVIndexFromSortV(triangles[index].GetVertex(1));
		v3 = GetVIndexFromSortV(triangles[index].GetVertex(2));

		triangle.SetVertexes(v1, v2, v3);

		return triangle;
	}

	// return triangle coordinates from the generated triangle array
	void GetTriangleCoords(long index, C2DVector *v1, C2DVector *v2, C2DVector *v3)
	{
		C2DTriangle triangle;

		*v1 = *Vertex(triangles[index].GetVertex(0));
		*v2 = *Vertex(triangles[index].GetVertex(1));
		*v3 = *Vertex(triangles[index].GetVertex(2));
	}

	// remove triangle from smalltri list, if in list
	void RemoveSmallTriangle(C2DEdgeTriangle &triangle)
	{
		C2DEdge *edge1, *edge2, *edge3;

		// get triangle edges
		edge1 = triangle.GetEdge(0);
		edge2 = triangle.GetEdge(1);
		edge3 = triangle.GetEdge(2);

		// the edges, belonging to the same triangle,
		// cannot have different 'small triangle' flags
		// (i.e., they must be all associated to a small
		// trangle, or not associated at all)
		ASSERT( edge1->GetSmallTriangleFlag() == edge2->GetSmallTriangleFlag() );
		ASSERT( edge1->GetSmallTriangleFlag() == edge3->GetSmallTriangleFlag() );
		ASSERT( edge2->GetSmallTriangleFlag() == edge3->GetSmallTriangleFlag() );

		if( edge1->GetSmallTriangleFlag() ) {
			// clear flags
			edge1->SetSmallTriangleFlag(false);
			edge2->SetSmallTriangleFlag(false);
			edge3->SetSmallTriangleFlag(false);
			// the edges, belonging to the same triangle in 'smalltri' list
			// cannot have different 'small triangle' pointers
			ASSERT( edge1->GetSmallTriangle() == edge2->GetSmallTriangle());
			ASSERT( edge1->GetSmallTriangle() == edge3->GetSmallTriangle() );
			ASSERT( edge2->GetSmallTriangle() == edge3->GetSmallTriangle() );
			// erase from list
			smalltri.erase(edge1->GetSmallTriangle());
		}
	}

	// remove triangle (identified from edge) from smalltri list, if in list
	void RemoveSmallTriangle(C2DEdge *edge)
	{
		C2DEdgeTriangle tri;
		C2DEdge *edge2, *edge3;

		edge2 = edge->Lnext();
		edge3 = edge2->Lnext();

		tri.SetEdges(edge, edge2, edge3);
		RemoveSmallTriangle(tri);
	}

	// add triangle to smalltri list
	void AddSmallTriangle(C2DEdgeTriangle &triangle)
	{
		C2DEdge *edge1, *edge2, *edge3;

		// get triangle edges
		edge1 = triangle.GetEdge(0);
		edge2 = triangle.GetEdge(1);
		edge3 = triangle.GetEdge(2);

		// set flags
		edge1->SetSmallTriangleFlag(true);
		edge2->SetSmallTriangleFlag(true);
		edge3->SetSmallTriangleFlag(true);

		// insert triangle in small triangles set
		smalltri.push_front(triangle);

		// set edge reference to small triangle
		edge1->SetSmallTriangle(smalltri.begin());
		edge2->SetSmallTriangle(smalltri.begin());
		edge3->SetSmallTriangle(smalltri.begin());
	}

	// insert triangle (identified from edge) in the small triangles list, if small-angled
	void SmallTriangleInsert(C2DEdge *edge)
	{
		C2DEdgeTriangle tri;
		C2DEdge *edge2, *edge3;

		edge2 = edge->Lnext();
		edge3 = edge2->Lnext();

		tri.SetEdges(edge, edge2, edge3);
		SmallTriangleInsert(tri);
	}

	// get triangle's angle cosines
	//
	//             GetEdge(2)
	// cosAngle[2] ________  cosAngle[1]
	//             \      /
	//   GetEdge(0) \    /  GetEdge(1)
	//               \  /
	//                \/
	//            cosAngle[0]
	//
	void ComputeCosAngles(C2DEdgeTriangle &triangle, double cosAngle[3])
	{
		C2DVector *v1, *v2, *v3, l1, l2, l3;
		double modl1, modl2, modl3;
		C2DEdge *edge1, *edge2, *edge3;

		// get triangle edges
		edge1 = triangle.GetEdge(0);
		edge2 = triangle.GetEdge(1);
		edge3 = triangle.GetEdge(2);
		// get vertexes positions
		v1 = Vertex(edge1);
		v2 = Vertex(edge2);
		v3 = Vertex(edge3);
		// find edge vectors
		l1 = *v2-*v1;
		l2 = *v3-*v2;
		l3 = *v1-*v3;
		// and their length
		modl1 = Mod(l1);
		modl2 = Mod(l2);
		modl3 = Mod(l3);
		// compute triangle's cosine of the angles
		cosAngle[0] = DotProd(l1,-l2) / (modl1 * modl2);
		cosAngle[1] = DotProd(l2,-l3) / (modl2 * modl3);
		cosAngle[2] = DotProd(l3,-l1) / (modl3 * modl1);
	}

	// verify the non-cycle assumption on given triangle
	bool NonCycleTriangle(C2DEdgeTriangle &triangle, C2DEdge *edges[3], unsigned char edgeLenOrder[3])
	{
		double cosAngle[3], smallestAngle;
		unsigned char angleOrder[3];
		bool boundaryEdges[3];
		C2DVector midpointP;
		C2DEdge *longestEdge, *mediumEdge;

		ComputeCosAngles(triangle, cosAngle);
		GetDoubleOrder(cosAngle, angleOrder);

		// if smallest angle is > 22 degrees and < 30 degrees,
		// cycle condition may happen
		// (note that smallest angle has largest cosine)
		smallestAngle = cosAngle[angleOrder[2]];
		if( smallestAngle > COS30DEG && smallestAngle < COS22DEG) {
			// if largest angle is <= 120 degrees, cycle condition may happen
			if( cosAngle[angleOrder[0]] >= COS120DEG ) {
				// if the small angled vertex belongs to the boundary, cycle condition may happen
				if( IsBoundaryVertex(triangle.GetEdge(angleOrder[2])->MateVertex()) ) {

					// if the circumcircle of triangle PBC (where ABC is given
					// triangle, AB >= BC >= AC, and P is midpoint of AB)
					// does not contain the third vertex Q of neighbor
					// triangle BCQ, then we have a cycle condition

					// if BCQ exists (egde BC is not on the boundary)
					FindBoundaryEdges(edges, boundaryEdges);
					if( boundaryEdges[edgeLenOrder[1]] == false ) {

						longestEdge = edges[edgeLenOrder[2]];
						midpointP = ( *Vertex(longestEdge) + *MateVertex(longestEdge) ) / 2.0;
						mediumEdge = edges[edgeLenOrder[1]];

						// if vertex Q is not contained in circumcircle
						if( MateVertex(mediumEdge->Mate()->Onext())->InCircle(*Vertex(mediumEdge), *MateVertex(mediumEdge), midpointP) != C2D_IN_CIRCLE) {
							// cycle condition!
							return false;
						}
					}
				}
			}
		}

		return true;
	}

protected:
	// The edge is swapped rotating it ccw inside the four-sided encompassing polygon;
	// the routine returns the swapped edge
	C2DEdge *SwapEdge(C2DEdge *edge)
	{
		C2DEdge *onext, *oprev;

		// store pointers (will be lost on deletion)
		onext = edge->Onext();
		oprev = edge->Oprev();

		// detach old edge (better not delete it: in CDT operations,
		// there could be some pointer left to this edge, so it's
		// ok if swaps but not if it is deleted
		DetachEdge(edge);

		// and create the swapped one
		JoinEdgesLeft(edge, oprev->Mate(), onext->Mate());

		return edge;
	}

	// find the length order of the edges of the given triangle,
	// from smallest to biggest; in the meanwhile, fill also
	// 'edges' and 'len' arrays
	void GetLenOrder(C2DEdgeTriangle &triangle, C2DEdge *edges[3], double len[3], unsigned char lenorder[3])
	{
		edges[0] = triangle.GetEdge(0);
		edges[1] = triangle.GetEdge(1);
		edges[2] = triangle.GetEdge(2);

		len[0] = Length(edges[0]);
		len[1] = Length(edges[1]);
		len[2] = Length(edges[2]);

		GetDoubleOrder(len, lenorder);
	}

	// find the value order of given double array of dim=3,
	// from smallest to biggest
	void GetDoubleOrder(double value[3], unsigned char valOrder[3])
	{
		if(value[0] <= value[1] && value[0] <= value[2]) {
			if(value[1] <= value[2]) {
				valOrder[0] = 0;
				valOrder[1] = 1;
				valOrder[2] = 2;
			}
			else {
				valOrder[0] = 0;
				valOrder[1] = 2;
				valOrder[2] = 1;
			}
		}

		if(value[1] <= value[0] && value[1] <= value[2]) {
			if(value[0] <= value[2]) {
				valOrder[0] = 1;
				valOrder[1] = 0;
				valOrder[2] = 2;
			}
			else {
				valOrder[0] = 1;
				valOrder[1] = 2;
				valOrder[2] = 0;
			}
		}

		if(value[2] <= value[0] && value[2] <= value[1]) {
			if(value[0] <= value[1]) {
				valOrder[0] = 2;
				valOrder[1] = 0;
				valOrder[2] = 1;
			}
			else {
				valOrder[0] = 2;
				valOrder[1] = 1;
				valOrder[2] = 0;
			}
		}
	}

	// check current triangle's edges to determine if they are
	// boundary edges
	void FindBoundaryEdges(C2DEdge *edges[3], bool boundaryedges[3])
	{
		unsigned char i;
		C2DEdge *tmpedge1, *tmpedge2;

		for(i= 0; i<3; i++) {
			// if this is a boundary edge
			if(HasEdgeRightTri(edges[i], &tmpedge1, &tmpedge2) == false) {
				boundaryedges[i] = true;
			}
			else {
				boundaryedges[i] = false;
			}
		}
	}

	// check current triangle's edges to determine if they are
	// boundary edges
	bool IsBoundaryVertex(long vertex)
	{
		C2DEdge *tmpedge1, *tmpedge2, *walkingedge, *startedge;

		startedge = GetVertexEdge(vertex);
		walkingedge = startedge;
		do {
			// if this is a boundary edge, 'vertex' is a boundary vertex
			if(HasEdgeRightTri(walkingedge, &tmpedge1, &tmpedge2) == false)
				return true;
			walkingedge = walkingedge->Onext();
		} while(walkingedge != startedge);

		return false;
	}

	double Length(C2DEdge *edge)
	{
		return Mod(*MateVertex(edge) - *Vertex(edge));
	}

protected:
	vector<C2DVertex> sortVertexes;
	vector<long> vToSortedV;
	long index;
	bool meshIsPresent;
	double cosEpsilon;
	list<C2DEdgeTriangle> smalltri;
	long m_lMark;

	void DivideAndMerge(long left, long right, C2DEdge **right_cw, C2DEdge **left_ccw);
	C2DEdge *Merge(C2DEdge *right_cw_left, C2DEdge *left_ccw_right);
	void FindLowestTangent(C2DEdge *right_cw_left, C2DEdge *left_ccw_right,
						   C2DEdge **lowest_left, C2DEdge **lowest_right);
	void ConstrainBoundary(C2DEdge *boundaryedge);
	long NewIndex() { return ++index; }
	void Retriangulate(C2DEdge *first);
	void DelaunayEdge(C2DEdge *first);
	void RefineTriangle(C2DEdgeTriangle &triangle);
	void SmallTriangleInsert(C2DEdgeTriangle &triangle);
	void LEPPPointInsert(C2DEdgeTriangle triangle);
	void BoundaryPointInsert(C2DEdgeTriangle &triangle, char edgeindex);
	void TerminalPointInsert(C2DEdgeTriangle &triangle, char edgeindex, C2DEdgeTriangle &triangle2);
	void LEPPDelaunayEdge(C2DEdge *edge);
	void DeleteEdges(C2DEdge *edge);
	void SortVertexes();
	void InitSortVertexes();
	long FindVertex(C2DVector point);

public:
	vector<C2DVector> vertexes;
	vector<C2DEdgeTriangle> triangles;

	C2DTriangulate() : index(0), meshIsPresent(false) {}
	~C2DTriangulate();
	void Triangulate();
	bool InsertConstrEdge(C2DVector p1, C2DVector p2);
	bool InsertConstrEdge(long a, long b);
	void GenerateTriangles();
	void Refine(double epsilon);
	bool DeleteCWEdges(C2DVector p1, C2DVector p2);
	bool DeleteCWEdges(long a, long b);
	void ClearMesh();
};


#endif //C2D_TRIANGULATE_DEFS

