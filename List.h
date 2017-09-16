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


// List.h : definition of list template - related classes
// Enrico Di Lorenzo 2003/03/12
//
// Modified on 2010/08/05 for commenting out ASSERT definition causing
// warning during compilation because of redefinition

#ifndef C3D_LIST_DEFS
#define C3D_LIST_DEFS

//#define ASSERT _ASSERT

// multiplatform condition for wxWidgets GCC vs. Visual Studio
#ifdef MS_VS
// must disable this warning because of operator->
// defined in C3DList<>. The operator can be used only
// with template arguments which provide an operator->
// on their side; e.g. C3Dlist<int*> cannot use
// operator-> (does not make sense).
// The warning is a reminder of this situation.
#pragma warning( disable : 4284 )
#endif

// double linked list base class
template<class T> class C3DList
{
    // the friend function is fully defined inside the host class template definition,
    // to avoid the warning about "friend declaration declares a non-template function":
    // as a matter of fact, SwapList() is not a template function in itself;
    // it just uses the class C3DList<t> as argument. Now, every time the C3DList class template
    // is instantiated, a new, ordinary function overload is created that takes an argument
    // of the current C3DList specialization (i.e. with 'T' replaced by the required type).
    // If we want SwapList() to be a friend template instead, we would have needed angle brackets
    // in the declaration of SwapList() inside C3DList (i.e. SwapList<>(C3DList... ). This is
    // necessary to tell the compiler that SwapList() is a template. Otherwise, the compiler
    // will look for an ordinary function named SwapList() and not find it. However, different
    // compilers behave differently in this respect; and MS VS w.r.t. MinGW does not work correctly.
    // Anyway, we define SwapList() as a friend template, we also need a forward declaration of the
    // function template SwapList() before the C3DList class definition: the language specifies that
    // friend function templates must be previously declared. However, to properly declare SwapList(),
    // C3DList must also have been declared, since SwapList() takes a C3DList argument, hence the
    // forward declaration of C3DList in the beginning
	friend void SwapLists(C3DList<T> &lista, C3DList<T> &listb)
	{
        class C3DList<T>::Cell *tmpHead, *tmpTail;
        long tmpElenum;

        tmpHead = lista.head;
        tmpTail = lista.tail;
        tmpElenum = lista.elenum;

        lista.head = listb.head;
        lista.tail = listb.tail;
        lista.elenum = listb.elenum;

        listb.head = tmpHead;
        listb.tail = tmpTail;
        listb.elenum = tmpElenum;
    }

// debug: must remove following public declaration for Cell
public:
	// the basic list cell
	class Cell
	{
	public:
		Cell *next;
		Cell *prev;
		T data;

		Cell(T dat, Cell *nxt, Cell *prv) : data(dat), next(nxt), prev(prv) {}
		Cell(T dat) : data(dat), next(NULL), prev(NULL) {}
	};

	// head of the list
	class Cell *head;
	// tail of the list
	class Cell *tail;
	// number of elements in list
	long elenum;

public:
	inline C3DList() : head(NULL), tail(NULL), elenum(0) {}
	inline ~C3DList() {clear();}

	// declaration required
	class iterator;
	// make it friend
	friend class iterator;
	// definition
	class iterator
	{
	public:
		Cell *p;
		// constructors
		iterator(Cell* cell) : p(cell) {}
		// copy-constructor
		iterator(const iterator& it) : p(it.p) {}
		// end sentinel iterator
		iterator() : p(NULL) {}

		inline T *getPointer()
		{
			if (!p)
				return NULL;
			else
				return (&(p->data));
		}

		inline void assign(T newdata)
		{
			p->data = newdata;
		}

		inline void operator++()
		{
			if(p != NULL)
				p = p->next;
		}

		inline void operator++(int)
		{
			operator++();
		}

		inline void operator--()
		{
			if(p != NULL)
				p = p->prev;
		}

		inline void operator--(int)
		{
			operator--();
		}

		inline T operator*()
		{
			ASSERT(p);
			return (p->data);
		}

		inline T* operator->() const
		{
			if (!p)
				return NULL;
			else
				return (&(p->data));
		}

		inline bool operator==(const iterator &it)
		{
			return (p == it.p);
		}

		inline bool operator!=(const iterator &it)
		{
			return (p != it.p);
		}
	};


	inline iterator begin()
	{
		iterator it(head);

		return it;
	}

	inline iterator end()
	{
		iterator it;

		// end iterator points to NULL
		return it;
	}

	inline iterator rbegin()
	{
		iterator it(tail);

		return it;
	}

	inline iterator rend()
	{
		iterator it;

		// end iterator points to NULL
		return it;
	}

	inline void push_front(T dat)
	{
		Cell *tmp;

		// if first element in list
		if( head == NULL ) {
			ASSERT(tail == NULL);
			head = tail = new Cell(dat);
		}
		else {
			tmp = head;
			// create new head cell, shifting old head ahead
			head = new Cell(dat, head, NULL);
			// adjust prev pointer of old head
			tmp->prev = head;
		}

		// increment the element counter
		elenum++;
	}

	inline void push_back(T dat)
	{
		Cell *tmp;

		// if first element in list
		if( head == NULL ) {
			ASSERT(tail == NULL);
			head = tail = new Cell(dat);
		}
		else {
			tmp = tail;
			// create new tail cell, shifting old tail backward
			tail = new Cell(dat, NULL, tail);
			// adjust next pointer of old head
			tmp->next = tail;
		}

		// increment the element counter
		elenum++;
	}

	inline T pop_front()
	{
		T tmpData;
		Cell *tmpCell;

//		if(head != NULL) {

			// save data to return
			tmpData = head->data;

			// if only one element in list
			if(elenum == 1) {
				ASSERT(head->next == NULL);
				ASSERT(head->prev == NULL);
				ASSERT(head == tail);
				delete head;
				head = tail = NULL;
			}
			else {
				tmpCell = head;
				// set new head
				head = head->next;
				// null prev link (is new head)
				head->prev = NULL;
				// delete old head
				delete tmpCell;
			}

			// decrement the element counter
			elenum--;

			return tmpData;
//		}
//		else
//			return NULL;
	}

	inline T pop_back()
	{
		T tmpData;
		Cell *tmpCell;

//		if(tail != NULL) {

			// save data to return
			tmpData = tail->data;

			// if only one element in list
			if(elenum == 1) {
				ASSERT(tail->next == NULL);
				ASSERT(tail->prev == NULL);
				ASSERT(head == tail);
				delete tail;
				head = tail = NULL;
			}
			else {
				tmpCell = tail;
				// set new head
				tail = tail->prev;
				// null prev link (is new head)
				tail->next = NULL;
				// delete old head
				delete tmpCell;
			}

			// decrement the element counter
			elenum--;

			return tmpData;
//		}
//		else
//			return NULL;
	}

	inline T pop(iterator &it)
	{
		T tmpData;
		Cell *tmpCell;

		tmpCell = it.p;

//		if(tmpCell != NULL) {
			tmpData = tmpCell->data;

			// if only one element in list
			if(elenum == 1) {
				ASSERT(tmpCell == head);
				ASSERT(head == tail);
				delete head;
				head = tail = NULL;
			}
			else if (tmpCell == head) {
				ASSERT(head->next);
				// shift head forward
				head = head->next;
				head->prev = NULL;
				// delete old head
				delete tmpCell;
			}
			else if (tmpCell == tail) {
				ASSERT(tail->prev);
				// shift tail backward
				tail = tail->prev;
				tail->next = NULL;
				// delete old tail
				delete tmpCell;
			}
			else {
				ASSERT(tmpCell->next);
				ASSERT(tmpCell->prev);
				// relink prev and next cells to skip cell to be popped
				tmpCell->prev->next = tmpCell->next;
				tmpCell->next->prev = tmpCell->prev;
				// delete old head
				delete tmpCell;
			}

			// decrement the element counter
			elenum--;

			return tmpData;
//		}
//		else
//			return NULL;
	}

	inline void clear()
	{
		Cell *cell, *tmp;

		// fast erase: no need to relink cells,
		// since are all deleted
		for(cell = head; cell != NULL; cell = tmp) {
			tmp = cell->next;
			delete cell;
		}

		head = tail = NULL;
		elenum = 0;
	}

	// return the length of the controlled sequence
	inline long size()
	{
		return elenum;
	}

	// return the length of the controlled sequence
	inline bool is_empty()
	{
		return (elenum == 0);
	}
/*
	// find an element in the controlled sequence
	inline iterator find(T dat)
	{
		iterator it;

		for(it = begin(); it != end(); it++)
			if( *it == dat)
				return it;

		// not found
		return NULL;
	}
*/
	// merge two sequences
	inline void Merge(C3DList<T> &listToBeMerged)
	{
		while( listToBeMerged.is_empty() == false )
			push_back( listToBeMerged.pop_front() );
	}

	// copy a sequence
	inline void Copy(C3DList<T> &listToBeCopied)
	{
		iterator it;

		clear();

		for(it = listToBeCopied.begin(); it != listToBeCopied.end(); it++)
			push_back(*it);
	}
};

#endif //C3D_LIST_DEFS

