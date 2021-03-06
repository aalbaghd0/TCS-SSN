/* this file implements the class Heap */

#include "tcs_ssn_h/heap.h"
#include "tcs_ssn_h/gendef.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include<iostream>

//------------------------------------------
HeapEntry::HeapEntry()
{
}
//------------------------------------------
HeapEntry::~HeapEntry()
{
}
//------------------------------------------
void HeapEntry::init_HeapEntry(int _dim)
{
	dim = _dim;
}
//------------------------------------------
void HeapEntry::copy(HeapEntry *_he)
{
	key = _he -> key;
	level = _he -> level;
	son1 = _he -> son1;
	son2 = _he -> son2;
}
//------------------------------------------
//------------------------------------------
//------------------------------------------
Heap::Heap()
{
	cont = NULL;
}
//------------------------------------------
Heap::~Heap()
{
	//printf("maximum heap entry=%d\n", maxused);
	delete [] cont;
	cont = NULL;
}
//------------------------------------------
void Heap::enter(HeapEntry *_he, int _pos)
//this function enters a new entry into the heap at position _pos
{

	for (int i = used - 1; i >= _pos; i --)
	{
		cont[i + 1].copy(&(cont[i]));
	}
	cont[_pos].copy(_he);
	used ++;

	if (maxused<used)
		maxused=used;
}
//------------------------------------------
void Heap::insert(HeapEntry *_he)
{
	int pos = used;  //pos is the position _he will be inserted

	enter(_he, pos);
	// now perform swapping
	pos++;
	int parent = pos;
	while (parent != 1)
	{
		int child = parent;
		parent /= 2;
		if (cont[parent - 1].key > cont[child - 1].key)
		{
			HeapEntry *the = new HeapEntry();
			the -> init_HeapEntry(cont[parent - 1].dim);
			the -> copy(&(cont[parent - 1]));
			cont[parent - 1].copy(&(cont[child - 1]));
			cont[child - 1].copy(the);
			delete the;
		}
		else 
			parent = 1;
	}

	if (used > hsize)  
		//this is why the heap size is initiated with one more place than hsize (to facilitate
		//coding)
	{
		std::cerr<<"heap exceeded...\n";
	}

	//if (!check())
	//	printf("testing...\n");
}


//------------------------------------------
bool Heap::remove(HeapEntry *_he)
//this function deheaps an entry. the return value indicates whether successful: false
//means heap is already empty
{
	if (used==0) 
		return false;
	_he -> copy(&(cont[0]));
	used--;
	cont[0].copy(&(cont[used]));
	int parent = 1;
	while (2 * parent <= used)
	{
		int child = 2 * parent;
		if (2 * parent + 1 > used)
			child = 2 * parent;
		else
			if (cont[2 * parent - 1].key < cont[2 * parent].key)
				child = 2 * parent;
			else 
				child = 2 * parent + 1;

		if (cont[parent - 1].key > cont[child - 1].key)
		{
			HeapEntry *the = new HeapEntry();
			the -> init_HeapEntry(cont[parent - 1].dim);
			the -> copy(&(cont[parent - 1]));
			cont[parent - 1].copy(&(cont[child - 1]));
			cont[child - 1].copy(the);
			delete the;
			parent = child; 
		}
		else
			parent = used;
	}

	//if (!check())
	//	printf("testing...\n");
	return true;
};
//------------------------------------------

bool Heap::deleteEntry(int a) {
	// this function delete a node with specific value passed a
	int pos;
	if (used == 0)
		return false;

	for (int i = 0; i < used; i++) {
		if (cont[i].son1 == a) {
			pos = i;
			break;
		}
		if (i == used - 1) {// the element is not there
			std::cerr << "THe element is not in the heap \n";
			return false;
		}
	}

	// copy the last element and put it in the deleted one
	// reuce the size of the heap --used
	cont[pos].copy(&(cont[used - 1]));
	used--;
	
label1:
	int parent = (pos - 1) / 2;
	int child1 = (pos * 2) + 1;
	int child2 = (pos * 2) + 1 + 1;

	if(cont[pos].key < cont[parent].key){
		HeapEntry* the = new HeapEntry();
		the->init_HeapEntry(cont[parent].dim);
				the->copy(&(cont[parent]));
				cont[parent].copy(&(cont[pos]));
				cont[pos].copy(the);
				delete the;
				goto label1;
	}
	else {
		if (parent == 0)
			parent = 1;
		while (2 * parent <= used) {
			int child = 2 * parent;
			if (2 * parent + 1 > used)
				child = 2 * parent;
			else
				if (cont[2 * parent - 1].key < cont[2 * parent].key)
					child = 2 * parent;
				else
					child = 2 * parent + 1;

			if (cont[parent - 1].key > cont[child - 1].key) {
				HeapEntry* the = new HeapEntry();
				the->init_HeapEntry(cont[parent - 1].dim);
				the->copy(&(cont[parent - 1]));
				cont[parent - 1].copy(&(cont[child - 1]));
				cont[child - 1].copy(the);
				delete the;
				parent = child;
			}
			else
				parent = used;
		}
	}
	
	// does not work
}
void Heap::clean(double _dist)
//this function cleans those entries with keys greater thab _dist
//in the memory
{
	for (int i = 0; i < used; i ++)
	{
		if (cont[i].key > _dist)
			used = i;
	}
}

/*****************************************************************
this function checks the integrity of the heap. it is an auxiliary
function for debugging

Coded by Yufei Tao 09/01/02
*****************************************************************/

bool Heap::check()
{
	for (int i = 0; i < used; i ++)
	{
		if (cont[i].son1<0 || cont[i].son2<0)
			return false;
	}
	return true;
}
//------------------------------------------
void Heap::init(int _dim, int _hsize)
{
	if (cont)
		delete [] cont;
	hsize = _hsize;
	cont = new HeapEntry [hsize + 1];   // +1 is to facilitate insert
	for (int i = 0; i < hsize + 1; i ++)
		cont[i].init_HeapEntry(_dim);
	used = 0;
	maxused=0;
}
//------------------------------------------
	

