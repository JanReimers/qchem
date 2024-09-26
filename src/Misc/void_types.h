// File: void_types.h  Determine void* types for STL containers.
#ifndef _void_types_h_
#define _void_types_h_

// Copyright (1994-2003), Jan N. Reimers

//
//  Need some template specialized classes to determine the void* type.
//
template <class T> struct VoidType;

template <class T> struct VoidType<T*>
{
  typedef void* void_type;
};

template <class T> struct VoidType<const T*>
{
  typedef const void* void_type;
};


#endif //_void_types_h_
