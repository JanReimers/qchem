// File: Unpickle.H  Template functions for unpickling objects.
#ifndef _Unpickle_H_
#define _Unpickle_H_

#include <string>
#include <fstream>
#include <iostream>
#include <cassert>

template <class T> bool UnPickle(T*& pointer,const char* filep, const char* name)
{
    std::string file(filep);
    bool file_error=true;
    if(file !="")
    {
        std::ifstream in(file.c_str());
        if(!in)
        {
            std::cerr << "Can't open " << name << " file :" << file << std::endl;
            file_error=false;
        }
        else
        {
            pointer = T::Factory(in);
            assert(pointer);
            in >> *pointer;
            assert(in);
        }
    }
    return file_error;
}



#endif // _Unpickle_H_
