// Graphical Interface (GUI).cpp : This file contains the 'main' function. Program execution begins and ends there.
// Remember to add Window->show() in cpp file

#include "WindowConstructor.h"
#include <iostream>


int main()
{  
#ifdef HAS_GRAPHINTERFACE

    WindowConstructor c;
    c.f_Open();

    return Fl::run();
#else 
    return 0;
#endif

}
