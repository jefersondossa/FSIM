// Graphical Interface (GUI).cpp : This file contains the 'main' function. Program execution begins and ends there.
// Remember to add Window->show() in cpp file

#include "WindowConstructor.h"
#include <iostream>
#include "PainelControle.h"

#ifdef HAS_GRAPHINTERFACE

int main()
{  
    WindowConstructor c;
    c.Open();
    //c.EP2D_f();
    // c.Show();

    // PainelControle controle;
    // controle.iniciar();

    return Fl::run();
}
#endif