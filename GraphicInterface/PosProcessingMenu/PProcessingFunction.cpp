/* This cpp file describes the following Pos-Processing Menu functions:

1.  fInput_vtk_cb

*/

#include "WindowConstructor.h"

/* ========================= fInput_vtk_cb ========================= */

void WindowConstructor::fInput_vtk_cb(){

    Fl_File_Chooser GP(".", "*.vtk", Fl_File_Chooser::SINGLE, "Select a .vtk file");

    string previous_vtk = vtk->value(); //Stores the previous value.

    GP.show(); //Opens the file chooser
    while (GP.shown()){
        Fl::wait();
        if (GP.value() != nullptr) {
            this->vtk->value(GP.value());
        }

        else{
            this->vtk->value(previous_vtk.c_str()); 
        }
    }
}