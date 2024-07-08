//This cpp file describes the fPProcessing function (Pos´Processing Menu from FSArl).

#include "WindowConstructor.h"

void WindowConstructor::fPProcessing(Fl_Double_Window *Window){

    PProcessing_Menu = new Fl_Group(1160, 75, 375, 290, "PProcessing");
    PProcessing_Menu->color(FL_LIGHT1);
    PProcessing_Menu->labelsize(18);
    PProcessing_Menu->hide();

    vtk = new Fl_Output(1190, 115, 325, 35, "Input .vtk Path: ");
    vtk->align(Fl_Align(FL_ALIGN_TOP_LEFT));
    vtk->value(""); //Sets vtk inicial value to nothing.
    vtk->callback((Fl_Callback*)fStatic_Input_vtk); //Sets the vtk Output callback

    PProcessing_Menu->end();
    Window->add(PProcessing_Menu);
 }