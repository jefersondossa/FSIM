/* This cpp file describes the following Pos-Processing Menu Inline and Static functions:

1.  From fInput_vtk_cb:
        fInline_Input_vtk;
        fStatic_Input_vtk;

*/

#include "WindowConstructor.h"

/* ========================= From fInput_vtk_cb ========================= */

void WindowConstructor::fInline_Input_vtk(Fl_Output*, void*) {
  this->fInput_vtk_cb();
}

void WindowConstructor::fStatic_Input_vtk(Fl_Output* o, void* v) {
  ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_Input_vtk(o,v);
}