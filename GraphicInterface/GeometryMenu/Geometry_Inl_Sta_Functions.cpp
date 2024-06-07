/* This cpp file describes the following Geometry Menu Inline and Static functions:

1.  From fWF_free_cb:
        fInline_WF_free
        fStatic_WF_free

2.  From fInput_geo_cb:
        fInline_Input_geo
        fStatic_Input_geo

3.  From fInput_msh_cb:
        fInline_Input_msh
        fStatic_Input_msh

4.  From fApply_geometry_cb:
        fInline_Apply_geometry
        fStatic_Apply_geometry

*/

#include "WindowConstructor.h"

/* ========================= From fWF_free_cb ========================= */

void WindowConstructor::fInline_WF_free(Fl_Output*, void*){
  this->fWF_free_cb();
}
void WindowConstructor::fStatic_WF_free(Fl_Output* o, void* v){
  ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_WF_free(o,v);

}

/* ========================= From fInput_geo_cb ========================= */

void WindowConstructor::fInline_Input_geo(Fl_Output*, void*) {
  this->fInput_geo_cb();
}

void WindowConstructor::fStatic_Input_geo(Fl_Output* o, void* v) {
  ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_Input_geo(o,v);
}

/* ========================= From fInput_msh_cb ========================= */

void WindowConstructor::fInline_Input_msh(Fl_Output*, void*) {
  this->fInput_msh_cb();
}

void WindowConstructor::fStatic_Input_msh(Fl_Output* o, void* v) {
  ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_Input_msh(o,v);
}

/* ========================= From fApply_geometry_cb ========================= */

void WindowConstructor::fInline_Apply_geometry(Fl_Button*, void*) {
  this->fApply_geometry_cb();
}

void WindowConstructor::fStatic_Apply_geometry(Fl_Button* o, void* v) {
  ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_Apply_geometry(o,v);
}
