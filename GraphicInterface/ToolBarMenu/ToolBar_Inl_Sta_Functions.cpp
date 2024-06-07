/* This cpp file describes the following Tool Bar Menu Inline and Static functions:

1.  From fButton_Geometry_cb:
        fInline_Button_Geometry
        fStatic_Button_Geometry

2.  From fButton_WeakForm_cb:
        fInline_Button_WeakForm
        fStatic_Button_WeakForm

3.  From fButton_Analysis_cb:
        fInline_Button_Analysis
        fStatic_Button_Analysis

4.  From fButton_PosProcressing_cb:
        fInline_Button_Analysis
        fStatic_Button_Analysis

*/

#include "WindowConstructor.h"

/* ========================= From fButton_Geometry_cb ========================= */

void WindowConstructor::fInline_Button_Geometry(Fl_Button*, void*){
    this->fButton_Geometry_cb();

}

void WindowConstructor::fStatic_Button_Geometry(Fl_Button* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_Button_Geometry(o,v);

}

/* ========================= From fButton_WeakForm_cb ========================= */

void WindowConstructor::fInline_Button_WeakForm(Fl_Button*, void*){
    this->fButton_WeakForm_cb();

}
void WindowConstructor::fStatic_Button_WeakForm(Fl_Button* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_Button_WeakForm(o,v);

}

/* ========================= From fButton_Analysis_cb ========================= */

void WindowConstructor::fInline_Button_Analysis(Fl_Button*, void*){
    this->fButton_Analysis_cb();

}
void WindowConstructor::fStatic_Button_Analysis(Fl_Button* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_Button_Analysis(o,v);

}

/* ========================= From fButton_PosProcressing_cb ========================= */

void WindowConstructor::fInline_Button_PProcessing(Fl_Button*, void*){
    this->fButton_PosProcressing_cb();

}
void WindowConstructor::fStatic_Button_PPospricessing(Fl_Button* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_Button_PProcessing(o,v);

}