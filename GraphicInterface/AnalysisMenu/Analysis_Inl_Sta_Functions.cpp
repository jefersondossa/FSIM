/* This cpp file describes the following Analysis Menu Inline and Static functions:

1.  From fLinear_cb:
        fInline_Linear;
        fStatic_Linear;

2.  From fNonLinear_cb:
        fInline_NonLinear;
        fStatic_NonLinear;

3.  From fIncrem_Transient_cb:
        fInline_Increm_Transient;
        fStatic_Increm_Transient;

4.  From fApply_Analysis_cb:
        fInline_Apply_Analysis;
        fStatic_Apply_Analysis;

5.  From fOF_Analysis_cb:
        fInline_OF_Analysis;
        fStatic_OF_Analysis;

*/

#include "WindowConstructor.h"

/* ========================= From fLinear_cb ========================= */

void WindowConstructor::fInline_Linear(Fl_Menu_*, void*){
    this->fLinear_cb();

}

void WindowConstructor::fStatic_Linear(Fl_Menu_* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_Linear(o,v);

}

/* ========================= From fNonLinear_cb ========================= */

void WindowConstructor::fInline_NonLinear(Fl_Menu_*, void*){
    this->fNonLinear_cb();

}

void WindowConstructor::fStatic_NonLinear(Fl_Menu_* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_NonLinear(o,v);

}

/* ========================= From fIncrem_Transient_cb ========================= */

void WindowConstructor::fInline_Increm_Transient(Fl_Menu_*, void*){
    this->fIncrem_Transient_cb();

}

void WindowConstructor::fStatic_Increm_Transient(Fl_Menu_* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_Increm_Transient(o,v);

}

/* ========================= From fOF_Analysis_cb ========================= */

void WindowConstructor::fInline_OF_Analysis(Fl_Menu_*, void*){
    this->fOFAnalysis_cb();

}

void WindowConstructor::fStatic_OF_Analysis(Fl_Menu_* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_OF_Analysis(o,v);

}

/* ========================= From fApply_Analysis_cb ========================= */

void WindowConstructor::fInline_Apply_Analysis(Fl_Button*, void*){
    this->fApply_Analysis_cb();

}

void WindowConstructor::fStatic_Apply_Analysis(Fl_Button* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->parent()->user_data()))->fInline_Apply_Analysis(o,v);

}
