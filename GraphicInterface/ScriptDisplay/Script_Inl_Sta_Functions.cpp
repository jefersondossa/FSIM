/* This cpp file describes the following Script Display Menu Inline and Static functions:

1.  From fButton_SaveScript_cb:
        fInline_Button_SaveScript
        fStatic_Button_SaveScript

2.  From fButton_TextSave_cb:
        fInline_Button_TextSave
        fStatic_Button_TextSave

3.  From fButton_TextCancel_cb:
        fInline_Button_TextCancel
        fStatic_Button_TextCancel

*/

#include "WindowConstructor.h"

/* ========================= From fButton_SaveScript_cb ========================= */

void WindowConstructor::fInline_Button_SaveScript(Fl_Button*, void*){
    this->fButton_SaveScript_cb();

}

void WindowConstructor::fStatic_Button_SaveScript(Fl_Button* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_Button_SaveScript(o,v);

}

/* ========================= From fButton_TextSave_cb ========================= */

void WindowConstructor::fInline_Button_TextSave(Fl_Button*, void*){
    this->fButton_TextSave_cb();

}

void WindowConstructor::fStatic_Button_TextSave(Fl_Button* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_Button_TextSave(o,v);

}

/* ========================= From fButton_TextCancel_cb ========================= */

void WindowConstructor::fInline_Button_TextCancel(Fl_Button*, void*){
    this->fButton_TextCancel_cb();

}

void WindowConstructor::fStatic_Button_TextCancel(Fl_Button* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_Button_TextCancel(o,v);


}