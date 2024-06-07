/* This cpp file describes the following Top Header Menu Inline and Static functions:

1.  From fOpen_File_cb:
        fInline_Open_File
        fStatic_Open_File

2.  From fNew_File_cb:
        fStatic_New_File
        fStatic_New_File

*/

#include "WindowConstructor.h"

/* ========================= From fOpen_File_cb ========================= */

void WindowConstructor::fInline_Open_File(Fl_Menu_*, void*){
    this->fOpen_File_cb();

}

void WindowConstructor::fStatic_Open_File(Fl_Menu_* o, void* v){
    ((WindowConstructor*)(o->parent()->user_data()))->fInline_Open_File(o,v);

}

/* ========================= From fNew_File_cb ========================= */

void WindowConstructor::fInline_New_File(Fl_Menu_*, void*) {
    this->fNew_File_cb();
}

void WindowConstructor::fStatic_New_File(Fl_Menu_* o, void* v) {
  ((WindowConstructor*)(o->parent()->user_data()))->fInline_New_File(o,v);
}