/* This cpp file describes the following Top Header Menu Inline and Static functions:

1.  From fOpen_File_cb:
        fInline_Open_File
        fStatic_Open_File

2.  From fNew_File_cb:
        fInline_New_File
        fStatic_New_File

3.  From fGenerateHeader_cb:
        fInline_G_HeaderFile
        fStatic_G_HeaderFile

3.  From fGenerateFoam_cb:
        fInline_G_FoamFile
        fStatic_G_FoamFile

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

/* ========================= From fGenerateHeader_cb ========================= */

void WindowConstructor::fInline_G_HeaderFile(Fl_Menu_*, void*) {
    this->fGenerateHeader_cb();
}

void WindowConstructor::fStatic_G_HeaderFile(Fl_Menu_* o, void* v) {
  ((WindowConstructor*)(o->parent()->user_data()))->fInline_G_HeaderFile(o,v);
}

/* ========================= From fGenerateFoam_cb ========================= */

void WindowConstructor::fInline_G_FoamFile(Fl_Menu_*, void*) {
    this->fGenerateFoam_cb();
}

void WindowConstructor::fStatic_G_FoamFile(Fl_Menu_* o, void* v) {
  ((WindowConstructor*)(o->parent()->user_data()))->fInline_G_FoamFile(o,v);
}