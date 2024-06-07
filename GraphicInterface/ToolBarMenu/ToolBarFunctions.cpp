/* This cpp file describes the following Tool Bar Menu functions:

1.  fButton_Geometry_cb
2.  fButton_WeakForm_cb
3.  fButton_Analysis_cb
4.  fButton_PosProcressing_cb

*/

#include "WindowConstructor.h"

/* ========================= fButton_Geometry_cb ========================= */

void WindowConstructor::fButton_Geometry_cb(){

    Geometry_Menu->show();
    WF_Menu->hide();
    Analysis_Menu->hide();
    PProcessing_Menu->hide();
}

/* ========================= fButton_WeakForm_cb ========================= */

void WindowConstructor::fButton_WeakForm_cb(){

    Geometry_Menu->hide();
    WF_Menu->show();
    Analysis_Menu->hide();
    PProcessing_Menu->hide();
}

/* ========================= fButton_Analysis_cb ========================= */

void WindowConstructor::fButton_Analysis_cb(){

    Geometry_Menu->hide();
    WF_Menu->hide();
    Analysis_Menu->show();
    PProcessing_Menu->hide();
}

/* ========================= fButton_PosProcressing_cb ========================= */

void WindowConstructor::fButton_PosProcressing_cb(){

    Geometry_Menu->hide();
    WF_Menu->hide();
    Analysis_Menu->hide();
    PProcessing_Menu->show();
}