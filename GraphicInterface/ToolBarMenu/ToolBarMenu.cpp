//This cpp file describes the fToolBar function (Tool Bar Menu from FSArl).

#include "WindowConstructor.h"

void WindowConstructor::fToolBar(Fl_Double_Window *window){

    ToolBar_Menu = new Fl_Group(0, 35, 1160, 45);
    ToolBar_Menu->box(FL_BORDER_BOX);
    ToolBar_Menu->color(FL_BACKGROUND2_COLOR);
    ToolBar_Menu->show();

    fButton_Geometry(ToolBar_Menu);
    fButton_WeakForm(ToolBar_Menu);
    fButton_Analysis(ToolBar_Menu);
    fButton_PosProcressing(ToolBar_Menu);

    ToolBar_Menu->end();

    window->add(ToolBar_Menu);
}