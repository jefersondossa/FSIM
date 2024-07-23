//This cpp file describes the fGeometry function (Geometry Menu from FSArl).

#include "WindowConstructor.h"

void WindowConstructor::fGeometry(Fl_Double_Window *Window){

    Geometry_Menu = new Fl_Group(1160, 75, 375, 290, "Geometry");
    Geometry_Menu->color(FL_LIGHT1);
    Geometry_Menu->labelsize(18);
    Geometry_Menu->show();

    geo = new Fl_Output(1190, 115, 325, 35, "Input .geo Path: ");
    geo->align(Fl_Align(FL_ALIGN_TOP_LEFT));
    geo->value(""); //Sets geo inicial value to nothing.
    geo->callback((Fl_Callback*)fStatic_Input_geo); //Sets the geo Output callback

    msh = new Fl_Output(1190, 190, 325, 35, "Input .msh Path:");
    msh->align(Fl_Align(FL_ALIGN_TOP_LEFT));
    this->msh->value(""); //Sets msh inicial value to nothing.
    msh->callback((Fl_Callback*)fStatic_Input_msh); //Sets the msh Output callback

    dimension = new Fl_Value_Output(1270, 250, 35, 25, "Dimension: ");
    dimension->color(FL_BACKGROUND2_COLOR);

    Apply_geometry = new Fl_Button(1425, 250, 90, 25, "Save paths");
    Apply_geometry->callback((Fl_Callback*)fStatic_Apply_geometry);

    Deactivate_Render = new Fl_Check_Button(1270, 300, 90, 25, "Render");

    Geometry_Menu->end();
    Window->add(Geometry_Menu);
 }