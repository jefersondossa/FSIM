//This cpp file describes the fIncrem_Transient function (Analysis Menu from FSArl).

#include "WindowConstructor.h"

void WindowConstructor::fIncrem_Transient(Fl_Group *group){

    Increm_Transient_Menu = new Fl_Group(1165, 175, 370, 240);
    Increm_Transient_Menu->hide();

    Nsteps = new Fl_Value_Input(1325, 185, 95, 25, "Steps Number: ");
    tolerance = new Fl_Value_Input(1325, 225, 95, 25, "Tolerance: ");
    maxInterations = new Fl_Value_Input(1325, 265, 95, 25, "Max. Interations: ");

    Apply_Analysis = new Fl_Button(1430, 305, 80, 25, "Apply");
    Apply_Analysis->callback((Fl_Callback*)fStatic_Apply_Analysis);

    Increm_Transient_Menu->end();

    group->add(Increm_Transient_Menu);
}