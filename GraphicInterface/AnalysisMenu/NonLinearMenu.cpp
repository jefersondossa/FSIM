//This cpp file describes the fNonLinear function (Analysis Menu from FSArl).

#include "WindowConstructor.h"

void WindowConstructor::fNonLinear(Fl_Group *group){

    NonLinear_Menu = new Fl_Group(1165, 175, 370, 240);
    NonLinear_Menu->hide();

    tolerance = new Fl_Value_Input(1325, 185, 95, 25, "Tolerance: ");
    maxInterations = new Fl_Value_Input(1325, 225, 95, 25, "Max. Interations: ");

    NonLinear_Menu->end();

    group->add(NonLinear_Menu);
}