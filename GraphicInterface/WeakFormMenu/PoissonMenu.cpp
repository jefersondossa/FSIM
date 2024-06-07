//This cpp file describes the fPoisson function (Poisson Submenu from Weak Form Menu).

#include "WindowConstructor.h"

void WindowConstructor::fPoisson(Fl_Group *group){

  Poisson_Menu = new Fl_Group(1175, 175, 360, 175);
  Poisson_Menu->hide();

  Nsteps = new Fl_Value_Input(1368, 230, 55, 25, "Number of State Variables: ");

  Poisson_Menu->end();

  group->add(Poisson_Menu);

  this->WeakForm->add("Poisson",0,(Fl_Callback*)WindowConstructor::fStatic_Poisson,0,0);

}