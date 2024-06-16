//This cpp file describes the fPoisson function (Poisson Submenu from Weak Form Menu).

#include "WindowConstructor.h"

void WindowConstructor::fPoisson(Fl_Group *group){

  Poisson_Menu = new Fl_Group(1175, 175, 360, 175);
  Poisson_Menu->hide();

  Nsteps = new Fl_Value_Input(1368, 230, 55, 25, "Number of State Variables: ");

  /*============================== Button_Apply_WeakForm ==============================*/

  Button_Apply_WeakForm = new Fl_Button(1410, 270, 80, 25, "Apply");
  Button_Apply_WeakForm->callback((Fl_Callback*)fStatic_Apply_WF);

  Poisson_Menu->end();

  group->add(Poisson_Menu);

  this->WeakForm->add("Poisson",0,(Fl_Callback*)WindowConstructor::fStatic_Poisson,0,0);

}