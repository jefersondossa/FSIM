//This cpp file describes the fNavierStokes function (NavierStokes Submenu from Weak Form Menu).

#include "WindowConstructor.h"

void WindowConstructor::fNavierStokes(Fl_Group *group){

  NavierStokes_Menu = new Fl_Group(1175, 175, 360, 175);
  NavierStokes_Menu->hide();

  density = new Fl_Value_Input(1310, 190, 95, 25, "Density: ");

  viscosity = new Fl_Value_Input(1310, 235, 95, 25, "Viscosity: ");

  /*============================== Button_Apply_WeakForm ==============================*/

  Button_Apply_WeakForm = new Fl_Button(1410, 275, 80, 25, "Apply");
  Button_Apply_WeakForm->callback((Fl_Callback*)fStatic_Apply_WF);

  NavierStokes_Menu->end();

  group->add(NavierStokes_Menu);

  this->WeakForm->add("Navier Stokes",0,(Fl_Callback*)WindowConstructor::fStatic_NavierStokes,0,0);
}