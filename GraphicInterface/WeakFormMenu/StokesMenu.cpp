//This cpp file describes the fStokes function (Stokes Submenu from Weak Form Menu).

#include "WindowConstructor.h"

void WindowConstructor::fStokes(Fl_Group *group){

  Stokes_Menu = new Fl_Group(1175, 175, 360, 175);
  Stokes_Menu->hide();

  density = new Fl_Value_Input(1310, 190, 95, 25, "Density: ");

  viscosity = new Fl_Value_Input(1310, 235, 95, 25, "Viscosity: ");

  /*============================== Button_Apply_WeakForm ==============================*/

  Button_Apply_WeakForm = new Fl_Button(1410, 275, 80, 25, "Apply");
  Button_Apply_WeakForm->callback((Fl_Callback*)fStatic_Apply_WF);

  Stokes_Menu->end();

  group->add(Stokes_Menu);

  this->WeakForm->add("Stokes",0,(Fl_Callback*)WindowConstructor::fStatic_Stokes,0,0);

}