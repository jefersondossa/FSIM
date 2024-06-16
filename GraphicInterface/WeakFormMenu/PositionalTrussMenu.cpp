//This cpp file describes the fPositionalTruss function (PositionalTruss Submenu from Weak Form Menu).

#include "WindowConstructor.h"

void WindowConstructor::fPositionalTruss(Fl_Group *group){

  PositionalTruss_Menu = new Fl_Group(1175, 175, 360, 175);
  PositionalTruss_Menu->hide();

  young = new Fl_Value_Input(1310, 190, 95, 25, "Elasticity Modulus: ");

  area = new Fl_Value_Input(1310, 235, 95, 25, "Area: ");

  /*============================== Button_Apply_WeakForm ==============================*/

  Button_Apply_WeakForm = new Fl_Button(1410, 275, 80, 25, "Apply");
  Button_Apply_WeakForm->callback((Fl_Callback*)fStatic_Apply_WF);

  PositionalTruss_Menu->end();

  group->add(PositionalTruss_Menu);

  this->WeakForm->add("Positional Truss",0,(Fl_Callback*)WindowConstructor::fStatic_PositionalTruss,0,0);

}