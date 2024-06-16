//This cpp file describes the fElasticity3D function (Elasticity2D Submenu from Weak Form Menu).

#include "WindowConstructor.h"

void WindowConstructor::fElasticity3D(Fl_Group *group){

  Elasticity3D_Menu = new Fl_Group(1175, 175, 360, 175);
  Elasticity3D_Menu->hide();

  young = new Fl_Value_Input(1310, 190, 95, 25, "Elasticity Modulus: ");

  poisson = new Fl_Value_Input(1310, 235, 95, 25, "Poisson Ratio: ");

  /*============================== Button_Apply_WeakForm ==============================*/

  Button_Apply_WeakForm = new Fl_Button(1410, 275, 80, 25, "Apply");
  Button_Apply_WeakForm->callback((Fl_Callback*)fStatic_Apply_WF);


  Elasticity3D_Menu->end();

  group->add(Elasticity3D_Menu);

  this->WeakForm->add("Elasticity 3D",0,(Fl_Callback*)WindowConstructor::fStatic_Elasticity3D,0,0);

}