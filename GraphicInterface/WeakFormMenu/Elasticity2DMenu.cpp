//This cpp file describes the fElasticity2D function (Elasticity2D Submenu from Weak Form Menu).

#include "WindowConstructor.h"

void WindowConstructor::fElasticity2D(Fl_Group *group){

  Elasticity2D_Menu = new Fl_Group(1175, 175, 360, 175);
  Elasticity2D_Menu->hide();

  young = new Fl_Value_Input(1310, 190, 95, 25, "Elasticity Modulus: ");

  poisson = new Fl_Value_Input(1310, 235, 95, 25, "Poisson Ratio: ");

  Fl_Check_Button *plane_stress = new Fl_Check_Button(1200, 275, 45, 25, "Plane Stress ");
  plane_stress->down_box(FL_DOWN_BOX);

  Elasticity2D_Menu->end();

  group->add(Elasticity2D_Menu);

  this->WeakForm->add("Elasticity 2D",0,(Fl_Callback*)WindowConstructor::fStatic_Elasticity2D,0,0);

}
