//This cpp file describes the fL2Projection function (L2Projection Submenu from Weak Form Menu).

#include "WindowConstructor.h"

Fl_Menu_Item WindowConstructor::menu_BoundaryCondition[] = {
        {"Set a Boundary Condition", 0,  0, 0, 1, (uchar)FL_NORMAL_LABEL, 0, 14, 8},
        {0,0,0,0,0,0,0,0,0}
    };

void WindowConstructor::fL2Projection(Fl_Group *group){

  L2Projection_Menu = new Fl_Group(1175, 175, 360, 175);
  L2Projection_Menu->hide();

  BoundaryCondition = new Fl_Choice(1320, 190, 205, 25, "Boundary Condition:");
  BoundaryCondition->box(FL_BORDER_BOX);
  BoundaryCondition->down_box(FL_BORDER_BOX);
  BoundaryCondition->menu(menu_BoundaryCondition);

  Values = new Fl_Box(1172, 240, 58, 25, "Values: ");

  x = new Fl_Value_Input(1255, 240, 55, 25, "X: ");

  y = new Fl_Value_Input(1350, 240, 60, 25, "Y: ");

  z = new Fl_Value_Input(1450, 240, 60, 25, "Z: ");

  //Adds items to Boundary Condition choice menu:
  this->BoundaryCondition->add("Dirichlet");
  this->BoundaryCondition->add("Neumann");
  this->BoundaryCondition->add("Directional Homogeneous Dirichlet");
  this->BoundaryCondition->add("Directional Non-Homogeneous Dirichlet");

 /*============================== Button_Apply_WeakForm ==============================*/

  Button_Apply_WeakForm = new Fl_Button(1410, 270, 80, 25, "Apply");
  Button_Apply_WeakForm->callback((Fl_Callback*)fStatic_Apply_WF);

  L2Projection_Menu->end();

  group->add(L2Projection_Menu);

   this->WeakForm->add("L2 Projection",0,(Fl_Callback*)WindowConstructor::fStatic_L2Projection,0,0);

}