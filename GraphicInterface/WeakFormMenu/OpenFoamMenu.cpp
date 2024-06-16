//This cpp file describes the fL2Projection function (L2Projection Submenu from Weak Form Menu).

#include "WindowConstructor.h"

Fl_Menu_Item WindowConstructor::menu_FluidSimulation[] = {
        {"Set a Simulation Type", 0,  0, 0, 1, (uchar)FL_NORMAL_LABEL, 0, 14, 8},
        {0,0,0,0,0,0,0,0,0}
    };

Fl_Menu_Item WindowConstructor::menu_FluidFlow[] = {
        {"Set a Flow Type", 0,  0, 0, 1, (uchar)FL_NORMAL_LABEL, 0, 14, 8},
        {0,0,0,0,0,0,0,0,0}
    };

Fl_Menu_Item WindowConstructor::menu_FluidVelocity[] = {
        {"Set a Velocity Type", 0,  0, 0, 1, (uchar)FL_NORMAL_LABEL, 0, 14, 8},
        {0,0,0,0,0,0,0,0,0}
    };

Fl_Menu_Item WindowConstructor::menu_FluidPressure[] = {
        {"Set a Pressure Type", 0,  0, 0, 1, (uchar)FL_NORMAL_LABEL, 0, 14, 8},
        {0,0,0,0,0,0,0,0,0}
    };

void WindowConstructor::fOpemFoam(Fl_Group *group){

  OpenFoam_Menu = new Fl_Group(1175, 175, 460, 420);
  OpenFoam_Menu->hide();

  FluidSimulation = new Fl_Choice(1320, 190, 205, 25, "Simulation Type:");
  FluidSimulation->box(FL_BORDER_BOX);
  FluidSimulation->down_box(FL_BORDER_BOX);
  FluidSimulation->menu(menu_FluidSimulation);

  this->FluidSimulation->add("Incompressible Fluid");

  FluidFlow = new Fl_Choice(1320, 230, 205, 25, "Flow Type:");
  FluidFlow->box(FL_BORDER_BOX);
  FluidFlow->down_box(FL_BORDER_BOX);
  FluidFlow->menu(menu_FluidFlow);

  this->FluidFlow->add("Laminar");

  FluidVelocity = new Fl_Choice(1320, 270, 205, 25, "Velocity Type:");
  FluidVelocity->box(FL_BORDER_BOX);
  FluidVelocity->down_box(FL_BORDER_BOX);
  FluidVelocity->menu(menu_FluidVelocity);

  this->FluidVelocity->add("zeroGradient",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
  this->FluidVelocity->add("noSlip",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
  this->FluidVelocity->add("slip",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
  this->FluidVelocity->add("symmetry",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
  this->FluidVelocity->add("empty",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
  this->FluidVelocity->add("fixedValue",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);

  Values = new Fl_Box(1172, 310, 58, 25, "U vector");

  x = new Fl_Value_Input(1255, 310, 60, 25, "X: ");

  y = new Fl_Value_Input(1350, 310, 60, 25, "Y: ");

  z = new Fl_Value_Input(1450, 310, 60, 25, "Z: ");

  x->deactivate();
  y->deactivate();
  z->deactivate();
  Values->deactivate();

  FluidPressure = new Fl_Choice(1320, 350, 205, 25, "Pressure Type:");
  FluidPressure->box(FL_BORDER_BOX);
  FluidPressure->down_box(FL_BORDER_BOX);
  FluidPressure->menu(menu_FluidVelocity);

  this->FluidPressure->add("zeroGradient");
  this->FluidPressure->add("noSlip");
  this->FluidPressure->add("slip");
  this->FluidPressure->add("symmetry");
  this->FluidPressure->add("empty");
  this->FluidPressure->add("fixedValue");

  viscosity = new Fl_Value_Input(1320, 390, 95, 25, "Viscosity: ");

  /*============================== Button_Apply_WeakForm ==============================*/

  Button_Apply_WeakForm = new Fl_Button(1410, 430, 80, 25, "Apply");
  Button_Apply_WeakForm->callback((Fl_Callback*)fStatic_Apply_WF);

  OpenFoam_Menu->end();

  group->add(OpenFoam_Menu);

  this->WeakForm->add("OpenFoam",0,(Fl_Callback*)WindowConstructor::fStatic_OpenFoam,0,0);

}