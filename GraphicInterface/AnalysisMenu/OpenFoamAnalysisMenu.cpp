//This cpp file describes the fOFAnalysis function (Analysis Menu from FSArl).

#include "WindowConstructor.h"

Fl_Menu_Item WindowConstructor::menu_FluidSimulation[] = {
        {"Set a Simulation Type", 0,  0, 0, 1, (uchar)FL_NORMAL_LABEL, 0, 14, 8},
        {0,0,0,0,0,0,0,0,0}
    };

Fl_Menu_Item WindowConstructor::menu_FluidFlow[] = {
        {"Set a Flow Type", 0,  0, 0, 1, (uchar)FL_NORMAL_LABEL, 0, 14, 8},
        {0,0,0,0,0,0,0,0,0}
    };

void WindowConstructor::fOFAnalysis(Fl_Group *group){

    OpenFoam_Analysis_Menu = new Fl_Group(1160, 175, 375, 600);
    OpenFoam_Analysis_Menu->hide();

    /* ========================= FluidSimulation ========================= */

    FluidSimulation = new Fl_Choice(1325, 185, 185, 25, "Simulation Type:");
    FluidSimulation->box(FL_BORDER_BOX);
    FluidSimulation->down_box(FL_BORDER_BOX);
    FluidSimulation->menu(menu_FluidSimulation);

    this->FluidSimulation->add("Incompressible Fluid");

    /* =========================  FluidFlow ========================= */

    FluidFlow = new Fl_Choice(1325, 225, 185, 25, "Flow Type:");
    FluidFlow->box(FL_BORDER_BOX);
    FluidFlow->down_box(FL_BORDER_BOX);
    FluidFlow->menu(menu_FluidFlow);

    this->FluidFlow->add("Laminar");

    /*============================== viscosity ==============================*/

    OFViscosity = new Fl_Value_Input(1325, 265, 95, 25, "Viscosity: ");

    /*============================== ControlDict ===================================*/

    endTime = new Fl_Value_Input(1325, 305, 100, 25, "endTime: ");
    deltaT = new Fl_Value_Input(1325, 345, 100, 25, "deltaT: ");
    writeInterval = new Fl_Value_Input(1325, 385, 100, 25, "writeInterval: ");
    writePrecision = new Fl_Value_Input(1325, 425, 100, 25, "writePrecision: ");
    timePrecision = new Fl_Value_Input(1325, 465, 100, 25, "timePrecision: ");

    endTime->value(0);
    deltaT->value(0.1);
    writeInterval->value(1);
    writePrecision->value(6);
    timePrecision->value(6);

    Apply_Analysis = new Fl_Button(1430, 505, 80, 25, "Apply");
    Apply_Analysis->callback((Fl_Callback*)fStatic_Apply_Analysis);

    OpenFoam_Analysis_Menu->end();

    group->add(OpenFoam_Analysis_Menu);
}