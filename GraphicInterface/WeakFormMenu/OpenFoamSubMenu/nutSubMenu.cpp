#include "WindowConstructor.h"

Fl_Menu_Item WindowConstructor::menu_BC_nut_OpenFoam[] = {
        {"Set a BC", 0,  0, 0, 1, (uchar)FL_NORMAL_LABEL, 0, 14, 8},
        {0,0,0,0,0,0,0,0,0}
    };

void WindowConstructor::fnut_OpenFoam(Fl_Group *group){

    nut_OpenFoam_Menu = new Fl_Group(1160, 240, 375, 175);
    nut_OpenFoam_Menu->box(FL_BORDER_BOX);
    nut_OpenFoam_Menu->resizable();
    nut_OpenFoam_Menu->hide();

    Fl_Box* o = new Fl_Box(1300, 260, 160, 25, "Turbulent kinematic viscosity");
    o->labelsize(18);
    o->align(FL_ALIGN_CENTER);

    BC_nut_OpenFoam = new Fl_Choice(1320, 300, 205, 25, "Boundary Condition: ");
    BC_nut_OpenFoam->box(FL_BORDER_BOX);
    BC_nut_OpenFoam->down_box(FL_BORDER_BOX);
    BC_nut_OpenFoam->menu(menu_BC_U_OpenFoam);

    this->BC_nut_OpenFoam->add("noSlip",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
    this->BC_nut_OpenFoam->add("slip",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
    this->BC_nut_OpenFoam->add("empty",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
    this->BC_nut_OpenFoam->add("fixedValue",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);  
    this->BC_nut_OpenFoam->add("zeroGradient",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,1);

    nut = new Fl_Value_Input(1255, 340, 60, 25, "value: ");   

    nut->deactivate();

    nut_OpenFoam_Menu->end();

    group->add(nut_OpenFoam_Menu);

}