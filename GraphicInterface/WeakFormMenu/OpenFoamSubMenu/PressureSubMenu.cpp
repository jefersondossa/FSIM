#include "WindowConstructor.h"

Fl_Menu_Item WindowConstructor::menu_BC_p_OpenFoam[] = {
        {"Set a BC", 0,  0, 0, 1, (uchar)FL_NORMAL_LABEL, 0, 14, 8},
        {0,0,0,0,0,0,0,0,0}
    };

void WindowConstructor::fp_OpenFoam(Fl_Group *group){

    p_OpenFoam_Menu = new Fl_Group(1160, 240, 375, 175);
    p_OpenFoam_Menu->box(FL_BORDER_BOX);
    p_OpenFoam_Menu->resizable();
    p_OpenFoam_Menu->hide();

    Fl_Box* o = new Fl_Box(1300, 260, 160, 25, "Pressure");
    o->labelsize(18);
    o->align(FL_ALIGN_CENTER);

    BC_p_OpenFoam = new Fl_Choice(1320, 300, 205, 25, "Boundary Condition: ");
    BC_p_OpenFoam->box(FL_BORDER_BOX);
    BC_p_OpenFoam->down_box(FL_BORDER_BOX);
    BC_p_OpenFoam->menu(menu_BC_U_OpenFoam);

    this->BC_p_OpenFoam->add("noSlip",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
    this->BC_p_OpenFoam->add("slip",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
    //this->BC_U_OpenFoam->add("slip",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
    //this->BC_U_OpenFoam->add("symmetry",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
    this->BC_p_OpenFoam->add("empty",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
    this->BC_p_OpenFoam->add("fixedValue",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);  
    this->BC_p_OpenFoam->add("zeroGradient",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,1);

    pressure = new Fl_Value_Input(1255, 340, 60, 25, "value: ");   

    pressure->deactivate();

    p_OpenFoam_Menu->end();

    group->add(p_OpenFoam_Menu);
}