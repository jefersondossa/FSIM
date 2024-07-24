#include "WindowConstructor.h"

Fl_Menu_Item WindowConstructor::menu_BC_nuTilda_OpenFoam[] = {
        {"Set a BC", 0,  0, 0, 1, (uchar)FL_NORMAL_LABEL, 0, 14, 8},
        {0,0,0,0,0,0,0,0,0}
    };

void WindowConstructor::fnuTilda_OpenFoam(Fl_Group *group){

    nuTilda_OpenFoam_Menu = new Fl_Group(1160, 240, 375, 175);
    nuTilda_OpenFoam_Menu->box(FL_BORDER_BOX);
    nuTilda_OpenFoam_Menu->resizable();
    nuTilda_OpenFoam_Menu->hide();

    Fl_Box* o = new Fl_Box(1300, 260, 160, 25, "Modified turbulent viscosity");
    o->labelsize(18);
    o->align(FL_ALIGN_CENTER);

    BC_nuTilda_OpenFoam = new Fl_Choice(1320, 300, 205, 25, "Boundary Condition: ");
    BC_nuTilda_OpenFoam->box(FL_BORDER_BOX);
    BC_nuTilda_OpenFoam->down_box(FL_BORDER_BOX);
    BC_nuTilda_OpenFoam->menu(menu_BC_U_OpenFoam);

    this->BC_nuTilda_OpenFoam->add("noSlip",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
    this->BC_nuTilda_OpenFoam->add("slip",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
    this->BC_nuTilda_OpenFoam->add("empty",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
    this->BC_nuTilda_OpenFoam->add("fixedValue",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);  
    this->BC_nuTilda_OpenFoam->add("zeroGradient",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,1);

    nuTilda = new Fl_Value_Input(1255, 340, 60, 25, "value: ");   

    nuTilda->deactivate();

    nuTilda_OpenFoam_Menu->end();

    group->add(nuTilda_OpenFoam_Menu);

}