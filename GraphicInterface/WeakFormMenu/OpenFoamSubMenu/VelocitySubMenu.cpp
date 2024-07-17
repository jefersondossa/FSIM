#include "WindowConstructor.h"

Fl_Menu_Item WindowConstructor::menu_BC1_OpenFoam[] = {
        {"Set a BC", 0,  0, 0, 1, (uchar)FL_NORMAL_LABEL, 0, 14, 8},
        {0,0,0,0,0,0,0,0,0}
    };

void WindowConstructor::fU_OpenFoam(Fl_Group *group){

    U_OpenFoam_Menu = new Fl_Group(1160, 240, 375, 175);
    U_OpenFoam_Menu->box(FL_BORDER_BOX);
    U_OpenFoam_Menu->show();

    Fl_Box* o = new Fl_Box(1300, 260, 160, 25, "Velocity");
    o->labelsize(18);
    o->align(FL_ALIGN_CENTER);

    BC1_OpenFoam = new Fl_Choice(1320, 300, 205, 25, "Boundary Condition: ");
    BC1_OpenFoam->box(FL_BORDER_BOX);
    BC1_OpenFoam->down_box(FL_BORDER_BOX);
    BC1_OpenFoam->menu(menu_BC1_OpenFoam);

    this->BC1_OpenFoam->add("noSlip",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
    this->BC1_OpenFoam->add("Slip",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
    //this->BC1_OpenFoam->add("slip",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
    //this->BC1_OpenFoam->add("symmetry",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
    this->BC1_OpenFoam->add("empty",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
    this->BC1_OpenFoam->add("fixedValue",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,0);
    this->BC1_OpenFoam->add("zeroGradient",0,(Fl_Callback*)WindowConstructor::fStatic_fixedValue,0,1);  

    Values1 = new Fl_Box(1172, 340, 58, 25, "Vector");  
    x = new Fl_Value_Input(1255, 340, 60, 25, "X: ");   
    y = new Fl_Value_Input(1350, 340, 60, 25, "Y: ");   
    z = new Fl_Value_Input(1450, 340, 60, 25, "Z: "); 

    x->deactivate();
    y->deactivate();
    z->deactivate();
    Values1->deactivate();

    U_OpenFoam_Menu->end();

    group->add(U_OpenFoam_Menu);

}