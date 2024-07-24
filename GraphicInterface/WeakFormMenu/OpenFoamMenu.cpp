//This cpp file describes the fL2Projection function (L2Projection Submenu from Weak Form Menu).

#include "WindowConstructor.h"

Fl_Menu_Item WindowConstructor::system_OpenFoam_Menu[] = {

    {"  Parameters  ", 0,  fStatic_Parameter_OF, 0, 128, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
    {0,0,0,0,0,0,0,0,0}
  };

Fl_Menu_Item WindowConstructor::parameter_OpenFoam_Menu[] = 
{
    {" U ", 0,  fStatic_U_OF, 0, 128, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
    {" p ", 0,  fStatic_p_OF, 0, 128, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
    {" k ", 0,  fStatic_k_OF, 0, 128, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
    {" epsilon ", 0,  fStatic_epsilon_OF, 0, 128, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
    {" nuTilda ", 0,  fStatic_nuTilda_OF, 0, 128, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
    {" nut ", 0,  fStatic_nut_OF, 0, 128, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
    {" omega ", 0,  fStatic_omega_OF, 0, 128, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
    {0,0,0,0,0,0,0,0,0}
};

void WindowConstructor::fOpemFoam(Fl_Group *group){

  OpenFoam_Menu = new Fl_Group(1160, 175, 375, 500);
  OpenFoam_Menu->hide();

  /*============================== System Menu ==============================*/
    
    System_OpenFoam_Menu = new Fl_Menu_Bar(1160, 190, 375, 25);
    System_OpenFoam_Menu->box(FL_BORDER_BOX);
    System_OpenFoam_Menu->color(FL_WHITE);
    System_OpenFoam_Menu->menu(system_OpenFoam_Menu);

    Parameter_OpenFoam_Menu = new Fl_Menu_Bar(1160, 215, 375, 25);
    Parameter_OpenFoam_Menu->box(FL_BORDER_BOX);
    Parameter_OpenFoam_Menu->color(FL_WHITE);
    Parameter_OpenFoam_Menu->menu(parameter_OpenFoam_Menu);

    fU_OpenFoam(OpenFoam_Menu);
    fp_OpenFoam(OpenFoam_Menu);
    fk_OpenFoam(OpenFoam_Menu);
    fepsilon_OpenFoam(OpenFoam_Menu);
    fnuTilda_OpenFoam(OpenFoam_Menu);
    fnut_OpenFoam(OpenFoam_Menu);
    fomega_OpenFoam(OpenFoam_Menu);

  /*============================== Button_Apply_WeakForm ==============================*/

  Button_Apply_WeakForm = new Fl_Button(1410, 435, 80, 25, "Apply");
  Button_Apply_WeakForm->callback((Fl_Callback*)fStatic_Apply_WF);

  OpenFoam_Menu->end();

  group->add(OpenFoam_Menu);

  this->WeakForm->add("OpenFoam",0,(Fl_Callback*)WindowConstructor::fStatic_OpenFoam,0,0);

}