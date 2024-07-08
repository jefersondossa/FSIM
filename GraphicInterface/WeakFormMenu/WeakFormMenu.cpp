/*This cpp file describes the fWF functions (Weak Form Menu from FSArl). It contains:

    01. WF_Menu;
    02. PhysicalGroups;
    03. WeakForm;
    04. Button_Apply_WeakForm.
*/

#include "WindowConstructor.h"



Fl_Menu_Item WindowConstructor::menu_PhysicalGroups[] = {
        {"Set a Physical Group", 0,  0, 0, 1, (uchar)FL_NORMAL_LABEL, 0, 14, 8},
        {0,0,0,0,0,0,0,0,0}
    };

Fl_Menu_Item WindowConstructor::menu_WeakForm[] = {
        {"Set a Weak Form", 0,  0, 0, 1, (uchar)FL_NORMAL_LABEL, 0, 14, 8},
        {0,0,0,0,0,0,0,0,0}
    };


void WindowConstructor::fWF(Fl_Double_Window *window){

/*============================== WF_Menu ==============================*/

    WF_Menu = new Fl_Group(1160, 75, 375, 700, "Weak Form");
    WF_Menu->color(FL_LIGHT1);
    WF_Menu->labelsize(18);
    WF_Menu->hide();

    /*============================== PhysicalGroups ==============================*/

    PhysicalGroups = new Fl_Choice(1305, 100, 185, 25, "Physical Groups: ");
    PhysicalGroups->box(FL_BORDER_BOX);
    PhysicalGroups->down_box(FL_BORDER_BOX);
    PhysicalGroups->menu(menu_PhysicalGroups);

    /*============================== WeakForm ==============================*/

    WeakForm = new Fl_Choice(1305, 145, 185, 25, "Weak Form: ");
    WeakForm->box(FL_BORDER_BOX);
    WeakForm->down_box(FL_BORDER_BOX);
    WeakForm->deactivate();
    WeakForm->menu(menu_WeakForm);

    fElasticity2D(WF_Menu);
    fElasticity3D(WF_Menu);
    fElasticityPositional2D(WF_Menu);
    fElasticTruss(WF_Menu);
    fL2Projection(WF_Menu);
    fNavierStokes(WF_Menu);
    fPoisson(WF_Menu);
    fPositionalTruss(WF_Menu);
    fStokes(WF_Menu);
    fOpemFoam(WF_Menu);

    WF_Menu->end();

    window->add(WF_Menu);

}