//This cpp file describes the fAnalysis function (Analysis Menu from FSArl).

#include "WindowConstructor.h"

Fl_Menu_Item WindowConstructor::menu_Analysis[] = {
        {"Set a Analysis Type", 0,  0, 0, 1, (uchar)FL_NORMAL_LABEL, 0, 14, 8},
        {0,0,0,0,0,0,0,0,0}
    };

Fl_Menu_Item WindowConstructor::menu_SolverType[] = {
        {"Set a Solver Type", 0,  0, 0, 1, (uchar)FL_NORMAL_LABEL, 0, 14, 8},
        {0,0,0,0,0,0,0,0,0}
    };

void WindowConstructor::fAnalysis(Fl_Double_Window *window){

    Analysis_Menu = new Fl_Group(1160, 75, 375, 500, "Analysis");
    Analysis_Menu->color(FL_LIGHT1);
    Analysis_Menu->labelsize(18);
    Analysis_Menu->hide();

    AnalysisType = new Fl_Choice(1325, 100, 185, 25, "Analysis Type: ");
    AnalysisType->box(FL_BORDER_BOX);
    AnalysisType->down_box(FL_BORDER_BOX);
    AnalysisType->menu(menu_Analysis);

    //Adds items to Analysis Type choice menu:
    AnalysisType->add("Linear",0,(Fl_Callback*)WindowConstructor::fStatic_Linear,0,0);
    AnalysisType->add("NonLinear",0,(Fl_Callback*)WindowConstructor::fStatic_NonLinear,0,0);
    AnalysisType->add("Incremental",0,(Fl_Callback*)WindowConstructor::fStatic_Increm_Transient,0,0);
    AnalysisType->add("Transient",0,(Fl_Callback*)WindowConstructor::fStatic_Increm_Transient,0,0);
    AnalysisType->add("OpenFoam Analysis",0,(Fl_Callback*)WindowConstructor::fStatic_OF_Analysis,0,0);

    SolverType = new Fl_Choice(1325, 140, 185, 25, "Solver Type: ");
    SolverType->box(FL_BORDER_BOX);
    SolverType->down_box(FL_BORDER_BOX);
    SolverType->menu(menu_SolverType);
    SolverType->deactivate();

    //Adds items to Solver Type choice menu:
    SolverType->add("EMumps");
    SolverType->add("EIterative");
    SolverType->add("EUmfpack");
    SolverType->add("ECholmod");
    SolverType->add("EKLU");
    SolverType->add("ESPQR");
    SolverType->add("ELU");
    SolverType->add("ELLt");
    SolverType->add("ELDLt");
    SolverType->add("EQR");
    SolverType->add("ECG");
    SolverType->add("EBiCGStab");
    SolverType->add("ELSCG");

    fNonLinear(Analysis_Menu);
    fIncrem_Transient(Analysis_Menu);
    fOFAnalysis(Analysis_Menu);

    Analysis_Menu->end();

    window->add(Analysis_Menu);
}
