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

    Analysis_Menu = new Fl_Group(1160, 75, 375, 340, "Analysis");
    Analysis_Menu->color(FL_LIGHT1);
    Analysis_Menu->labelsize(18);
    Analysis_Menu->hide();

    AnalysisType = new Fl_Choice(1325, 100, 185, 25, "Analysis Type: ");
    AnalysisType->box(FL_BORDER_BOX);
    AnalysisType->down_box(FL_BORDER_BOX);
    AnalysisType->menu(menu_Analysis);

    //Adds items to Analysis Type choice menu:
    this->AnalysisType->add("Linear",0,(Fl_Callback*)WindowConstructor::fStatic_Linear,0,0);
    this->AnalysisType->add("NonLinear",0,(Fl_Callback*)WindowConstructor::fStatic_NonLinear,0,0);
    this->AnalysisType->add("Incremental",0,(Fl_Callback*)WindowConstructor::fStatic_Increm_Transient,0,0);
    this->AnalysisType->add("Transient",0,(Fl_Callback*)WindowConstructor::fStatic_Increm_Transient,0,0);

    SolverType = new Fl_Choice(1325, 140, 185, 25, "Solver Type: ");
    SolverType->box(FL_BORDER_BOX);
    SolverType->down_box(FL_BORDER_BOX);
    SolverType->menu(menu_SolverType);

    //Adds items to Solver Type choice menu:
    this->SolverType->add("EMumps");
    this->SolverType->add("EIterative");
    this->SolverType->add("EUmfpack");
    this->SolverType->add("ECholmod");
    this->SolverType->add("EKLU");
    this->SolverType->add("ESPQR");
    this->SolverType->add("ELU");
    this->SolverType->add("ELLt");
    this->SolverType->add("ELDLt");
    this->SolverType->add("EQR");
    this->SolverType->add("ECG");
    this->SolverType->add("EBiCGStab");
    this->SolverType->add("ELSCG");

    fNonLinear(Analysis_Menu);
    fIncrem_Transient(Analysis_Menu);

    Apply_Analysis = new Fl_Button(1430, 350, 80, 25, "Apply");
    Apply_Analysis->callback((Fl_Callback*)fStatic_Apply_Analysis);

    Analysis_Menu->end();

    window->add(Analysis_Menu);
}
