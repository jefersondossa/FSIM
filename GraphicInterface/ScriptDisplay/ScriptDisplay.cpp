//This cpp file describes the fScript_Menu function (Script Display from FSArl).

#include "WindowConstructor.h"

void WindowConstructor::fScript_Menu(Fl_Double_Window *window){

    Script_Menu = new Fl_Group(0, 79, 445, 705);
    Script_Menu->box(FL_BORDER_BOX);
    Script_Menu->color(FL_LIGHT1);
    Script_Menu->labelsize(18);

    ScriptDisplay = new Fl_Text_Display(0, 108, 445, 617, "Script");
    ScriptDisplay->box(FL_BORDER_BOX);
    ScriptDisplay->labelsize(18);

    Button_SaveScript = new Fl_Button(335, 738, 95, 32, "Save Script");
    Button_SaveScript->callback((Fl_Callback*)fStatic_Button_SaveScript);

    Script_Menu->end();

    window->add(Script_Menu);
}