//This cpp file describes the fScript_Popup_Menu function (Save Script Name Popup from FSArl).

#include "WindowConstructor.h"

void WindowConstructor::fScript_Popup_Menu(Fl_Double_Window *window){

    Script_Popup_Menu = new Fl_Group(10, 580, 425, 135);
    Script_Popup_Menu->box(FL_SHADOW_BOX);
    Script_Popup_Menu->color(FL_LIGHT1);
    Script_Popup_Menu->hide();

    Script_Popup_Text = new Fl_Box(50, 635, 355, 30, "The Script File will be saved at ScriptFiles folder as .txt.");

    Script_Text_Name = new Fl_Input(75, 601, 337, 28, "Name: ");

    Fl_Button* Button_TextSave = new Fl_Button(30, 675, 70, 25, "Save");
    Button_TextSave->callback((Fl_Callback*)fStatic_Button_TextSave);

    Fl_Button* Button_TextCancel = new Fl_Button(343, 675, 70, 25, "Cancel");
    Button_TextCancel->callback((Fl_Callback*)fStatic_Button_TextCancel);


    Script_Popup_Menu->end();

    window->add(Script_Popup_Menu);
}