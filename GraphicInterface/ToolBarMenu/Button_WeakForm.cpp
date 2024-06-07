//This cpp file describes the fButton_WeakForm function (Tool Bar Menu from FSArl).

#include "WindowConstructor.h"

void WindowConstructor::fButton_WeakForm(Fl_Group *group){

    png_image_WF_act = new Fl_PNG_Image("../GraphicInterface/Icons/1x/WeakForm_Icon.png");
    png_image_WF_desact = new Fl_PNG_Image("../GraphicInterface/Icons/1x/WeakFormInactive_Icon.png");

    Button_WeakForm = new Fl_Button(70, 35, 45, 45);
    Button_WeakForm->box(FL_NO_BOX);
    Button_WeakForm->image(png_image_WF_act);
    Button_WeakForm->deimage(png_image_WF_desact);
    Button_WeakForm->callback((Fl_Callback*)fStatic_Button_WeakForm);
    Button_WeakForm->deactivate();
    Button_WeakForm->show();

    group->add(Button_WeakForm);
}