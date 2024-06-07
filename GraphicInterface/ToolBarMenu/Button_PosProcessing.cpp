//This cpp file describes the fButton_PosProcressing function (Tool Bar Menu from FSArl).

#include "WindowConstructor.h"

void WindowConstructor::fButton_PosProcressing(Fl_Group *group){

    png_image_PP = new Fl_PNG_Image("../GraphicInterface/Icons/1x/ParaView_Icon.png");

    Button_PosProcressing = new Fl_Button(200, 35, 45, 45);
    Button_PosProcressing->box(FL_NO_BOX);
    Button_PosProcressing->image(png_image_PP);
    Button_PosProcressing->callback((Fl_Callback*)fStatic_Button_PPospricessing);
    Button_PosProcressing->show();

    group->add(Button_PosProcressing);
    
}