//This cpp file describes the fButton_Analysis function (Tool Bar Menu from FSArl).

#include "WindowConstructor.h"

void WindowConstructor::fButton_Analysis(Fl_Group *group){

    png_image_A = new Fl_PNG_Image("../GraphicInterface/Icons/1x/Analysis_Icon.png");

    Button_Analysis = new Fl_Button(135, 35, 45, 45);
    Button_Analysis->box(FL_NO_BOX);
    Button_Analysis->image(png_image_A);
    Button_Analysis->callback((Fl_Callback*)fStatic_Button_Analysis);
    Button_Analysis->redraw();
    
    group->add(Button_Analysis);
}

