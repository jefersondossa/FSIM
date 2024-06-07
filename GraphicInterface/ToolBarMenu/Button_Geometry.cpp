//This cpp file describes the fButton_Geometry function (Tool Bar Menu from FSArl).

#include "WindowConstructor.h"

void WindowConstructor::fButton_Geometry(Fl_Group *group){

    png_image_G = new Fl_PNG_Image("../GraphicInterface/Icons/1x/Geometry_Icon.png");

    Button_Geometry = new Fl_Button(10, 35, 45, 45);
    Button_Geometry->box(FL_NO_BOX);
    Button_Geometry->image(png_image_G);
    Button_Geometry->callback((Fl_Callback*)fStatic_Button_Geometry);
    Button_Geometry->redraw();

    group->add(Button_Geometry);
}