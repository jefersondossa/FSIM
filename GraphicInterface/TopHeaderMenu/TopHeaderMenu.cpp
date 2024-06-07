//This cpp file describes the Top Header Menu function (Top Header Menu from FSArl).

#include "WindowConstructor.h"

Fl_Menu_Item WindowConstructor::Topheader_Menu[] = 
{
    {"  File  ", 0,  0, 0, 192, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
    {"Open...    ", 0,  (Fl_Callback*)WindowConstructor::fStatic_Open_File, 0, 128, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
    {"New...      ", 0,  (Fl_Callback*)WindowConstructor::fStatic_New_File, 0, 0, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
    {0,0,0,0,0,0,0,0,0},
    {"  Options  ", 0x4006f,  0, 0, 128, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
    {"  Help  ", 0,  0, 0, 64, (uchar)FL_NORMAL_LABEL, 0, 14, 0},
    {0,0,0,0,0,0,0,0,0}
};

void WindowConstructor::fTopHeader_Menu(Fl_Double_Window *window){


    //Creates the Top Header Menu
    TopHeader_Menu = new Fl_Menu_Bar(0, 0, 1535, 25);
    TopHeader_Menu->box(FL_FLAT_BOX);
    TopHeader_Menu->color(FL_WHITE);
    TopHeader_Menu->menu(Topheader_Menu);

    window->add(TopHeader_Menu);

}

