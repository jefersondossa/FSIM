#include "WindowConstructor.h"

void WindowConstructor::f_Open(){

    Window->show();

    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);
    if(ScriptMemory_txt.is_open()){
      ScriptMemory_txt.close();
    }
}