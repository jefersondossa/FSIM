/* This cpp file describes the following Script Display functions:

1.  fButton_SaveScript_cb
2.  fButton_TextSave_cb
3.  fButton_TextCancel_cb

*/

#include "WindowConstructor.h"

/* ========================= fButton_SaveScript_cb() ========================= */

void WindowConstructor::fButton_SaveScript_cb(){

    this->Script_Popup_Menu->show();

}

/* ========================= fButton_TextSave_cb() ========================= */

void WindowConstructor::fButton_TextSave_cb(){

//**********************************************************************************************
/*This function is separated in 3 parts:

  1. Reads the txt input file name and creates the name.txt;
  2. Reads the txt input file name and creates the name.h;

//**********************************************************************************************

/* ========================= PART 1. ========================= */

  string name = Script_Text_Name->value();
  string name_txt = name + ".txt";
  string name_path_txt = "../GraphicInterface/ScriptFiles/" + name_txt;

  fstream new_script;

  new_script.open(name_path_txt,ios::out); //Allows to edit new_script.
  ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in);//Allow to read ScriptMemory.

  string Str;
  string Line;

  if(ScriptMemory_txt.is_open()){ //Opens ScriptMemory.
      Str.clear();
      while(getline(ScriptMemory_txt, Line)){
        Str += Line;
        Str += "\n";
      }
      ScriptMemory_txt.close();//Closes ScriptMemory.
      if(new_script.is_open()){
        new_script << Str; //Copy text to new_script

        new_script.close();
      }
    }

  ScriptMemory_txt.close(); // Fecha ScriptMemory.

this->Script_Popup_Menu->hide();
}

/* ========================= fButton_TextCancel_cb() ========================= */

void WindowConstructor::fButton_TextCancel_cb(){
    this->Script_Popup_Menu->hide();
}