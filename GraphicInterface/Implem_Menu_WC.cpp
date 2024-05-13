#include "WindowConstructor.h"
#include <FL/Fl_File_Chooser.H>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

string str;
string file_path;
fstream ScriptMemory_3;
fstream ScriptTxt;
ifstream file;
Fl_Text_Buffer *buff3 = new Fl_Text_Buffer;


void WindowConstructor::File_Open_f(){

  //Clears ScriptMemory.txt text
  ScriptMemory_3.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);
  if(ScriptMemory_3.is_open()){
    ScriptMemory_3.close();
  }

  // 1. Set User to choose a txt file path:

  Fl_File_Chooser txt_file("../GraphicInterface/ScriptFiles/", "*.txt", Fl_File_Chooser::SINGLE, "Select a .txt file to open");

  txt_file.show();
  while (txt_file.shown()){
    Fl::wait();
  }

  if (txt_file.value() != nullptr & txt_file.count() > 0){
    file_path = txt_file.value();//Saves the path of the file that the user has choosen.

    //Clears Script Display
    buff3->text("");
    scriptdisplay->buffer(buff3);
  }
  else{
    cout << "User canceled the file chooser dialog.";
    
    
  }

  txt_file.hide();

  // 2. Copy text from the file choosen to string str.

  //Since we don´t want to accumlate text inside the str, lets clear it:
  str = "";

  ScriptTxt.open(file_path,ios::in);

  if(ScriptTxt.is_open()){ 
    string line;
    while(getline(ScriptTxt, line)){
      str += line;
      str += "\n"; //Saves the text of the file choosen in str string (str = ScriptTxt)
    }
    ScriptTxt.close();
  }

  // 3. Copy text from the file choosen to ScriptMemory_3

  ScriptMemory_3.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);

  if(ScriptMemory_3.is_open()){
    ScriptMemory_3 << str;
    ScriptMemory_3.close(); // (ScriptMemory_3 = ScriptTxt).
  }

  // 4. Search for geo. and msh. paths.

  ScriptTxt.open(file_path,ios::in);

  if(ScriptTxt.is_open()){

    string line;

    while(getline(ScriptTxt,line)){
      size_t pos = line.find(": ");
      if(pos != string::npos){
      string key = line.substr(0,pos);
        string path = line.substr(pos + 2);
        if(key == "geo path"){
          this->geo_path->value(path.c_str());
        }
        
        else if(key == "msh path"){
          this->msh_path->value(path.c_str());
        }
      }
    
  }
  ScriptTxt.close();
  buff3->text();
  scriptdisplay->buffer(buff3);
  this->Save_paths->do_callback();
}
}

void WindowConstructor::File_New_f(){
  
  //Clears ScriptMemory
  ScriptMemory_3.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);
  if(ScriptMemory_3.is_open()){
    ScriptMemory_3.close();
  }

  //Clears Container_1 (geo and msh paths)
  this->geo_path->value("");
  this->msh_path->value("");

  //Clears str
  str = "";

  //buff3->remove(0,buff3->length());
  buff3->text("");
  scriptdisplay->buffer(buff3);

  //Lock Weak Form (Container_2)
  this->WF_Icon->deactivate();
  this->Container_2->hide();
  this->Container_5->hide();
  this->Container_1->show();

}


