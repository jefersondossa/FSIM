/* This cpp file describes the following Top Header Menu functions:

1.  fOpen_File_cb;
2.  fNew_File_cb.
3.  fGenerate_cb;

*/

#include "WindowConstructor.h"

/* ========================= fOpen_File_cb() ========================= */

void WindowConstructor::fOpen_File_cb(){

  string Str;

    //Set the first item as the default even after change the file open
  this->PhysicalGroups->value(0);

  //Lock Weak Form (Container_2)
  this->Button_WeakForm->deactivate();
  this->WF_Menu->hide();

  //Clears ScriptMemory.txt text
  ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);
  if(ScriptMemory_txt.is_open()){
    ScriptMemory_txt.close();
  }

  // 1. Set User to choose a txt file path:

  Fl_File_Chooser txt_file("../GraphicInterface/ScriptFiles/", "*.txt", Fl_File_Chooser::SINGLE, "Select a .txt file to open");

  txt_file.show();
  while (txt_file.shown()){
    Fl::wait();
  }

    string file_path; //Local Object that stores the path of the file that the user has choosen.

  if (txt_file.value() != nullptr & txt_file.count() > 0){
    file_path = txt_file.value();

    //Clears Script Display
    Buffer.text("");
    ScriptDisplay->buffer(Buffer);
  }

  else{
    cout << "User canceled the file chooser dialog.";
    
  }

  txt_file.hide();

  // 2. Copy text from the file choosen to string str.

  //Since we don´t want to accumlate text inside the Str, lets clear it:
  Str = "";

    fstream ScriptTxt;
  ScriptTxt.open(file_path,ios::in);

  if(ScriptTxt.is_open()){ 
    string line;
    while(getline(ScriptTxt, line)){
      Str += line;
      Str += "\n"; //Saves the text of the file choosen in str string (str = ScriptTxt)
    }
    ScriptTxt.close();
  }

  // 3. Copy text from the file choosen to ScriptMemory_3

  ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);

  if(ScriptMemory_txt.is_open()){
    ScriptMemory_txt << Str;
    ScriptMemory_txt.close(); // (ScriptMemory_3 = ScriptTxt).
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
          this->geo->value(path.c_str());
        }
        
        else if(key == "msh path"){
          this->msh->value(path.c_str());
        }
      }
    
  }

  ScriptTxt.close();
  Buffer.text(); //Clears the Buffer
  ScriptDisplay->buffer(Buffer); //Clears the
  this->Apply_geometry->do_callback();
  Geometry_Menu->show();
  WF_Menu->hide();
  Analysis_Menu->hide();
  PProcessing_Menu->hide();
  zeroFiles.clear();
  systemFiles.clear();
  constantFiles.clear();
  
}

}

/* ========================= fNew_File_cb() ========================= */

void WindowConstructor::fNew_File_cb(){

//Clears OpenGL window
glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

//Clears ScriptMemory
ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);
if(ScriptMemory_txt.is_open()){
  ScriptMemory_txt.close();
}

//Clears Container_1 (geo and msh paths)
this->geo->value("");
this->msh->value("");

//Clears Dimension
this->dimension->value(0);

//Clears str
string Str;
ostringstream Oss;

//buff3->remove(0,buff3->length());
Buffer.text("");
ScriptDisplay->buffer(Buffer);

//Lock Weak Form (Container_2)
this->Button_WeakForm->deactivate();
this->WF_Menu->hide();
this->Analysis_Menu->hide();
this->PProcessing_Menu->hide();
this->Geometry_Menu->show(); 

//Set the first item as the default even after change the file open
this->PhysicalGroups->value(0);

zeroFiles.clear();
systemFiles.clear();
constantFiles.clear();
}

/* ========================= fGenerateHeader_cb() ========================= */

void WindowConstructor::fGenerateHeader_cb(){

  fHeaderFile();
}

/* ========================= fGenerateFoam_cb() ========================= */

void WindowConstructor::fGenerateFoam_cb(){

  fFoamFile();
}