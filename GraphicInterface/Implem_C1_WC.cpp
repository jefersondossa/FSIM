/*
This file describes the Container_2 and Container_3 created in WindowConstructor.h. The Window is divided in Containers, being:

-> Container_1: Contains the geo. and msh. files input and output (File Path (Container_3)). Besides, reads and interprets the geo. file,
granting the user, in the Container_2, to choose all analysis type and variant. 

-> Container_2: Allows the user to input the desired variants values. This values can be verified at Script (Container_3).

-> Container_3: It's a text display which allows the user to see the entries that he has placed.

-> Container_4: Allerts erros for the user.

-> Container_5: Contains the Analysis.

*/

#include "WindowConstructor.h"
#include "OpenGL_WC.h"
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Native_File_Chooser.H>
#include <string>
#include <vector>
#include <iostream>

//Scanner Libraries:
#include <fstream>
#include <regex>

using namespace std;
//Stores variables values
//The Fl_text_Buffer stores text that will be displayed at FL_Text_Display

Fl_Text_Buffer *buff2 = new Fl_Text_Buffer; //buff2 for files path;

std::ostringstream oss; /*Since the variables are contained in objects and, since Fl_buffer function only accepts "text", it's necessary to save the entire
text ("text" << objects) inside the oss object. */

//Define a txt file where will be all Script Memory information
fstream ScriptMemory_2; //fstream library: out to write mode (overwrites), app (append), in (read)

//Create the variable line and the reading pattern for the .geo file:
std::string line;
std::regex pattern("(Physical (Curve|Surface|Volume))\\(\"(\\w+)\", (\\d+)\\)");

//Parameters for repetition count (curve, surface and volume):
int curveCount = 0;
int surfaceCount = 0;
int volumeCount = 0;
int PhysicalGroupCount;

//Parameter for the Dimension:
int d;

//Structure in vector format named "names" to store the name = match[3] values:
vector<string> names;
string name;

//Structure in vector format named "matid" to store the matid = match[4] values:

int matid;

//________________________________________________________________________________

//Opens Window:
void WindowConstructor::Open_f() {
  this->Window->show();
  ScriptMemory_2.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);
  if(ScriptMemory_2.is_open()){
        ScriptMemory_2.close();
      }

  //Default values
  this->geo_path->value("");
  this->msh_path->value("");
}

//Closes geo Warning Window:
void WindowConstructor::Close1_f(){
  this->Container_4_1->hide();
}

//Closes msh Warning Window:
void WindowConstructor::Close2_f(){
  this->Container_4_2->hide();
}

//Closes msh Warning Window:
void WindowConstructor::Close3_f(){
  this->Container_4_3->hide();
}

//Allows the user to choose the weak form only if has choosen a Physical Group:
void WindowConstructor::free_WF_f(){
  this->WeakForm_Menu->activate();
}

//_________________________________________________________________________________
//Container_1

//Add .geo file path and set it to "gp":
void WindowConstructor::Input_geo_f() {

Fl_File_Chooser gpf(".", "*.geo", Fl_File_Chooser::SINGLE, "Select a .geo file");

//Stores the previous value.
string previous_gp = geo_path->value();
//Opens the file chooser
gpf.show();
    while (gpf.shown())
        Fl::wait();
    if (gpf.value() != nullptr) {
        this->geo_path->value(gpf.value());
        string gp_value = gpf.value();
    }
    else{
      /*Since the previous_gp is a string type class and the Fl_File_Output expects a C-style (a null-terminated array of chacarters,
      therefore, the last chacarter is "\O"), c_str() is necessary here to convert the std::string to a C-style string. */
      this->geo_path->value(previous_gp.c_str()); 
    }
}

//Add .msh file path and set it to "mp":
void WindowConstructor::Input_msh_f() {
  
Fl_File_Chooser mpf(".", "*.msh", Fl_File_Chooser::SINGLE, "Select a .msh file");

//Stores the previous value.
  string previous_mp = msh_path->value();

 //Opens the file chooser
mpf.show();
    while (mpf.shown())
        Fl::wait();
    if (mpf.value() != nullptr) {
        this->msh_path->value(mpf.value());
    }
    else{
      /*Since the previous_gp is a string type class and the Fl_File_Output expects a C-style (a null-terminated array of chacarters,
      therefore, the last chacarter is "\O"), c_str() is necessary here to convert the std::string to a C-style string. */
      this->msh_path->value(previous_mp.c_str()); 
    }
}


//Button to save path values and allows the user to use the Container_2
void WindowConstructor::SavePath_f(){

  gp = this->geo_path->value();
  mp = this->msh_path->value();
  playback->SetGeoPath(gp);
  playback->SetMshPath(mp);
  //playback->redraw();
  

  //Finds "gp" path file:
  std::ifstream file(gp);

  //If the file path is valid:
  if (file.is_open()) {
    while (getline(file, line)) {

      // Verifys if the line beggins with "//"
      if (line.find("//") == 0) {

        //Ignores the line and follows to the next interation loop              
        continue;
      }

      //Declares the "match" name for each pattern found in the "gp" file:
      smatch match;

      //PhysicalGroups Identifier:
      //match[1] = Physical Group (Ex: Physical Curve)
      //match[3] = The Physical Group names (Ex: Left, Right...)
      //match[4] = matid 
      //match[2] = (Curve, Surface or Volume)

      if (std::regex_search(line, match, pattern)) {

        //Defines the "count" variable being equal to the pattern found by match[2]:
        std::string count = match[2];

        if (count == "Curve") {
          curveCount++;
        }

        else if (count == "Surface") {
          surfaceCount++;
        }

        else if (count == "Volume") {
          volumeCount++;
        }

        //Defines the "type" variable being equal to the pattern found by match[1]:
        std::string type = match[1];

        //Defines the "name" variable being equal to the pattern found by match[3]:
        string name = match[3];

        //Saves the "name" variables found in the "gp" file inside the "names" vector
        names.push_back(name);

        //Defines the "matid" variable being equal to the pattern found by match[4]:
        int matid = std::stoi(match[4]);

        //Saves the "matid" variables found in the "gp" file inside the "matids" vector
        matids.push_back(matid);

        //Adds the "name" variables saved in the vector to the Physical Groups Menu
        for (const auto& label : names) {
          this->PhysicalGroups->add(label.c_str(),0,(Fl_Callback*)WindowConstructor::freewf,0,0);

        }
      }
    }

    file.close();

    //___________________________________________________________________________________________________________________

    if (curveCount > 0 && surfaceCount == 0 && volumeCount == 0) {
      d =1;//Dimension = 1
      this->WF_Icon->activate();
    }

    else if (curveCount > 0 && surfaceCount > 0 && volumeCount == 0) {
      d = 2;//Dimension = 2
      this->WF_Icon->activate();
    }

    else if (curveCount > 0 && surfaceCount > 0 && volumeCount > 0) {
      d = 3;//Dimension = 3
      this->WF_Icon->activate();
    }

    else { 
      this->Container_4_1->show();
      return;
    }

    if(mp == "\0" & gp == "\0"){
      this->WF_Icon->deactivate();
      this->Container_4_3->show();
      buff2->text("");
    }

    else if(mp == "\0" & gp != "\0"){
      this->WF_Icon->deactivate();
      this->Container_4_2->show();
    }

    else if(gp == "\0" & mp != "\0"){
      this->WF_Icon->deactivate();
      this->Container_4_1->show();
    }

    //Dimension and output paths

    this->dimension->value(d); //Sets the output dimension value;

    //Defines a txt file where will be all script information
    fstream text_file; //fstream library: out to write mode (overwrites), app (append), in (read)

    string str; //string that saves a copie from text_file.

    scriptdisplay->buffer(buff2);
    std::ostringstream oss;

    //_________________________________________________________________________________________________________________

    //Read File and find gp (geo_path) and mp (msh_path) first character position

    ScriptMemory_2.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory_2.

    //Copies ScriptMemory_2 text to str.
    if(ScriptMemory_2.is_open()){ //Opens ScriptMemory_2.
      string line;
      while(getline(ScriptMemory_2, line)){
        str += line;
        str += "\n";
      }
      ScriptMemory_2.close();//Closes ScriptMemory_2.
    }

    size_t posg = str.find("geo path: ",0); //Finds and returns the first "Input_geo: " character position. If was not found, returns npos. pos = position.


    //Reads string str and finds the next NULL line position after "Input_geo: ".
    if(posg != string::npos){
      size_t pos_null_g = str.find("\n\n",posg); //find and return the NULL line position. If not found, returns npos.
      str.erase(posg, pos_null_g - posg+2); //Delete the pos_g and pos_null_g gap text.
    
      //Firt clear the ScriptMemory_2
      ScriptMemory_2.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
      if(ScriptMemory_2.is_open()){
        ScriptMemory_2.close();
      }

      ScriptMemory_2.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
      if(ScriptMemory_2.is_open()){

        ScriptMemory_2 << str; //Append str to ScriptMemory_2
        ScriptMemory_2 << "geo path: " << gp << "\n" << "msh path: " << mp << "\n" << "\n";
        ScriptMemory_2.close(); //Closes ScriptMemory_2.
        }
        buff2->remove(0,buff2->length());
        oss << str;
        oss << "geo path: " << gp << "\n" << "msh path: " << mp << "\n" << "\n";
        buff2->append(oss.str().c_str());

    }

    else{
      ScriptMemory_2.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
      ScriptMemory_2 << "geo path: " << gp << "\n" << "msh path: " << mp << "\n" << "\n";
      ScriptMemory_2.close(); //Closes ScriptMemory_2.
      
      buff2->remove(0,buff2->length());
      oss << str;
      oss << "geo path: " << gp << "\n" << "msh path: " << mp << "\n" << "\n";
      buff2->append(oss.str().c_str());
    }
  }

  playback->Setmatid(matids);
  //playback->DrawGeometry();
  
  playback->draw();
  playback->take_focus();
}

//____________________________________________________________________________________________________________
//Menu

