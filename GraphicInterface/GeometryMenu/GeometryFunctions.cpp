/* This cpp file describes the following Geometry Menu functions:

1.  fOpen_Main_Window
3.  fInput_geo_cb
4.  fInput_msh_cb
5.  fApply_geometry_cb

*/

#include "WindowConstructor.h"

using namespace std;

/* ========================= fOpen_Main_Window() ========================= */

void WindowConstructor::fOpen_Main_Window(){

    this->Window->show(); //Opens the FSArl window.

    //This is necessary to every time clean the ScriptMemory.txt file.
    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);
        if(ScriptMemory_txt.is_open()){
        ScriptMemory_txt.close();
    }
}

/* ========================= fInput_geo() ========================= */

void WindowConstructor::fInput_geo_cb(){

    Fl_File_Chooser GP(".", "*.geo", Fl_File_Chooser::SINGLE, "Select a .geo file");

    string previous_geo = geo->value(); //Stores the previous value.

    GP.show(); //Opens the file chooser
    while (GP.shown()){
        Fl::wait();
        if (GP.value() != nullptr) {
            this->geo->value(GP.value());
        }

        else{
            /*Since the previous_geo is a string type class and the Fl_File_Output expects a C-style (a null-terminated array of chacarters,
            therefore, the last chacarter is "\O"), c_str() is necessary here to convert the std::string to a C-style string. */
            this->geo->value(previous_geo.c_str()); 
        }
    }
}

/* ========================= fInput_msh() ========================= */

void WindowConstructor::fInput_msh_cb(){

    Fl_File_Chooser MP(".", "*.msh", Fl_File_Chooser::SINGLE, "Select a .msh file");

    string previous_msh = msh->value(); //Stores the previous value.

    MP.show(); //Opens the file chooser
    while (MP.shown()){
        Fl::wait();
        if (MP.value() != nullptr) {
            this->msh->value(MP.value());
        }

        else{

          this->msh->value(previous_msh.c_str()); 
        }
    }
}

/* ========================= fApply_geometry() ========================= */

void WindowConstructor::fApply_geometry_cb(){


//**********************************************************************************************
/*This function is separated in 3 parts:

  1. Reads and writes the paths in the script and ScriptMemory_txt;
  2. Verifies the dimension problem and sets the Dimension output problem. 
     It also creates new Physical Groups names itens in the menu and frees the WeakForm use.
  3. Sent to the OpenGl Window the Physical Groups id's (matids) and the .geo and .msh paths. */

//**********************************************************************************************
  //Stores the .geo and .msh file paths.
  gp = geo->value();
  mp = msh->value();

  // if((gp == "" || mp != "") || (!gp.empty() && gp.substr(gp.size() - 4) != ".geo") || 
  //   (!mp.empty() && mp.substr(mp.size() - 4) != ".msh")){

  //   return;
  // }

  // else{

/* ========================= PART 1. ========================= */

  //Sends the text to the Script Display GUI
  ScriptDisplay->buffer(Buffer);

  //Read File and find gp (geo_path) and mp (msh_path) first character position
  ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory_txt.

  string Str;
  string Line;
  ostringstream Oss;

  if(ScriptMemory_txt.is_open()){ //Opens ScriptMemory_txt.
    while(getline(ScriptMemory_txt, Line)){
      Str += Line;
      Str += "\n";
    }
    ScriptMemory_txt.close();//Closes ScriptMemory_txt.
  }
  size_t posg = Str.find("geo path: ",0); //Finds and returns the first "Input_geo: " character position. If was not found, returns npos. pos = position.
  //Reads string str and finds the next NULL line position after "Input_geo: ".
  if(posg != string::npos){
    size_t pos_null_g = Str.find("\n\n",posg); //find and return the NULL line position. If not found, returns npos.
    Str.erase(posg, pos_null_g - posg+2); //Delete the pos_g and pos_null_g gap text.
  
    //Firt clears the ScriptMemory_txt
    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
    if(ScriptMemory_txt.is_open()){
      ScriptMemory_txt.close();
    }
    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
    if(ScriptMemory_txt.is_open()){
      ScriptMemory_txt << Str; //Append Str to ScriptMemory_txt
      ScriptMemory_txt << "geo path: " << gp << "\n" << "msh path: " << mp << "\n" << "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory_txt.
      }
      Buffer.remove(0,Buffer.length());
      Oss << Str;
      Oss << "geo path: " << gp << "\n" << "msh path: " << mp << "\n" << "\n";
      Buffer.append(Oss.str().c_str());

  }

  else{
    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
    ScriptMemory_txt << "geo path: " << gp << "\n" << "msh path: " << mp << "\n" << "\n";
    ScriptMemory_txt.close(); //Closes ScriptMemory_txt.
    
    Buffer.remove(0,Buffer.length());
    Oss << Str;
    Oss << "geo path: " << gp << "\n" << "msh path: " << mp << "\n" << "\n";
    Buffer.append(Oss.str().c_str());
  }

/* ========================= PART 2. ========================= */

  //Parameters for repetition count (curve, surface and volume):
  int curveCount = 0;
  int surfaceCount = 0;
  int volumeCount = 0;

  //First, finds geo path file:
  ifstream file(gp);

  //Clears the vector names to add new Physical Groups names when changing the files.
  names.clear();

  //Since the goal is to find the Physical Groups and store then:

    //Create the Physical Group pattern at geo file.
    regex PhysGroup_pattern(R"(Physical ((Curve|Surface|Volume))\(\"([^\"]+)\", (\d+)\))");

    //Physical Groups Identifier by regex library:
      // 1. match[1] = Physical Group (Ex: Physical Curve);
      // 2. match[2] = (Curve, Surface or Volume);
      // 3. match[3] = The Physical Group names (Ex: Left, Right...);
      // 4. match[4] = matid.

  //If the file path is valid:
  if (file.is_open()) {
    while (getline(file, Line)) {

      // Verifys if the line beggins with "//".
      if (Line.find("//") == 0) {

        //Ignores the line and follows to the next interation loop.              
        continue;
      }

      smatch Match;
      if (regex_search(Line, Match, PhysGroup_pattern)) {

        //Physical Groups by Objects:
          // 1. Defines the Physical Group "type" variable being equal to the pattern found by match[1]:
          string type = Match[1];

          // 2. Defines the "count" variable being equal to the pattern found by match[2]:
          string count = Match[2];

          // 3. Defines the "name" variable being equal to the pattern found by match[3]:
          string name = Match[3];

          // 4. Defines the "matid" variable being equal to the pattern found by match[4]:
          int matid = std::stoi(Match[4]);
          //curveCount = 6;
           if (count == "Curve") {
             curveCount++;
           }

          if (count == "Surface") {
            surfaceCount++;
          }

          else if (count == "Volume") {
            volumeCount++;
          }

          //Saves the "name" variables found in the "gp" file inside the "names" vector
          names.push_back(name);

          //Saves the "matid" variables found in the "gp" file inside the "matids" vector
          matids.push_back(matid);

          PhysGroup[name] = matid;
          PhysGroup2[name] = type;

          //Clears the Physical Group menu to add new itens when changing the files.
          for(int i = 1; i <= PhysicalGroups->size(); i++){
            PhysicalGroups->remove(i);
          }
          //Adds the "name" variables saved in the vector to the Physical Groups Menu
          for (const auto& label : names) {
            this->PhysicalGroups->add(label.c_str(),0,(Fl_Callback*)WindowConstructor::fStatic_WF_free,0,0);
        }
      }
    }

    file.close();
  }

  //___________________________________________________________________________________________________________________

    //Parameter for the Dimension:
    int d = 0;

    if (curveCount > 0 && surfaceCount == 0 && volumeCount == 0) {
      d = 1;//Dimension = 1
      this->Button_WeakForm->activate();
    }

    else if (surfaceCount > 0 && volumeCount == 0) {
      d = 2;//Dimension = 2
      this->Button_WeakForm->activate();
    }

    else if (volumeCount > 0) {
      d = 3;//Dimension = 3
      this->Button_WeakForm->activate();
    }

    //Dimension and output paths

    this->dimension->value(d); //Sets the output dimension value;
    Dimension = d;

/* ========================= PART 3. ========================= */

  playback->SetGeoPath(gp);
  playback->SetMshPath(mp);

  playback->Setmatid(matids);
  
  playback->draw();
  playback->redraw();
  playback->take_focus();
}


