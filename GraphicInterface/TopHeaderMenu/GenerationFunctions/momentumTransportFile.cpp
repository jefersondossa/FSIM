#include "WindowConstructor.h"

void WindowConstructor::fmomentumTransportFile(){

  fstream ScriptMemory_txt;
  string Line;
  ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in);

  regex rff("Flow Type: (\\w+)");
  string ff;
  smatch Match;

  if (ScriptMemory_txt.is_open()) {
    while (getline(ScriptMemory_txt, Line)) {
        if (regex_search(Line, Match, rff)) { 
            ff = Match[1];
        }
    }
  }

    fstream mT_file;
    string mTheader;
    fstream headerFrommT;

    if(ff == "Laminar"){
      headerFrommT.open("../GraphicInterface/TopHeaderMenu/GenerationFunctions/Parameters_Header/LaminarMomentumTransport.txt",ios::in);
      if(headerFrommT.is_open()){ //Copies text from headerFrommT to a string mTheader.
        while(getline(headerFrommT, Line)){
          mTheader += Line + "\n";
        }
        headerFrommT.close();
      }

      mT_file.open((constantFiles[0]),ios::app);
      if (mT_file.is_open()){

          mT_file << mTheader;
          mT_file.close();
      }
    }

    if(ff == "Turbulent"){
      headerFrommT.open("../GraphicInterface/TopHeaderMenu/GenerationFunctions/Parameters_Header/TurbulentMomentumTransport.txt",ios::in);
      if(headerFrommT.is_open()){ //Copies text from headerFrommT to a string mTheader.
        while(getline(headerFrommT, Line)){
          mTheader += Line + "\n";
        }
        headerFrommT.close();
      }
  
      mT_file.open((constantFiles[0]),ios::app);
      if (mT_file.is_open()){
      
          mT_file << mTheader;
          mT_file.close();
      }
    }
}