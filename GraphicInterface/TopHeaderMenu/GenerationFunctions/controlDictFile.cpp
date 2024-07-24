#include "WindowConstructor.h"

void WindowConstructor::fcontrolDictFile(){

fstream cD_file;
string Line;

//Lets open the controlDictfile and write his header

    string cDheader;
    fstream headerFromcD;

    headerFromcD.open("../GraphicInterface/TopHeaderMenu/GenerationFunctions/Parameters_Header/headerfromcontrolDict.txt",ios::in);
    if(headerFromcD.is_open()){ //Copies text from headerFromcD to a string cDheader.
      while(getline(headerFromcD, Line)){
        cDheader += Line + "\n";
      }
      headerFromcD.close();
    }

    cD_file.open((systemFiles[0]),ios::app);
    if (cD_file.is_open()){

        cD_file << cDheader;
        cD_file.close();
    }

//Colocar os parametros editáveis do controlDict (endTime, deltaT, writeInterval, writePrecision, timePrecision)

    //Creation of matches that will store the goals values:
    regex rendTime(R"(endTime:\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?))");
    regex rdeltaT(R"(deltaT:\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?))");
    regex rwriteInterval(R"(writeInterval:\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?))");
    regex rwritePrecision(R"(writePrecision:\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?))");
    regex rtimePrecision(R"(timePrecision:\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?))");
    smatch Match;

    //Now, lets read the Script Memory and write the Physical Properties file:

    vector <string> vendTime;
    vector <string> vdeltaT;
    vector <string> vwriteInterval;
    vector <string> vwritePrecision;
    vector <string> vtimePrecision;

    fstream ScriptMemory_txt;
    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in);
    cD_file.open((systemFiles[0]),ios::app);

    if (ScriptMemory_txt.is_open()) {
        while (getline(ScriptMemory_txt, Line)) {
            if (regex_search(Line, Match, rendTime)) { 
                vendTime.push_back(Match[1]);
            }
            if (regex_search(Line, Match, rdeltaT)) { 
                vdeltaT.push_back(Match[1]);
            }
            if (regex_search(Line, Match, rwriteInterval)) { 
                vwriteInterval.push_back(Match[1]);
            }
            if (regex_search(Line, Match, rwritePrecision)) { 
                vwritePrecision.push_back(Match[1]);
            }
            if (regex_search(Line, Match, rtimePrecision)) { 
                vtimePrecision.push_back(Match[1]);
            }
                      
        }
    }

  if (cD_file.is_open()) {
      for(int i = 0; i < vendTime.size(); i++){
          cD_file  << "endTime         "<< vendTime[i] << ";" << "\n\n"
          << "deltaT          " << vdeltaT[i] << ";" << "\n\n"
          << "writeControl    timeStep;" << "\n\n"
          << "writeInterval   " << vwriteInterval[i] << "\n\n"
          << "purgeWrite      0;" << "\n\n"
          << "writeFormat     ascii;" << "\n\n"
          << "writePrecision  "  << vwritePrecision[i] << "\n\n"
          << "writeCompression off;" << "\n\n"
          << "timeFormat      general;" << "\n\n"
          << "timePrecision   " << vtimePrecision[i] << "\n\n"
          << "runTimeModifiable true;"
        //   << "\nfunctions {" 
        //   << "\nsurfaceFieldValue1 {"
        //   << "\ntype            surfaces;"
        //   << "\nlibs            (\"libfieldFunctionObjects.so\");"
        //   << "\nwriteControl    timeStep;"
        //   << "\nlog             true;"
        //   << "\nfields          (U);"
        //   << "\nwriteFields     yes;"
        //   << "\nregionType      patch;"
        //   << "\nname            edificio;"
        //   << "\noperation       areaIntegrate;"
        //   << "\nuseUserTime     true;"
        //   << "\ntimeStart       0;"
        //   << "\ntimeEnd         100;"
        //   << "}}"
          << "\n\n\n"
          << "// ************************************************************************* //" << "\n";
      }
  }
}