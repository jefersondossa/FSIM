#include "WindowConstructor.h"

void WindowConstructor::fphysicalPropertiesFile(){

    fstream pP_file;
    string Line;

    //Lets open the pPfile and write his header

    string pPheader;
    fstream headerfromPp;

    headerfromPp.open("../GraphicInterface/TopHeaderMenu/GenerationFunctions/Parameters_Header/headerfromPhysicalProperties.txt",ios::in);
    if(headerfromPp.is_open()){ //Copies text from headerfromPp to a string pPheader.
      while(getline(headerfromPp, Line)){
        pPheader += Line + "\n";
      }
      headerfromPp.close();
    }

    pP_file.open((constantFiles[1]),ios::app);
    if (pP_file.is_open()){

        pP_file << pPheader;
        pP_file.close();
    }

    //Creation of patterns that will store the goals values:
    regex rNu(R"(Viscosity:\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?))");
    smatch Match;

    //Now, lets read the Script Memory and write the Physical Properties file:
    fstream ScriptMemory_txt; //GUI -> Just a fstream that is used all over the functions;
    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in);
    pP_file.open((constantFiles[1]),ios::app);

    vector <string> vNu;

    if (ScriptMemory_txt.is_open()) {

        while (getline(ScriptMemory_txt, Line)) {
            if (regex_search(Line, Match, rNu)) { 

                vNu.push_back(Match[1]);
            }
        }
    }

    if (pP_file.is_open()) {
        for(int i = 0; i < vNu.size(); i++){

            pP_file  << "nu              " << vNu[i] << ";";

        }

    pP_file << "\n\n" << "// ************************************************************************* //";
    pP_file.close();

    }
}

