#include "WindowConstructor.h"

void WindowConstructor::fUFile(){

    fstream U_file;

    //Lets open the Ufile and write his header

    string Uheader;
    fstream headerFromU;

    headerFromU.open("../GraphicInterface/TopHeaderMenu/GenerationFunctions/headerfromU.txt",ios::in);
    if(headerFromU.is_open()){ //Copies text from headerfromU to a string Uheader.
      while(getline(headerFromU, Line)){
        Uheader += Line + "\n";
      }
      headerFromU.close();
    }

    U_file.open((zeroFiles[0]),ios::app);
    if (U_file.is_open()){

        U_file << Uheader;
        U_file.close();
    }

    //Creation of patterns that will store the goals values:
    regex rPhysicalGroup("Physical Group: (\\w+)");
    regex rBoundaryCondition("Velocity BC: (\\w+)");
    regex rVector(R"(Uvector:\s*\(([^)]+)\))");

    //Now, lets read the Script Memory and write the U file:

    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in);
    U_file.open((zeroFiles[0]),ios::app);

    vector <string> vPhysicalGroup;
    vector <string> vBoundaryCondition;
    vector <string> vVector;
    vector <string> vCorrectVector; //Corrects the vVector, turning in to the real Uvector values.


    if (ScriptMemory_txt.is_open() && U_file.is_open()) {

        while (getline(ScriptMemory_txt, Line)) {
            if (regex_search(Line, Match, rPhysicalGroup)) { 

                vPhysicalGroup.push_back(Match[1]);
            }

            if (regex_search(Line, Match, rBoundaryCondition)) { 

                vBoundaryCondition.push_back(Match[1]);
            }
            if (regex_search(Line, Match, rVector)){

                vVector.push_back(Match[1]);
            }
        }
    }

    for(int i = 0; i < vBoundaryCondition.size(); i++){

            if(vBoundaryCondition[i] != "fixedValue"){

                vCorrectVector.push_back(" ");
            }

            else{

                int j = 0;
                vCorrectVector.push_back(vVector[j]);
                j++;
            }
        
    }

    if (U_file.is_open()) {
        for(int i = 0; i < vCorrectVector.size(); i++){
            U_file  << "    "<< vPhysicalGroup[i] << "\n" << "    {" << "\n"
            << "        type            " << vBoundaryCondition[i] << ";";

            if(vBoundaryCondition[i] == "fixedValue"){

                U_file << "\n" << "        value           uniform "
                << "(" << vCorrectVector[i] << ")" << ";" << "\n" << "    }" << "\n\n";
            }

            else{

                U_file << "\n" << "    }" << "\n\n";
            }

        }

        U_file << "}" << "\n\n" << "// ************************************************************************* //";
        U_file.close();
    }

}
