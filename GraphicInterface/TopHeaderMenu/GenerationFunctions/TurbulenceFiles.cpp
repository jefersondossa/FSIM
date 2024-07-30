#include "WindowConstructor.h"

 void WindowConstructor::fkFile(){

    fstream k_file;
    string Line;

    string kheader;
    fstream headerFrommk;
    fstream ScriptMemory_txt;

    headerFrommk.open("../GraphicInterface/TopHeaderMenu/GenerationFunctions/Parameters_Header/headerfromK.txt",ios::in);
    if(headerFrommk.is_open()){ //Copies text from headerFrommT to a string mTheader.
      while(getline(headerFrommk, Line)){
        kheader += Line + "\n";
      }
      headerFrommk.close();
    }

    k_file.open((zeroFiles[2]),ios::app);
    if (k_file.is_open()){

        k_file << kheader;
        k_file.close();
    }

    //Creation of patterns that will store the goals values:
    regex rK("Turbulent kinetic energy: ([+-]?([0-9]*[.])?[0-9]+)");
    smatch Match;
    string K;
    vector <string> vK_BC;

    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in);
    if (ScriptMemory_txt.is_open()) {

        while (getline(ScriptMemory_txt, Line)) {
            if (regex_search(Line, Match, rK)) { 

                K = Match[1];
            }
        }
    }

    for(int i = 0; i < vBoundaryCondition.size(); i++){

        vK_BC.push_back(vBoundaryCondition[i]);
    }

    for(int i = 0; i < vK_BC.size(); i++){

        if(vK_BC[i] != "fixedValue" && vK_BC[i] != "zeroGradient" && vK_BC[i] != "empty"){
            vK_BC[i] = "kqRWallFunction";
        }
    }
    
    k_file.open((zeroFiles[2]),ios::app);
    if (k_file.is_open()) {

        k_file << "internalField   uniform " << K << ";"
               << "\nboundaryField" << "\n{\n";

        for(int i = 0; i < vK_BC.size(); i++){
            k_file << "    "<< vPhysicalGroup[i] << "\n" << "    {" << "\n"
            << "        type            " << vK_BC[i] << ";";
            

            if(vK_BC[i] != "zeroGradient"){

                k_file << "\n" << "        value           uniform "
                << K << ";" << "\n" << "    }" << "\n\n";
            }

            else{

                k_file << "\n" << "    }" << "\n\n";
            }

        }

        k_file << "}" << "\n\n" << "// ************************************************************************* //";
        k_file.close();
    }

    
 }

 void WindowConstructor::fEpsilonFile(){

    fstream E_file;
    string Line;

    string Eheader;
    fstream headerFrommE;
    fstream ScriptMemory_txt;

    headerFrommE.open("../GraphicInterface/TopHeaderMenu/GenerationFunctions/Parameters_Header/headerfromEpsilon.txt",ios::in);
    if(headerFrommE.is_open()){ //Copies text from headerFrommT to a string mTheader.
      while(getline(headerFrommE, Line)){
        Eheader += Line + "\n";
      }
      headerFrommE.close();
    }

    E_file.open((zeroFiles[3]),ios::app);
    if (E_file.is_open()){

        E_file << Eheader;
        E_file.close();
    }

    //Creation of patterns that will store the goals values:
    regex rE("Turbulent kinetic energy dissipation rate: ([+-]?([0-9]*[.])?[0-9]+)");
    smatch Match;
    string E;
    vector <string> vE_BC;

    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in);
    if (ScriptMemory_txt.is_open()) {

        while (getline(ScriptMemory_txt, Line)) {
            if (regex_search(Line, Match, rE)) { 

                E = Match[1];
            }
        }
    }

    for(int i = 0; i < vBoundaryCondition.size(); i++){

        vE_BC.push_back(vBoundaryCondition[i]);
    }

    for(int i = 0; i < vE_BC.size(); i++){

        if(vE_BC[i] != "fixedValue" && vE_BC[i] != "zeroGradient" && vE_BC[i] != "empty"){
            vE_BC[i] = "epsilonWallFunction";
        }
    }
    
    E_file.open((zeroFiles[3]),ios::app);
    if (E_file.is_open()) {

        E_file << "internalField   uniform " << E << ";"
               << "\nboundaryField" << "\n{\n";

        for(int i = 0; i < vE_BC.size(); i++){
            E_file << "    "<< vPhysicalGroup[i] << "\n" << "    {" << "\n"
            << "        type            " << vE_BC[i] << ";";
            

            if(vE_BC[i] != "zeroGradient"){

                E_file << "\n" << "        value           uniform "
                << E << ";" << "\n" << "    }" << "\n\n";
            }

            else{

                E_file << "\n" << "    }" << "\n\n";
            }

        }

        E_file << "}" << "\n\n" << "// ************************************************************************* //";
        E_file.close();
    }
 }

 void WindowConstructor::fnutFile(){

    fstream nut_file;
    string Line;

    string nutheader;
    fstream headerFrommnut;
    fstream ScriptMemory_txt;

    headerFrommnut.open("../GraphicInterface/TopHeaderMenu/GenerationFunctions/Parameters_Header/headerfromnut.txt",ios::in);
    if(headerFrommnut.is_open()){ //Copies text from headerFrommT to a string mTheader.
      while(getline(headerFrommnut, Line)){
        nutheader += Line + "\n";
      }
      headerFrommnut.close();
    }

    nut_file.open((zeroFiles[4]),ios::app);
    if (nut_file.is_open()){

        nut_file << nutheader;
        nut_file.close();
    }

    //Creation of patterns that will store the goals values:
    string nut = to_string(0);
    vector <string> vnut_BC;

    for(int i = 0; i < vBoundaryCondition.size(); i++){

        vnut_BC.push_back(vBoundaryCondition[i]);
    }

    for(int i = 0; i < vnut_BC.size(); i++){

        if(vnut_BC[i] != "fixedValue" && vnut_BC[i] != "zeroGradient" && vnut_BC[i] != "empty"){
            vnut_BC[i] = "nutkWallFunction";
        }

        if(vnut_BC[i] == "fixedValue" || vnut_BC[i] == "zeroGradient"){
            vnut_BC[i] = "calculated";
        }
    }
    
    nut_file.open((zeroFiles[4]),ios::app);
    if (nut_file.is_open()) {

        nut_file << "internalField   uniform " << nut << ";"
               << "\nboundaryField" << "\n{\n";

        for(int i = 0; i < vnut_BC.size(); i++){
            nut_file << "    "<< vPhysicalGroup[i] << "\n" << "    {" << "\n"
            << "        type            " << vnut_BC[i] << ";";
            nut_file << "\n" << "        value           uniform "
            << nut << ";" << "\n" << "    }" << "\n\n";
            
        }

        nut_file << "}" << "\n\n" << "// ************************************************************************* //";
        nut_file.close();
    }
 }