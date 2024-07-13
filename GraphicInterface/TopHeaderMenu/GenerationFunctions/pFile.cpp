#include "WindowConstructor.h"

void WindowConstructor::fpFile(){

    fstream p_file;

    //Lets open the Ufile and write his header

    string pheader;
    fstream headerFromp;

    headerFromp.open("../GraphicInterface/TopHeaderMenu/GenerationFunctions/headerfromp.txt",ios::in);
    if(headerFromp.is_open()){ //Copies text from headerFromp to a string pheader.
      while(getline(headerFromp, Line)){
        pheader += Line + "\n";
      }
    }

    p_file.open((zeroFiles[1]),ios::app);
    if (p_file.is_open()){

        p_file << pheader;
        p_file.close();
    }

    //Creation of matches that will store the goals values:
    regex rPhysicalGroup("Physical Group: (\\w+)");
    regex rBoundaryCondition("Pressure BC: (\\w+)");
    regex rValue("pValue: ([+-]?([0-9]*[.])?[0-9]+)");

    vector <string> vPhysicalGroup;
    vector <string> vBoundaryCondition;
    vector <string> vValue;
    vector <string> vCorrectValue; //Corrects the vValue, turning in to the real p values.

    //Now, lets read the Script Memory and write the p file:
    fstream ScriptMemory_txt;
    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in);
    p_file.open((zeroFiles[1]),ios::app);

    if (ScriptMemory_txt.is_open() && p_file.is_open()) {
    while (getline(ScriptMemory_txt, Line)) {
        if (regex_search(Line, Match, rPhysicalGroup)) { 
            vPhysicalGroup.push_back(Match[1]);
        }
        if (regex_search(Line, Match, rBoundaryCondition)) { 
            vBoundaryCondition.push_back(Match[1]);
        }
        if (regex_search(Line, Match, rValue)){
            vValue.push_back(Match[1]);
        }
    }
}
for(int i = 0; i < vBoundaryCondition.size(); i++){
        if(vBoundaryCondition[i] != "fixedValue"){
            vCorrectValue.push_back(" ");
        }
        else{
            int j = 0;
            vCorrectValue.push_back(vValue[j]);
            j++;
        }
    
}
if (p_file.is_open()) {
    for(int i = 0; i < vCorrectValue.size(); i++){
        p_file  << "    "<< vPhysicalGroup[i] << "\n" << "    {" << "\n"
        << "        type            " << vBoundaryCondition[i] << ";";
        if(vBoundaryCondition[i] == "fixedValue"){
            p_file << "\n" << "        value           uniform "
            << "(" << vCorrectValue[i] << ")" << ";" << "\n" << "    }" << "\n\n";
        }
        else{
            p_file << "\n" << "    }" << "\n\n";
        }
    }
    p_file << "}" << "\n\n" << "// ************************************************************************* //";
    p_file.close();
}
}