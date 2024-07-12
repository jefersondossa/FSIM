#include "WindowConstructor.h"

void WindowConstructor::fRunFoam(){

    string MSH_PATH = msh->value();
    string cpMSH = "cd GeneratedFiles && cp -r " + MSH_PATH + " ." + " && msh_file=$(find . -type f -name '*.msh') && gmsh \"$msh_file\" -2 -format msh2 && gmshToFoam \"$msh_file\"";
    system(cpMSH.c_str());

    fstream RunFoam_file;

    RunFoam_file.open("../build/GeneratedFiles/constant/polyMesh/boundary");
    string boundary_str;

    //Copies ScriptMemory text to Str.
    if(RunFoam_file.is_open()){ //Opens boundary.
      while(getline(RunFoam_file, Line)){
        boundary_str += Line;
        boundary_str += "\n";
      }
      RunFoam_file.close();//Closes RunFoam
    }

    vector <string> vCorrectboundary;


        for(int i = 0; i < vBoundaryCondition.size(); i++){
            if(vBoundaryCondition[i] == "fixedValue" | vBoundaryCondition[i] =="zeroGradient" ){
                vCorrectboundary.push_back("patch;");
            }

            if(vBoundaryCondition[i] != "noSlip"){
                vCorrectboundary.push_back("wall;");
            }
            else{
  	             vCorrectboundary.push_back("empty;");              
            }
        }

        for(int i = 0; i <vCorrectboundary.size(); i++){
            size_t pos = boundary_str.find(vPhysicalGroup[i],0);
            size_t pos_editavel = boundary_str.find("patch",pos);
            if(pos_editavel!= string::npos){
                boundary_str.replace(pos_editavel, 6, vCorrectboundary[i]);
            }

        }

    RunFoam_file.open("../build/GeneratedFiles/constant/polyMesh/boundary");
    if(RunFoam_file.is_open(),ios::out){
    RunFoam_file << boundary_str;
    RunFoam_file.close();
    }

    string foamRun = "cd GeneratedFiles && foamRun";
    system(foamRun.c_str());
}