#include "WindowConstructor.h"

void WindowConstructor::fRunFoam(){

    string MSH_PATH = msh->value();
    string cpMSH = "cd GeneratedFiles && cp -r " + MSH_PATH + " ." + " && msh_file=$(find . -type f -name '*.msh') && gmsh \"$msh_file\" -2 -format msh2 && gmshToFoam \"$msh_file\"";
    system(cpMSH.c_str());

    fstream RunFoam_file;
    string Line;

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
    vector <string> vinGroups;


        for(int i = 0; i < vBoundaryCondition.size(); i++){
            if(vBoundaryCondition[i] == "fixedValue" || vBoundaryCondition[i] == "zeroGradient" ){
                vCorrectboundary.push_back("type            patch");
                vinGroups.push_back("");
            }

            if(vBoundaryCondition[i] == "noSlip"){
                vCorrectboundary.push_back("type            wall");
                vinGroups.push_back("inGroups        List<word> 1(wall);");
            }
            if(vBoundaryCondition[i] == "empty"){
  	             vCorrectboundary.push_back("type            empty");
                 vinGroups.push_back("inGroups        List<word> 1(empty);");       
            }
        }

        for(int i = 0; i < vCorrectboundary.size(); i++){

            size_t pos = boundary_str.find(vPhysicalGroup[i],0);
            size_t pos_type = boundary_str.find("type            patch",pos);
            size_t pos_inGroups = boundary_str.find("physicalType    patch;",pos_type);

            if(pos_type != string::npos && pos_inGroups != string::npos){

                boundary_str.replace(pos_type, 21, vCorrectboundary[i]);

                if(vCorrectboundary[i] == "type            empty"){
                    boundary_str.replace(pos_inGroups, 22, vinGroups[i]);
                }

                if(vCorrectboundary[i] == "type            wall"){
                    boundary_str.replace(pos_inGroups-1, 22, vinGroups[i]);
                }

                if(vCorrectboundary[i] == "type            patch"){
                    boundary_str.replace(pos_inGroups, 22, vinGroups[i]);
                }
            }
            
        }

    RunFoam_file.open("../build/GeneratedFiles/constant/polyMesh/boundary");
    if(RunFoam_file.is_open(),ios::out){
    RunFoam_file << boundary_str;
    RunFoam_file.close();
    }

    string foamRun = "cd GeneratedFiles && foamRun";
    system(foamRun.c_str());

    string foam = "cd GeneratedFiles && touch project.foam";
    system(foam.c_str());
}