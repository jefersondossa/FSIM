#include "WindowConstructor.h"

void WindowConstructor::fRunFoam(){

    string MSH_PATH = msh->value();

    //Copys the path of .msh file to inside the "GeneratedFiles" directory and defines the string code to put on Ubunto(system comand) and convert the .msh file to format 2
    string cpMSH = "cd GeneratedFiles && cp -r " + MSH_PATH + " ." + " && msh_file=$(find . -type f -name '*.msh') && gmsh \"$msh_file\" -2 -format msh2 && gmshToFoam \"$msh_file\"";
    system(cpMSH.c_str());

    //Mesh verification
    //string mesh = "cd GeneratedFiles && checkMesh";
    //system(mesh.c_str());

    fstream RunFoam_file;
    string Line;

    //Open and edit the boundary file
    RunFoam_file.open("../build/GeneratedFiles/constant/polyMesh/boundary");
    string boundary_str;

    //Copies bondary text to Str.
    if(RunFoam_file.is_open()){ //Opens boundary.
      while(getline(RunFoam_file, Line)){
        boundary_str += Line;
        boundary_str += "\n";
      }
      RunFoam_file.close();//Closes RunFoam
    }

    vector <string> vCorrectboundary;
    vector <string> vinGroups;
    vector <string> vSlip; //Stores terms if "Slip" was found;

    //Include on vectors vCorrectboundary and vinGroups the correlated values since the vBoundaryCondition (Global Vector with zeroGradient, fixedValue, empty...)
        for(int i = 0; i < vBoundaryCondition.size(); i++){
            if(vBoundaryCondition[i] == "fixedValue" || vBoundaryCondition[i] == "zeroGradient" ){
                vCorrectboundary.push_back("type            patch");
                vinGroups.push_back("");
            }

            if(vBoundaryCondition[i] == "noSlip" || vBoundaryCondition[i] == "slip"){
                vCorrectboundary.push_back("type            wall");
                vinGroups.push_back("inGroups        List<word> 1(wall);");
            }

            if(vBoundaryCondition[i] == "slip"){
                vSlip.push_back("slip");
            }

            if(vBoundaryCondition[i] == "empty"){
  	             vCorrectboundary.push_back("type            empty");
                 vinGroups.push_back("inGroups        List<word> 1(empty);");       
            }
        }

        for(int i = 0; i < vCorrectboundary.size(); i++){

            size_t pos = boundary_str.find(vPhysicalGroup[i],0); //Function to find the position of PhysicalGroups include on vector vPhysicalGroup start on position 0

            size_t pos_type = boundary_str.find("type            patch",pos); //Function to find the position of "type            patch" after the position of PhysicalGroups

            size_t pos_inGroups = boundary_str.find("physicalType    patch;",pos_type); //Function to find the position of "physicalType    patch;" after the position of "type            patch" 

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

    if(vSlip.size() == 0){

        string foamRun = "cd GeneratedFiles && foamRun";
        system(foamRun.c_str());
    }

    if(vSlip.size() > 0){

        // string foamRun = "cd GeneratedFiles && foamRun";
        // system(foamRun.c_str());

         fstream Slip;
         string decompose;
         string Line;

         Slip.open("../GraphicInterface/TopHeaderMenu/GenerationFunctions/Parameters_Header/decomposeParDict.txt",ios::in);
         if(Slip.is_open()){ //Copies text from headerfromfv to a string fvheader.
             while(getline(Slip, Line)){
                 decompose += Line + "\n";
             }

             Slip.close();
         }

         system("cd GeneratedFiles && cd system && touch decomposeParDict");

         Slip.open(("../build/GeneratedFiles/system/decomposeParDict"),ios::app);
         if (Slip.is_open()){
             Slip << decompose;
             Slip.close();
     }
         string decomposePar = "cd GeneratedFiles && decomposePar";
         //Here, it's necessary edit the number of processors (NUMBERP):
         //1. string foamRun: "cd GeneratedFiles && mpirun -np NUMBERP foamRun -parallel"
         //2. decomposeParDict.txt: numberOfSubdomains NUMBERP and simpleCoeffs {n               (1 2 2);} | 1*2*2 = NUMBERP;
         string foamRunparallel = "cd GeneratedFiles && mpirun -np 4 foamRun -parallel";
         system(decomposePar.c_str());
         system(foamRunparallel.c_str());
    }

    string foam = "cd GeneratedFiles && touch project.foam";
    system(foam.c_str());

    if(vSlip.size() > 0){

        string reconstruct = "cd GeneratedFiles && reconstructPar";
        system(reconstruct.c_str());
    }

    //Post-Process

    string p_cell = "cd GeneratedFiles && foamPostProcess -solver incompressibleFluid -func wallShearStress";
    system(p_cell.c_str()); //Calculates the Shear Stress


    string vv_cell = "cd GeneratedFiles && foamPostProcess -solver incompressibleFluid -func vorticity";
    system(vv_cell.c_str());

    string u_cell = "cd GeneratedFiles && foamPostProcess -solver incompressibleFluid -func \"cellMaxMag(field = U)\"";
    system(u_cell.c_str());

    string um_cell = "cd GeneratedFiles && foamPostProcess -solver incompressibleFluid -func \"cellMinMag(field = U)\"";
    system(um_cell.c_str());

    

}
