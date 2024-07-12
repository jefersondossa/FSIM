#include "WindowConstructor.h"

void WindowConstructor::fRunFoam(){

    string MSH_PATH = msh->value();
    string cpMSH = "cd GeneratedFiles && cp -r " + MSH_PATH + " ." + " && msh_file=$(find . -type f -name '*.msh') && gmsh \"$msh_file\" -2 -format msh2 && gmshToFoam \"$msh_file\"";
    system(cpMSH.c_str());


    // string 
    // system("msh_file=$(find . -type f -name '*.msh')");
    // system("gmsh \"$msh_file\" -2 -format msh2");
    // system("gmshToFoam \"$msh_file\"");
}