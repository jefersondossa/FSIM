#include "WindowConstructor.h"



void WindowConstructor::fOFfolderGenerator(){

    fstream files;
    //**********************************************************************************************
    /*This function is separated in 4 parts (the default folders creations):

    1. Main OpenFoam Project Folder;
    2. 0 Folder; 
        2.1. U Folder;
        2.2. p Folder;
    3. system Folder;
        3.1. ControlDict Folder;
    4. constant Folder;
        4.1. momentumTransport;
        4.2. physicalProperties;
        4.3. polymesh

    //**********************************************************************************************/

/* ========================= PART 1. ========================= */


    system("mkdir GeneratedFiles");
    //Stores the txt file name
    string name = Script_Text_Name->value();

    //Sets the new folderName as OF_name
    string NewFolderName = "OF_" + name;

    //Path where the new folder will be created
    string defaultFolder = "GeneratedFiles";

    //Full path from new folder
    fs::path basePath = defaultFolder;
    fs::path OFPath = basePath / NewFolderName; //full path = OFPath

    
    // Checks if the new folder already exists
    if (fs::exists(OFPath)) {
        // Removes the existing folder and its contents
        fs::remove_all(OFPath);

        // Create the new directory
        fs::create_directory(OFPath);
    }

    else{
    // Creates the new directory
    fs::create_directory(OFPath);
    }
    
/* ========================= PART 2. ========================= */

    string zeroFolder = "0";
    fs::path zeroPath = OFPath / zeroFolder;

    {
    // Checks if the new folder already exists
    if (fs::exists(zeroPath)) {
        // Removes the existing folder and its contents
        fs::remove_all(zeroPath);

        // Create the new directory
        fs::create_directory(zeroPath);
    }

    else{
    // Creates the new directory
    fs::create_directory(zeroPath);
    }
    }
/* ========================= PART 2.1 & 2.2 ========================= */

    string zero = zeroPath.string(); //Stores the ZeroPath Folder in to a string
    vector<string> zeroNames; zeroNames = {"U", "p"};
    string zeroNamespath;

    for(int i = 0; i < zeroNames.size(); i++){

        zeroNamespath = zero + "/" + zeroNames[i];
        zeroFiles.push_back(zeroNamespath);
    }

    for(int i = 0; i <= zeroFiles.size(); i++){

        files.open(zeroFiles[i],ios::out);
        files.close();
    }

/* ========================= PART 3. ========================= */

    string systemFolder = "system";
    fs::path systemPath = OFPath / systemFolder;

    {
    // Checks if the new folder already exists
    if (fs::exists(systemPath)) {
        // Removes the existing folder and its contents
        fs::remove_all(systemPath);

        // Create the new directory
        fs::create_directory(systemPath);
    }

    else{
    // Creates the new directory
    fs::create_directory(systemPath);
    }
    }
///////////////////////////////////////////////////////////////

    string system = systemPath.string(); //Stores the ZeroPath Folder in to a string
    vector<string> systemNames; systemNames = {"controlDict"};
    string systemNamespath;

    for(int i = 0; i < systemNames.size(); i++){

        systemNamespath = system + "/" + systemNames[i];
        systemFiles.push_back(systemNamespath);
    }

    for(int i = 0; i <= systemFiles.size(); i++){

        files.open(systemFiles[i],ios::out);
        files.close();
    }

/* ========================= PART 4. ========================= */

    string constantFolder = "constant";
    fs::path constantPath = OFPath / constantFolder;

    {
    // Checks if the new folder already exists
    if (fs::exists(constantPath)) {
        // Removes the existing folder and its contents
        fs::remove_all(constantPath);

        // Create the new directory
        fs::create_directory(constantPath);
    }

    else{
    // Creates the new directory
    fs::create_directory(constantPath);
    }
    }
///////////////////////////////////////////////////////////////

    string constant = constantPath.string(); //Stores the ZeroPath Folder in to a string
    vector<string> constantNames; constantNames = {"momentumTransport", "physicalProperties"}; 
    string constantNamespath;
   
    for(int i = 0; i < constantNames.size(); i++){

        constantNamespath = constant + "/" + constantNames[i];
        constantFiles.push_back(constantNamespath);
    }

    for(int i = 0; i <= constantFiles.size(); i++){

        files.open(constantFiles[i],ios::out);
        files.close();
    }

}
