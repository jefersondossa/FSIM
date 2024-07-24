#include "WindowConstructor.h"

void WindowConstructor::ffv(){

    fstream fv_file;

    //Lets open the Ufile and write his header
    //smoothSolver

    string fvheader;
    fstream headerfromfv;
    string Line;

    headerfromfv.open("../GraphicInterface/TopHeaderMenu/GenerationFunctions/Parameters_Header/headerfromfvSchemes.txt",ios::in);
    if(headerfromfv.is_open()){ //Copies text from headerfromfv to a string fvheader.
      while(getline(headerfromfv, Line)){
        fvheader += Line + "\n";
      }
      headerfromfv.close();
    }

    fv_file.open((systemFiles[1]),ios::app);
    if (fv_file.is_open()){

        fv_file << fvheader;
        fv_file.close();
    }

    fvheader.clear();

    

    headerfromfv.open("../GraphicInterface/TopHeaderMenu/GenerationFunctions/Parameters_Header/headerfromfvSolution.txt",ios::in);
    if(headerfromfv.is_open()){ //Copies text from headerfromfv to a string fvheader.
      while(getline(headerfromfv, Line)){
        fvheader += Line + "\n";
      }
      headerfromfv.close();
    }

    fv_file.open((systemFiles[2]),ios::app);
    if (fv_file.is_open()){

        fv_file << fvheader;
        fv_file.close();
    }
}