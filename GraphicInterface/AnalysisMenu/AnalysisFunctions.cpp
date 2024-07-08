/* This cpp file describes the following Analysis Menu functions:

1.  fLinear_cb();
2.  fNonLinear_cb(); 
3.  fIncrem_Transient_cb();
4.  fApply_Analysis_cb(); 

*/

#include "WindowConstructor.h"

/* ========================= fLinear_cb ========================= */

void WindowConstructor::fLinear_cb(){

    NonLinear_Menu->hide();
    Increm_Transient_Menu->hide();
    OpenFoam_Analysis_Menu->hide();
    SolverType->activate();
}

/* ========================= fNonLinear_cb ========================= */

void WindowConstructor::fNonLinear_cb(){

    NonLinear_Menu->show();
    Increm_Transient_Menu->hide();
    OpenFoam_Analysis_Menu->hide();
    SolverType->activate();
}

/* ========================= fIncrem_Transient_cb ========================= */

void WindowConstructor::fIncrem_Transient_cb(){

    NonLinear_Menu->hide();
    Increm_Transient_Menu->show();
    OpenFoam_Analysis_Menu->hide();
    SolverType->activate();
}

/* ========================= fOFAnalysis_cb ========================= */

void WindowConstructor::fOFAnalysis_cb(){

    NonLinear_Menu->hide();
    Increm_Transient_Menu->hide();
    OpenFoam_Analysis_Menu->show();
    SolverType->deactivate();
}

/* ========================= fApply_Analysis_cb ========================= */

void WindowConstructor::fApply_Analysis_cb(){

    Str = "";
    Oss.str("");

     string at = this->AnalysisType->text();
     string st  = this->SolverType->text();
     string sit = FluidSimulation->text();
     string ft = FluidFlow->text();
     double a; 
     double b; 
     double c; 
     double d; 
     double e; 
     float f; 

     ScriptDisplay->buffer(Buffer);

     if(NonLinear_Menu->visible() == 0 & Increm_Transient_Menu->visible() == 0 & OpenFoam_Analysis_Menu->visible() == 0){

         a = NULL;
         b = NULL;
         c = NULL;

         ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in);

         if(ScriptMemory_txt.is_open()){ //Opens ScriptMemory.
             while(getline(ScriptMemory_txt, Line)){
                 Str += Line;
                 Str += "\n";
             }   
             ScriptMemory_txt.close();//Closes ScriptMemory.
         }

         size_t pos = Str.find("Analysis Type",0);

         if(pos != string::npos){
             size_t pos_null = Str.find("\n\n",pos);
             Str.erase(pos, pos_null - pos +2); 
      

             //Firt clear the ScriptMemory
             ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
             if(ScriptMemory_txt.is_open()){
                 ScriptMemory_txt.close();
             }

              ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
              if(ScriptMemory_txt.is_open()){
      
                  ScriptMemory_txt << Str; //Append Str to ScriptMemory
                  ScriptMemory_txt << "Analysis Type: " << at << "\n" << "Solver Type: " << st << "\n"<< "\n";
                  ScriptMemory_txt.close(); //Closes ScriptMemory.
              }

              Buffer.remove(0,Buffer.length());
              Oss<< Str;
              Oss<<"Analysis Type: " << at << "\n" << "Solver Type: "<< st << "\n"<< "\n";
              Buffer.append(Oss.str().c_str());
         }

         else{
             ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
             ScriptMemory_txt << "Analysis Type: " << at << "\n" << "Solver Type: "<< st << "\n"<< "\n";
             ScriptMemory_txt.close(); //Closes ScriptMemory.

             Buffer.remove(0,Buffer.length());
             Oss<< Str;
             Oss<<"Analysis Type: " << at << "\n" << "Solver Type: "<< st << "\n"<< "\n";
             Buffer.append(Oss.str().c_str());
         }
     }

     if(this->NonLinear_Menu->visible() != 0){

         a = this->tolerance->value();
         b = this->maxInterations->value();
         c = NULL;

         ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in);

         if(ScriptMemory_txt.is_open()){ //Opens ScriptMemory.
             while(getline(ScriptMemory_txt, Line)){
                 Str += Line;
                 Str += "\n";
             }   
             ScriptMemory_txt.close();//Closes ScriptMemory.
         }

         size_t pos = Str.find("Analysis Type",0);

         if(pos != string::npos){
             size_t pos_null = Str.find("\n\n",pos);
             Str.erase(pos, pos_null - pos +2); 
      

             //Firt clear the ScriptMemory
             ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
             if(ScriptMemory_txt.is_open()){
                 ScriptMemory_txt.close();
             }

              ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
              if(ScriptMemory_txt.is_open()){
      
                  ScriptMemory_txt << Str; //Append Str to ScriptMemory
                  ScriptMemory_txt << "Analysis Type: " << at << "\n" << "Solver Type: " << st << "\n"<< "Tolerance: " << a << "\n"<< "Maximum Interations: " << b << "\n"<< "\n";
                  ScriptMemory_txt.close(); //Closes ScriptMemory.
              }

              Buffer.remove(0,Buffer.length());
              Oss<< Str;
              Oss<<"Analysis Type: " << at << "\n" << "Solver Type: " << st << "\n" << "Tolerance: " << a << "\n"<< "Maximum Interations: " << b << "\n"<< "\n";
              Buffer.append(Oss.str().c_str());
         }

         else{
             ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
             ScriptMemory_txt << "Analysis Type: " << at << "\n" << "Solver Type: " << st << "\n" << "Tolerance: " << a << "\n"<< "Maximum Interations: " << b << "\n"<< "\n";
             ScriptMemory_txt.close(); //Closes ScriptMemory.

             Buffer.remove(0,Buffer.length());
             Oss<< Str;
             Oss<<"Analysis Type: " << at << "\n" << "Solver Type: " << st << "\n" << "Tolerance: " << a << "\n"<< "Maximum Interations: " << b << "\n"<< "\n";
             Buffer.append(Oss.str().c_str());
         }
     }

     if(this->Increm_Transient_Menu->visible() != 0){

         a = this->tolerance->value();
         b = this->maxInterations->value();
         c = this->Nsteps->value();

         ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in);

         if(ScriptMemory_txt.is_open()){ //Opens ScriptMemory.
         string line;
             while(getline(ScriptMemory_txt, line)){
                 Str += line;
                 Str += "\n";
             }   
             ScriptMemory_txt.close();//Closes ScriptMemory.
         }

         size_t pos = Str.find("Analysis Type",0);

         if(pos != string::npos){
             size_t pos_null = Str.find("\n\n",pos);
             Str.erase(pos, pos_null - pos +2); 
      

             //Firt clear the ScriptMemory
             ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
             if(ScriptMemory_txt.is_open()){
                 ScriptMemory_txt.close();
             }

              ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
              if(ScriptMemory_txt.is_open()){
      
                  ScriptMemory_txt << Str; //Append Str to ScriptMemory
                  ScriptMemory_txt << "Analysis Type: " << at << "\n" << "Solver Type: " << st << "\n"  << "Steps Number: " << c << "\n"<< "Tolerance: " << a << "\n"<< "Maximum Interations: " << b << "\n"<< "\n";
                  ScriptMemory_txt.close(); //Closes ScriptMemory.
              }

              Buffer.remove(0,Buffer.length());
              Oss<< Str;
              Oss<<"Analysis Type: " << at << "\n" << "Solver Type: " << st << "\n"  <<"Steps Number: " << c << "\n"<< "Tolerance: " << a << "\n"<< "Maximum Interations: " << b << "\n"<< "\n";
              Buffer.append(Oss.str().c_str());
         }

         else{
             ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
             ScriptMemory_txt << "Analysis Type: " << at << "\n" << "Solver Type: " << st << "\n"  << "Steps Number: " << c << "\n"<< "Tolerance: " << a << "\n"<< "Maximum Interations: " << b << "\n"<< "\n";
             ScriptMemory_txt.close(); //Closes ScriptMemory.

             Buffer.remove(0,Buffer.length());
             Oss<< Str;
             Oss<<"Analysis Type: " << at << "\n" << "Solver Type: " << st << "\n"  << "Steps Number: " << c << "\n"<< "Tolerance: " << a << "\n"<< "Maximum Interations: " << b << "\n"<< "\n";
             Buffer.append(Oss.str().c_str());
         }
     }

     if(this->OpenFoam_Analysis_Menu->visible() != 0){

        f = OFViscosity->value();
        a = endTime->value();
        b = deltaT->value();
        c = writeInterval->value();
        d = writePrecision->value();
        e = timePrecision->value();

         ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in);

         if(ScriptMemory_txt.is_open()){ //Opens ScriptMemory.
         string line;
             while(getline(ScriptMemory_txt, line)){
                 Str += line;
                 Str += "\n";
             }   
             ScriptMemory_txt.close();//Closes ScriptMemory.
         }

         size_t pos = Str.find("Analysis Type",0);

         if(pos != string::npos){
             size_t pos_null = Str.find("\n\n",pos);
             Str.erase(pos, pos_null - pos +2); 
      

             //Firt clear the ScriptMemory
             ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
             if(ScriptMemory_txt.is_open()){
                 ScriptMemory_txt.close();
             }

              ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
              if(ScriptMemory_txt.is_open()){
      
                  ScriptMemory_txt << Str; //Append Str to ScriptMemory
                  ScriptMemory_txt << "Analysis Type: " << at << "\n" << "Simulation Type: " << sit << "\n" 
                                   << "Flow Type: " << ft << "\n" << "Viscosity: " << f << "\n"
                                   << "endTime: " << a << "\n" << "deltaT: " << b << "\n" << "writeInterval: " << c << "\n" 
                                   << "writePrecision: " << d << "\n" << "timePrecision: " << e << "\n" << "\n";
                  ScriptMemory_txt.close(); //Closes ScriptMemory.
              }

              Buffer.remove(0,Buffer.length());
              Oss<< Str;
              Oss<< "Analysis Type: " << at << "\n" << "Simulation Type: " << sit << "\n" 
                 << "Flow Type: " << ft << "\n" << "Viscosity: " << f << "\n"
                 << "endTime: " << a << "\n" << "deltaT: " << b << "\n" << "writeInterval: " << c << "\n"
                 << "writePrecision: " << d << "\n" << "timePrecision: " << e << "\n" << "\n";
              Buffer.append(Oss.str().c_str());
         }

         else{
             ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
             ScriptMemory_txt << "Analysis Type: " << at << "\n" << "Simulation Type: " << sit << "\n" 
                              << "Flow Type: " << ft << "\n" << "Viscosity: " << f << "\n"
                              << "endTime: " << a << "\n" << "deltaT: " << b << "\n" << "writeInterval: " << c << "\n"
                              << "writePrecision: " << d << "\n" << "timePrecision: " << e << "\n" << "\n";
             ScriptMemory_txt.close(); //Closes ScriptMemory.

             Buffer.remove(0,Buffer.length());
             Oss<< Str;
             Oss<< "Analysis Type: " << at << "\n" << "Simulation Type: " << sit << "\n" 
                << "Flow Type: " << ft << "\n" << "Viscosity: " << f << "\n"
                << "endTime: " << a << "\n" << "deltaT: " << b << "\n" << "writeInterval: " << c << "\n"
                << "writePrecision: " << d << "\n" << "timePrecision: " << e << "\n" << "\n";
             Buffer.append(Oss.str().c_str());
         }
     }
}
