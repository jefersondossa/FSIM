/* This cpp file describes the following Analysis Menu functions:

1.  fLinear_cb();
2.  fNonLinear_cb(); 
3.  fIncrem_Transient_cb();
4.  fFluidFlow_cb();
5.  fApply_Analysis_cb(); 

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

/* ========================= fFluidFlow_cb ========================= */

void WindowConstructor::fFluidFlow_cb(){

    string ff = FluidFlow->text();

    if(ff == "Turbulent"){

        Reynolds_Number->activate();
        Ref_Length->activate();
    }

    if(ff == "Laminar"){

        Reynolds_Number->deactivate();
        Ref_Length->deactivate();
    }
}

/* ========================= fApply_Analysis_cb ========================= */

void WindowConstructor::fApply_Analysis_cb(){

    string Str;
    string Line;
    ostringstream Oss;

    string at = AnalysisType->text();
    string st = SolverType->text();
    string sit = FluidSimulation->text();
    string ft = FluidFlow->text();
    double a; 
    double b; 
    double c; 
    double d; 
    double e; 
    float f;

    ScriptDisplay->buffer(Buffer);

    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in);
        
         if(ScriptMemory_txt.is_open()){ //Opens ScriptMemory.
             while(getline(ScriptMemory_txt, Line)){
                 Str += Line;
                 Str += "\n";
             }   
             ScriptMemory_txt.close();//Closes ScriptMemory.
         }

     if(NonLinear_Menu->visible() == 0 & Increm_Transient_Menu->visible() == 0 & OpenFoam_Analysis_Menu->visible() == 0){

         a = NULL;
         b = NULL;
         c = NULL;

         
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
        string ff = FluidFlow->text();

        //Turbulence Variables:
        double re = Reynolds_Number->value();
        double l = Ref_Length->value();
        double k; // Turbulent kinetic energy [ m2s−2]
        double I; // Turbulence intensity [%]
        double er; // Turbulent kinetic energy dissipation rate [ m2s−3]

            //Find velocity uref:

            ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in);
            regex U("Velocity BC:\\s*\\((\\d+),\\s*(\\d+),\\s*(\\d+)\\)");
            smatch Match;
            double UX = 0, UY = 0, UZ = 0; //(x,y,z) Velocity

            if (ScriptMemory_txt.is_open()){
                while(getline(ScriptMemory_txt, Line)){
                    if (regex_search(Line, Match, U)) { 
                        
                        UX += stod(Match[1]);
                        UY += stod(Match[2]);
                        UZ += stod(Match[3]);
                    }
                }
                ScriptMemory_txt.close();
            }

        string turbulence;

        if(ff == "Laminar"){

            turbulence = "";
        }

        if(ff == "Turbulent"){

            double Uref = pow((UX*UX + UY*UY + UZ*UZ),1/2);
            I = 0.16*pow(re, -0.125);
            k = 1.5*pow((I*Uref),2);
            er = pow(0.09,0.75)*pow(k,1.5)/l;

            string Re = to_string(re);
            string L = to_string(l);
            string K = to_string(k);
            string E = to_string(er);

            turbulence =    "\nReynolds Number: " + Re + "\nReference Length: " + L + " [m]" +
                            "\nTurbulent kinetic energy: " + K + " [m2s-2]" + "\nTurbulent kinetic energy dissipation rate: " + E + " [m2s-3]";
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
                                   << "writePrecision: " << d << "\n" << "timePrecision: " << e << turbulence
                                   << "\n" << "\n";
                  ScriptMemory_txt.close(); //Closes ScriptMemory.
              }

              Buffer.remove(0,Buffer.length());
              Oss<< Str;
              Oss<< "Analysis Type: " << at << "\n" << "Simulation Type: " << sit << "\n" 
                 << "Flow Type: " << ft << "\n" << "Viscosity: " << f << "\n"
                 << "endTime: " << a << "\n" << "deltaT: " << b << "\n" << "writeInterval: " << c << "\n"
                 << "writePrecision: " << d << "\n" << "timePrecision: " << e << turbulence
                 << "\n" << "\n";
              Buffer.append(Oss.str().c_str());
         }

         else{
             ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
             ScriptMemory_txt << "Analysis Type: " << at << "\n" << "Simulation Type: " << sit << "\n" 
                              << "Flow Type: " << ft << "\n" << "Viscosity: " << f << "\n"
                              << "endTime: " << a << "\n" << "deltaT: " << b << "\n" << "writeInterval: " << c << "\n"
                              << "writePrecision: " << d << "\n" << "timePrecision: " << e << turbulence
                              << "\n" << "\n";
             ScriptMemory_txt.close(); //Closes ScriptMemory.

             Buffer.remove(0,Buffer.length());
             Oss<< Str;
             Oss<< "Analysis Type: " << at << "\n" << "Simulation Type: " << sit << "\n" 
                << "Flow Type: " << ft << "\n" << "Viscosity: " << f << "\n"
                << "endTime: " << a << "\n" << "deltaT: " << b << "\n" << "writeInterval: " << c << "\n"
                << "writePrecision: " << d << "\n" << "timePrecision: " << e << turbulence
                << "\n" << "\n";
             Buffer.append(Oss.str().c_str());
         }
     }
}
