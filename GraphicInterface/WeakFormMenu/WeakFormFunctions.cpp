/* This cpp file describes the following Weak Form Menu functions:

01. fApply_WF_cb
02. fElasticity2D_cb(); 
03. fElasticityPositional2D_cb();
04. fElasticity3D_cb(); 
05. fElasticTruss_cb(); 
06. fPositionalTruss_cb(); 
07. fStokes_cb(); 
08. fNavierStokes_cb(); 
09. fPoisson_cb();
10. fL2Projection_cb(); 

*/

#include "WindowConstructor.h"
#include "WeakFormMenu/Elasticity2DMenu.cpp"

using namespace std;

/* ========================= fApply_WF_cb ========================= */

void WindowConstructor::fApply_WF_cb(){

    string name = WeakForm->text();

    //Clears str
    Str = "";
    Oss.str("");

    double a; 
    double b; 
    double c; 
    double d; 
  
    string wf; 
    string pg; 
    string bc; 
  
  ScriptDisplay->buffer(Buffer);

  //Elasticity2D Variables
  if(name == "Elasticity 2D"){

    //Variables Values
    a = this->young->value();
    b = this->poisson->value();
    c = NULL;

    d= this->plane_stress->value();
    //Read File and find PhysicalGroups->text() first character position

    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory.

    //Copies ScriptMemory text to Str.
    if(ScriptMemory_txt.is_open()){ //Opens ScriptMemory.
      while(getline(ScriptMemory_txt, Line)){
        Str += Line;
        Str += "\n";
      }
      ScriptMemory_txt.close();//Closes ScriptMemory.
    }

    size_t pos = Str.find(PhysicalGroups->text(),0); //Finds and returns the first Pyisical Groups character position. If was not found, returns npos. pos = position.

    //Reads string Str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = Str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      Str.erase(pos-17, pos_null - pos +19); //Delete the pos and pos_null gap text (-17 because it's necessary to erase the "Physical Group: ").

      /*
      At this moment, Str is modified. Since we want the text to be copied to the ScriptMemory (for future modifications in Str and header),
      besides the script print, the Str is copied to ScriptMemory and Oss(Oss-> buff -> scriptdisplay).
      */
      
      //Firt clear the ScriptMemory
      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
      if(ScriptMemory_txt.is_open()){
        ScriptMemory_txt.close();
      }

      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
      if(ScriptMemory_txt.is_open()){

      ScriptMemory_txt << Str; //Append Str to ScriptMemory
      ScriptMemory_txt << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n" << "Plane Stress: " << d << "\n"<< "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory.
      }
      Buffer.remove(0,Buffer.length());
      Oss<< Str;
      Oss<<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n" << "Plane Stress: " << d << "\n"<< "\n";
      Buffer.append(Oss.str().c_str());

    }
    else{
      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
      ScriptMemory_txt << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n" << "Plane Stress: " << d << "\n"<< "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory.
      
      Buffer.remove(0,Buffer.length());
      Oss<< Str;
      Oss<<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n" << "Plane Stress: " << d << "\n"<< "\n";
      Buffer.append(Oss.str().c_str());
    }
  }

  //ElasticityPositional2D Variables
  if(name == "Elasticity Positional 2D"){

    //Variables Values
    a = this->young->value();
    b = this->poisson->value();
    c = NULL;
    d = this->plane_stress->value();

    wf = this->WeakForm->label();
    pg = this->PhysicalGroups->text();
    bc = "";

    //Read File and find PhysicalGroups->text() first character position

    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory.

    //Copies ScriptMemory text to Str.
    if(ScriptMemory_txt.is_open()){ //Opens ScriptMemory.
      
      while(getline(ScriptMemory_txt, Line)){
        Str += Line;
        Str += "\n";
      }
      ScriptMemory_txt.close();//Closes ScriptMemory.
    }

    size_t pos = Str.find(PhysicalGroups->text(),0); //Finds and returns the first Pyisical Groups character position. If was not found, returns npos. pos = position.

    //Reads string Str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = Str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      Str.erase(pos-17, pos_null - pos +19); //Delete the pos and pos_null gap text (-17 because it's necessary to erase the "Physical Group: ").

      /*
      At this moment, Str is modified. Since we want the text to be copied to the ScriptMemory (for future modifications in Str and header),
      besides the script print, the Str is copied to ScriptMemory and Oss(Oss-> buff -> scriptdisplay).
      */
      
      //Firt clear the ScriptMemory
      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
      if(ScriptMemory_txt.is_open()){
        ScriptMemory_txt.close();
      }

      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
      if(ScriptMemory_txt.is_open()){

      ScriptMemory_txt << Str; //Append Str to ScriptMemory
      ScriptMemory_txt << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n" << "Plane Stress: " << d << "\n"<< "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory.
      }
      Buffer.remove(0,Buffer.length());
      Oss<< Str;
      Oss<<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n" << "Plane Stress: " << d << "\n"<< "\n";
      Buffer.append(Oss.str().c_str());

    }
    else{
      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
      ScriptMemory_txt << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n" << "Plane Stress: " << d << "\n"<< "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory.
      
      Buffer.remove(0,Buffer.length());
      Oss<< Str;
      Oss<<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n" << "Plane Stress: " << d << "\n"<< "\n";
      Buffer.append(Oss.str().c_str());
    }
  }

  //Elasticity3D Variables
  if(name == "Elasticity 3D"){

    a = this->young->value();
    b = this->poisson->value();
    c = NULL;

    wf = this->WeakForm->label();
    pg = this->PhysicalGroups->text();
    bc = "";

  //Read File and find PhysicalGroups->text() first character position

    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory.

    //Copies ScriptMemory text to Str.
    if(ScriptMemory_txt.is_open()){ //Opens ScriptMemory.
      
      while(getline(ScriptMemory_txt, Line)){
        Str += Line;
        Str += "\n";
      }
      ScriptMemory_txt.close();//Closes ScriptMemory.
    }


    size_t pos = Str.find(PhysicalGroups->text(),0); //Finds and returns the first Pyisical Groups character position. If was not found, returns npos. pos = position.

    //Reads string Str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = Str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      Str.erase(pos-17, pos_null - pos +19); //Delete the pos and pos_null gap text (-17 because it's necessary to erase the "Physical Group: ").

      /*
      At this moment, Str is modified. Since we want the text to be copied to the ScriptMemory (for future modifications in Str and header),
      besides the script print, the Str is copied to ScriptMemory and Oss(Oss-> buff -> scriptdisplay).
      */
      
      //Firt clear the ScriptMemory
      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
      if(ScriptMemory_txt.is_open()){
        ScriptMemory_txt.close();
      }

      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
      if(ScriptMemory_txt.is_open()){

      ScriptMemory_txt << Str; //Append Str to ScriptMemory
      ScriptMemory_txt << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n"<< "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory.
      }
      Buffer.remove(0,Buffer.length());
      Oss<< Str;
      Oss<<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n"<< "\n";
      Buffer.append(Oss.str().c_str());

    }
    else{
      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
      ScriptMemory_txt << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n"<< "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory.
      
      Buffer.remove(0,Buffer.length());
      Oss<< Str;
      Oss<<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n"<< "\n";
      Buffer.append(Oss.str().c_str());
    }
  }

  //ElasticTruss Variables
  if(name == "Elastic Truss"){

    a = this->young->value();
    b = this->area->value(); 
    c = NULL;

    wf = this->WeakForm->label();
    pg = this->PhysicalGroups->text();
    bc = "";

   //Read File and find PhysicalGroups->text() first character position

    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory.

    //Copies ScriptMemory text to Str.
    if(ScriptMemory_txt.is_open()){ //Opens ScriptMemory.
      
      while(getline(ScriptMemory_txt, Line)){
        Str += Line;
        Str += "\n";
      }
      ScriptMemory_txt.close();//Closes ScriptMemory.
    }

    size_t pos = Str.find(PhysicalGroups->text(),0); //Finds and returns the first Pyisical Groups character position. If was not found, returns npos. pos = position.

    //Reads string Str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = Str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      Str.erase(pos-17, pos_null - pos +19); //Delete the pos and pos_null gap text (-17 because it's necessary to erase the "Physical Group: ").

      /*
      At this moment, Str is modified. Since we want the text to be copied to the ScriptMemory (for future modifications in Str and header),
      besides the script print, the Str is copied to ScriptMemory and Oss(Oss-> buff -> scriptdisplay).
      */
      
      //Firt clear the ScriptMemory
      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
      if(ScriptMemory_txt.is_open()){
        ScriptMemory_txt.close();
      }

      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
      if(ScriptMemory_txt.is_open()){

      ScriptMemory_txt << Str; //Append Str to ScriptMemory
      ScriptMemory_txt << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Area: " << b << "\n"<< "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory.
      }
      Buffer.remove(0,Buffer.length());
      Oss<< Str;
      Oss<<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Area: " << b << "\n"<< "\n";
      Buffer.append(Oss.str().c_str());

    }
    else{
      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
      ScriptMemory_txt << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Area: " << b << "\n"<< "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory.
      
      Buffer.remove(0,Buffer.length());
      Oss<< Str;
      Oss<<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Area: " << b << "\n"<< "\n";
      Buffer.append(Oss.str().c_str());
    }
  }  
  

  //PositionalTruss Variables
  if(name == "Positional Truss"){

    a = this->young->value();
    b = this->area->value();
    c = NULL;

    wf = this->WeakForm->label();
    pg = this->PhysicalGroups->text();
    bc = "";

  //Read File and find PhysicalGroups->text() first character position

    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory.

    //Copies ScriptMemory text to Str.
    if(ScriptMemory_txt.is_open()){ //Opens ScriptMemory.
  
      while(getline(ScriptMemory_txt, Line)){
        Str += Line;
        Str += "\n";
      }
      ScriptMemory_txt.close();//Closes ScriptMemory.
    }


    size_t pos = Str.find(PhysicalGroups->text(),0); //Finds and returns the first Pyisical Groups character position. If was not found, returns npos. pos = position.

    //Reads string Str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = Str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      Str.erase(pos-17, pos_null - pos +19); //Delete the pos and pos_null gap text (-17 because it's necessary to erase the "Physical Group: ").

      /*
      At this moment, Str is modified. Since we want the text to be copied to the ScriptMemory (for future modifications in Str and header),
      besides the script print, the Str is copied to ScriptMemory and Oss(Oss-> buff -> scriptdisplay).
      */
      
      //Firt clear the ScriptMemory
      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
      if(ScriptMemory_txt.is_open()){
        ScriptMemory_txt.close();
      }

      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
      if(ScriptMemory_txt.is_open()){

      ScriptMemory_txt << Str; //Append Str to ScriptMemory
      ScriptMemory_txt << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Area: " << b << "\n"<< "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory.
      }
      Buffer.remove(0,Buffer.length());
      Oss<< Str;
      Oss<<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Area: " << b << "\n"<< "\n";
      Buffer.append(Oss.str().c_str());

    }
    else{
      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
      ScriptMemory_txt << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Area: " << b << "\n"<< "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory.
      
      Buffer.remove(0,Buffer.length());
      Oss<< Str;
      Oss<<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Area: " << b << "\n"<< "\n";
      Buffer.append(Oss.str().c_str());
    }
  }

  //Stokes Variables
  if(name == "Stokes"){

    a = this->density->value();
    b = this->viscosity->value();
    c = NULL;

    wf = this->WeakForm->label();
    pg = this->PhysicalGroups->text();
    bc = "";

  //Read File and find PhysicalGroups->text() first character position

    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory.

    //Copies ScriptMemory text to Str.
        if(ScriptMemory_txt.is_open()){ //Opens ScriptMemory.
      while(getline(ScriptMemory_txt, Line)){
        Str += Line;
        Str += "\n";
      }
      ScriptMemory_txt.close();//Closes ScriptMemory.
    }

    size_t pos = Str.find(PhysicalGroups->text(),0); //Finds and returns the first Pyisical Groups character position. If was not found, returns npos. pos = position.

    //Reads string Str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = Str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      Str.erase(pos-17, pos_null - pos +19); //Delete the pos and pos_null gap text (-17 because it's necessary to erase the "Physical Group: ").

      /*
      At this moment, Str is modified. Since we want the text to be copied to the ScriptMemory (for future modifications in Str and header),
      besides the script print, the Str is copied to ScriptMemory and Oss(Oss-> buff -> scriptdisplay).
      */
      
      //Firt clear the ScriptMemory
      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
      if(ScriptMemory_txt.is_open()){
        ScriptMemory_txt.close();
      }

      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
      if(ScriptMemory_txt.is_open()){

      ScriptMemory_txt << Str; //Append Str to ScriptMemory
      ScriptMemory_txt << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Density: " << a << "\n" << "Viscosity: " << b << "\n"<< "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory.
      }
      Buffer.remove(0,Buffer.length());
      Oss<< Str;
      Oss<<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Density: " << a << "\n" << "Viscosity: " << b << "\n"<< "\n";
      Buffer.append(Oss.str().c_str());

    }
    else{
      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
      ScriptMemory_txt << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Density: " << a << "\n" << "Viscosity: " << b << "\n"<< "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory.
      
      Buffer.remove(0,Buffer.length());
      Oss<< Str;
      Oss<<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Density: " << a << "\n" << "Viscosity: " << b << "\n"<< "\n";
      Buffer.append(Oss.str().c_str());
    }
  }

  //NavierStokes Variables
  if(name == "Navier Stokes"){

    a = this->density->value();
    b = this->viscosity->value();
    c = NULL;

    wf = this->WeakForm->label();
    pg = this->PhysicalGroups->text();
    bc = "";

  //Read File and find PhysicalGroups->text() first character position

    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory.

    //Copies ScriptMemory text to Str.
    if(ScriptMemory_txt.is_open()){ //Opens ScriptMemory.
      
      while(getline(ScriptMemory_txt, Line)){
        Str += Line;
        Str += "\n";
      }
      ScriptMemory_txt.close();//Closes ScriptMemory.
    }

    size_t pos = Str.find(PhysicalGroups->text(),0); //Finds and returns the first Pyisical Groups character position. If was not found, returns npos. pos = position.

    //Reads string Str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = Str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      Str.erase(pos-17, pos_null - pos +19); //Delete the pos and pos_null gap text (-17 because it's necessary to erase the "Physical Group: ").

      /*
      At this moment, Str is modified. Since we want the text to be copied to the ScriptMemory (for future modifications in Str and header),
      besides the script print, the Str is copied to ScriptMemory and Oss(Oss-> buff -> scriptdisplay).
      */
      
      //Firt clear the ScriptMemory
      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
      if(ScriptMemory_txt.is_open()){
        ScriptMemory_txt.close();
      }

      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
      if(ScriptMemory_txt.is_open()){

      ScriptMemory_txt << Str; //Append Str to ScriptMemory
      ScriptMemory_txt << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Density: " << a << "\n" << "Viscosity: " << b << "\n"<< "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory.
      }
      Buffer.remove(0,Buffer.length());
      Oss<< Str;
      Oss<<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Density: " << a << "\n" << "Viscosity: " << b << "\n"<< "\n";
      Buffer.append(Oss.str().c_str());

    }
    else{
      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
      ScriptMemory_txt << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Density: " << a << "\n" << "Viscosity: " << b << "\n"<< "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory.
      
      Buffer.remove(0,Buffer.length());
      Oss<< Str;
      Oss<<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Density: " << a << "\n" << "Viscosity: " << b << "\n"<< "\n";
      Buffer.append(Oss.str().c_str());
    }
  }

  //Poisson Variables
  if(name == "Poisson"){

    a = this->Nsteps->value();
    b = NULL;
    c = NULL;

    wf = this->WeakForm->label();
    pg = this->PhysicalGroups->text();
    bc = "";

  //Read File and find PhysicalGroups->text() first character position

    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory.

    //Copies ScriptMemory text to Str.
    if(ScriptMemory_txt.is_open()){ //Opens ScriptMemory.
      
      while(getline(ScriptMemory_txt, Line)){
        Str += Line;
        Str += "\n";
      }
      ScriptMemory_txt.close();//Closes ScriptMemory.
    }

    size_t pos = Str.find(PhysicalGroups->text(),0); //Finds and returns the first Pyisical Groups character position. If was not found, returns npos. pos = position.

    //Reads string Str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = Str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      Str.erase(pos-17, pos_null - pos +19); //Delete the pos and pos_null gap text (-17 because it's necessary to erase the "Physical Group: ").

      /*
      At this moment, Str is modified. Since we want the text to be copied to the ScriptMemory (for future modifications in Str and header),
      besides the script print, the Str is copied to ScriptMemory and Oss(Oss-> buff -> scriptdisplay).
      */
      
      //Firt clear the ScriptMemory
      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
      if(ScriptMemory_txt.is_open()){
        ScriptMemory_txt.close();
      }

      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
      if(ScriptMemory_txt.is_open()){

      ScriptMemory_txt << Str; //Append Str to ScriptMemory
      ScriptMemory_txt << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Number of State Variables: " << a << "\n"<< "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory.
      }
      Buffer.remove(0,Buffer.length());
      Oss<< Str;
      Oss<<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Number of State Variables: " << a << "\n"<< "\n";
      Buffer.append(Oss.str().c_str());

    }
    else{
      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
      ScriptMemory_txt << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Number of State Variables: " << a << "\n"<< "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory.
      
      Buffer.remove(0,Buffer.length());
      Oss<< Str;
      Oss<<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Number of State Variables: " << a << "\n"<< "\n";
      Buffer.append(Oss.str().c_str());
    }
  }

  //L2Projection Variables
  if(name == "L2 Projection"){

    a = this->x->value();
    b = this->y->value();
    c = this->z->value();

    wf = this->WeakForm->label();
    pg = this->PhysicalGroups->text();
    bc = this->BoundaryCondition->text();

   //Read File and find PhysicalGroups->text() first character position

    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory.

    //Copies ScriptMemory text to Str.
    if(ScriptMemory_txt.is_open()){ //Opens ScriptMemory.
      while(getline(ScriptMemory_txt, Line)){
        Str += Line;
        Str += "\n";
      }
      ScriptMemory_txt.close();//Closes ScriptMemory.
    }

    size_t pos = Str.find(PhysicalGroups->text(),0); //Finds and returns the first Pyisical Groups character position. If was not found, returns npos. pos = position.

    //Reads string Str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = Str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      Str.erase(pos-17, pos_null - pos +19); //Delete the pos and pos_null gap text (-17 because it's necessary to erase the "Physical Group: ").

      /*
      At this moment, Str is modified. Since we want the text to be copied to the ScriptMemory (for future modifications in Str and header),
      besides the script print, the Str is copied to ScriptMemory and Oss(Oss-> buff -> scriptdisplay).
      */
      
      //Firt clear the ScriptMemory
      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
      if(ScriptMemory_txt.is_open()){
        ScriptMemory_txt.close();
      }

      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
      if(ScriptMemory_txt.is_open()){

      ScriptMemory_txt << Str; //Append Str to ScriptMemory
      ScriptMemory_txt << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Boundary Condition: " << bc << "\n" <<"X: " << a << "\n" << "Y: " << b << "\n" << "Z: " << c << "\n"<< "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory.
      }
      Buffer.remove(0,Buffer.length());
      Oss<< Str;
      Oss<<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Boundary Condition: " << bc << "\n" <<"X: " << a << "\n" << "Y: " << b << "\n" << "Z: " << c << "\n"<< "\n";
      Buffer.append(Oss.str().c_str());

    }
    else{
      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
      ScriptMemory_txt << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Boundary Condition: " << bc << "\n" <<"X: " << a << "\n" << "Y: " << b << "\n" << "Z: " << c << "\n"<< "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory.
      
      Buffer.remove(0,Buffer.length());
      Oss<< Str;
      Oss<<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Boundary Condition: " << bc << "\n" <<"X: " << a << "\n" << "Y: " << b << "\n" << "Z: " << c << "\n"<< "\n";
      Buffer.append(Oss.str().c_str());
    }
  }
}

/* ========================= fElasticity2D_cb ========================= */

void WindowConstructor::fElasticity2D_cb(){

    this->Elasticity2D_Menu->show();
    this->ElasticityPositional2D_Menu->hide();
    this->Elasticity3D_Menu->hide();
    this->ElasticTruss_Menu->hide();
    this->PositionalTruss_Menu->hide();
    this->Stokes_Menu->hide();
    this->NavierStokes_Menu->hide();
    this->Poisson_Menu->hide();
    this->L2Projection_Menu->hide();
}

/* ========================= fElasticityPositional2D_cb ========================= */

void WindowConstructor::fElasticityPositional2D_cb(){

    this->Elasticity2D_Menu->hide();
    this->ElasticityPositional2D_Menu->show();
    this->Elasticity3D_Menu->hide();
    this->ElasticTruss_Menu->hide();
    this->PositionalTruss_Menu->hide();
    this->Stokes_Menu->hide();
    this->NavierStokes_Menu->hide();
    this->Poisson_Menu->hide();
    this->L2Projection_Menu->hide();
}

/* ========================= fElasticity3D_cb ========================= */

void WindowConstructor::fElasticity3D_cb(){

    this->Elasticity2D_Menu->hide();
    this->ElasticityPositional2D_Menu->hide();
    this->Elasticity3D_Menu->show();
    this->ElasticTruss_Menu->hide();
    this->PositionalTruss_Menu->hide();
    this->Stokes_Menu->hide();
    this->NavierStokes_Menu->hide();
    this->Poisson_Menu->hide();
    this->L2Projection_Menu->hide();
}

/* ========================= fElasticTruss_cb ========================= */

void WindowConstructor::fElasticTruss_cb(){

    this->Elasticity2D_Menu->hide();
    this->ElasticityPositional2D_Menu->hide();
    this->Elasticity3D_Menu->hide();
    this->ElasticTruss_Menu->show();
    this->PositionalTruss_Menu->hide();
    this->Stokes_Menu->hide();
    this->NavierStokes_Menu->hide();
    this->Poisson_Menu->hide();
    this->L2Projection_Menu->hide();
}

/* ========================= fPositionalTruss_cb ========================= */

void WindowConstructor::fPositionalTruss_cb(){

    this->Elasticity2D_Menu->hide();
    this->ElasticityPositional2D_Menu->hide();
    this->Elasticity3D_Menu->hide();
    this->ElasticTruss_Menu->hide();
    this->PositionalTruss_Menu->show();
    this->Stokes_Menu->hide();
    this->NavierStokes_Menu->hide();
    this->Poisson_Menu->hide();
    this->L2Projection_Menu->hide();
}

/* ========================= fStokes_cb ========================= */

void WindowConstructor::fStokes_cb(){

    this->Elasticity2D_Menu->hide();
    this->ElasticityPositional2D_Menu->hide();
    this->Elasticity3D_Menu->hide();
    this->ElasticTruss_Menu->hide();
    this->PositionalTruss_Menu->hide();
    this->Stokes_Menu->show();
    this->NavierStokes_Menu->hide();
    this->Poisson_Menu->hide();
    this->L2Projection_Menu->hide();
}

/* ========================= fNavierStokes_cb ========================= */

void WindowConstructor::fNavierStokes_cb(){

    this->Elasticity2D_Menu->hide();
    this->ElasticityPositional2D_Menu->hide();
    this->Elasticity3D_Menu->hide();
    this->ElasticTruss_Menu->hide();
    this->PositionalTruss_Menu->hide();
    this->Stokes_Menu->hide();
    this->NavierStokes_Menu->show();
    this->Poisson_Menu->hide();
    this->L2Projection_Menu->hide();
}

/* ========================= fPoisson_cb ========================= */

void WindowConstructor::fPoisson_cb(){

    this->Elasticity2D_Menu->hide();
    this->ElasticityPositional2D_Menu->hide();
    this->Elasticity3D_Menu->hide();
    this->ElasticTruss_Menu->hide();
    this->PositionalTruss_Menu->hide();
    this->Stokes_Menu->hide();
    this->NavierStokes_Menu->hide();
    this->Poisson_Menu->show();
    this->L2Projection_Menu->hide();
}

/* ========================= fL2Projection_cb ========================= */

void WindowConstructor::fL2Projection_cb(){

    this->Elasticity2D_Menu->hide();
    this->ElasticityPositional2D_Menu->hide();
    this->Elasticity3D_Menu->hide();
    this->ElasticTruss_Menu->hide();
    this->PositionalTruss_Menu->hide();
    this->Stokes_Menu->hide();
    this->NavierStokes_Menu->hide();
    this->Poisson_Menu->hide();
    this->L2Projection_Menu->show();
}