/* This cpp file describes the following Weak Form Menu functions:

00. fWF_free_cb();
01. fApply_WF_cb();
02. fElasticity2D_cb(); 
03. fElasticityPositional2D_cb();
04. fElasticity3D_cb(); 
05. fElasticTruss_cb(); 
06. fPositionalTruss_cb(); 
07. fStokes_cb(); 
08. fNavierStokes_cb(); 
09. fPoisson_cb();
10. fL2Projection_cb(); 
11. fOpenFoam_cb();
12. fU_OpenFoam_cb();
13. fp_OpenFoam_cb();

*/

#include "WindowConstructor.h"
#include "WeakFormMenu/Elasticity2DMenu.cpp"

using namespace std;

/* ========================= fWF_free_cb() ========================= */

void WindowConstructor::fWF_free_cb(){

    this->WeakForm->activate();
    //int matid = PhysGroup[this->WeakForm_Menu->text()];
    int matid = PhysGroup[this->PhysicalGroups->text()];
    playback->setInt(matid);

    //Blocks the user to use the OpenFoam for Curve and Volume Physical Groups.

    for (const auto& pair : PhysGroup2){

      int index = WeakForm->find_index("OpenFoam");
      cout << "Name: " << pair.first << ", Phisycal Group Type: " << pair.second << std::endl;

      if(Dimension == 3 && pair.first == PhysicalGroups->text() && pair.second == "Surface"){

        WeakForm->remove(index);
        WeakForm->add("OpenFoam",0,(Fl_Callback*)WindowConstructor::fStatic_OpenFoam,0,0);
        OpenFoam_Menu->activate();
      }

      if(Dimension != 3 || (pair.first == PhysicalGroups->text() && pair.second != "Surface")){

        WeakForm->remove(index);
        WeakForm->add("OpenFoam",0,(Fl_Callback*)WindowConstructor::fStatic_OpenFoam,0,1);
        OpenFoam_Menu->deactivate();
      }
    }
}

/* ========================= fApply_WF_cb ========================= */

void WindowConstructor::fApply_WF_cb(){

    string name = WeakForm->text();

    string Str;
    string Line;
    ostringstream Oss;

    double a; 
    double b; 
    double c; 
    double d; 
    double e;
    double f;
    double g;

    string ea;
    string eb;
    string ec;
    string ed;
    string ee;
    string ef;
    string wf; 
    string pg; 
    string bc; 
    string bc2;
    string st;
    string ft;
 
  ScriptDisplay->buffer(Buffer);

  //Elasticity2D Variables
  if(name == "Elasticity 2D"){

    //Variables Values
    a = this->young->value();
    b = this->poisson->value();
    c = NULL;

    d= this->plane_stress->value();

    wf = this->WeakForm->text();
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

    size_t pos = Str.find(PhysicalGroups->text(),0); //Finds and returns the first Physical Groups character position. If was not found, returns npos. pos = position.

    //Reads string Str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = Str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      Str.erase(pos-16, pos_null - pos +18); //Delete the pos and pos_null gap text (-16 because it's necessary to erase the "Physical Group: ").

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

    wf = this->WeakForm->text();
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

    size_t pos = Str.find(PhysicalGroups->text(),0); //Finds and returns the first Physical Groups character position. If was not found, returns npos. pos = position.

    //Reads string Str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = Str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      Str.erase(pos-16, pos_null - pos +18); //Delete the pos and pos_null gap text (-16 because it's necessary to erase the "Physical Group: ").

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

    wf = this->WeakForm->text();
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


    size_t pos = Str.find(PhysicalGroups->text(),0); //Finds and returns the first Physical Groups character position. If was not found, returns npos. pos = position.

    //Reads string Str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = Str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      Str.erase(pos-16, pos_null - pos +18); //Delete the pos and pos_null gap text (-16 because it's necessary to erase the "Physical Group: ").

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

    wf = this->WeakForm->text();
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

    size_t pos = Str.find(PhysicalGroups->text(),0); //Finds and returns the first Physical Groups character position. If was not found, returns npos. pos = position.

    //Reads string Str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = Str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      Str.erase(pos-16, pos_null - pos +18); //Delete the pos and pos_null gap text (-16 because it's necessary to erase the "Physical Group: ").

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

    wf = this->WeakForm->text();
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


    size_t pos = Str.find(PhysicalGroups->text(),0); //Finds and returns the first Physical Groups character position. If was not found, returns npos. pos = position.

    //Reads string Str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = Str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      Str.erase(pos-16, pos_null - pos +18); //Delete the pos and pos_null gap text (-16 because it's necessary to erase the "Physical Group: ").

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

    wf = this->WeakForm->text();
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

    size_t pos = Str.find(PhysicalGroups->text(),0); //Finds and returns the first Physical Groups character position. If was not found, returns npos. pos = position.

    //Reads string Str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = Str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      Str.erase(pos-16, pos_null - pos +18); //Delete the pos and pos_null gap text (-16 because it's necessary to erase the "Physical Group: ").

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

    wf = this->WeakForm->text();
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

    size_t pos = Str.find(PhysicalGroups->text(),0); //Finds and returns the first Physical Groups character position. If was not found, returns npos. pos = position.

    //Reads string Str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = Str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      Str.erase(pos-16, pos_null - pos +18); //Delete the pos and pos_null gap text (-16 because it's necessary to erase the "Physical Group: ").

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

    wf = this->WeakForm->text();
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

    size_t pos = Str.find(PhysicalGroups->text(),0); //Finds and returns the first Physical Groups character position. If was not found, returns npos. pos = position.

    //Reads string Str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = Str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      Str.erase(pos-16, pos_null - pos +18); //Delete the pos and pos_null gap text (-16 because it's necessary to erase the "Physical Group: ").

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

    wf = this->WeakForm->text();
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

    size_t pos = Str.find(PhysicalGroups->text(),0); //Finds and returns the first Physical Groups character position. If was not found, returns npos. pos = position.

    //Reads string Str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = Str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      Str.erase(pos-16, pos_null - pos +18); //Delete the pos and pos_null gap text (-16 because it's necessary to erase the "Physical Group: ").

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

  //OpenFoam Variables
  if(name == "OpenFoam"){

    ea = to_string(x->value());
    eb = to_string(y->value());
    ec = to_string(z->value());
    ed = to_string(x2->value());

    wf = this->WeakForm->text();
    pg = this->PhysicalGroups->text();
    
    bc = this->BC1_OpenFoam->text();
    bc2 = this->BC2_OpenFoam->text();

    //Stores U and p conditions
    string Uvector;
    string pvector;
    /*_____________________________Condition 1__________________________________*/

    if(bc == "fixedValue"){
      Uvector = "\nUvector: (" + ea + " " + eb + " " + ec + ")";
    }

    else{
      Uvector = "";
    }

    if(bc2 == "fixedValue"){
      pvector = "\npValue: " + ed;
    }

    else{
      pvector = "";
    }

   //Reads File and finds PhysicalGroups->text() first character position

    ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory.

    //Copies ScriptMemory text to Str.
    if(ScriptMemory_txt.is_open()){ //Opens ScriptMemory.
      while(getline(ScriptMemory_txt, Line)){
        Str += Line;
        Str += "\n";
      }
      ScriptMemory_txt.close();//Closes ScriptMemory.
    }

    size_t pos = Str.find(PhysicalGroups->text(),0); //Finds and returns the first Physical Groups character position. If was not found, returns npos. pos = position.

    //Reads string Str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = Str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      Str.erase(pos-16, pos_null - pos +18); //Delete the pos and pos_null gap text (-16 because it's necessary to erase the "Physical Group: ").

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
      ScriptMemory_txt << "Physical Group: " << pg << "\n" << "Weak Form: " << wf << "\n"
                       << "Velocity BC: " << bc << Uvector << "\n" <<"Pressure BC: " << bc2 << pvector << "\n" << "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory.
      }
      Buffer.remove(0,Buffer.length());
      Oss<< Str;
      Oss << "Physical Group: " << pg << "\n" << "Weak Form: " << wf << "\n"
          << "Velocity BC: " << bc << Uvector << "\n" <<"Pressure BC: " << bc2 << pvector << "\n" << "\n";
      Buffer.append(Oss.str().c_str());

    }
    else{
      ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
      ScriptMemory_txt << "Physical Group: " << pg << "\n" << "Weak Form: " << wf << "\n"
                       << "Velocity BC: " << bc << Uvector << "\n" <<"Pressure BC: " << bc2 << pvector << "\n" << "\n";
      ScriptMemory_txt.close(); //Closes ScriptMemory.

      Buffer.remove(0,Buffer.length());
      Oss<< Str;
      Oss << "Physical Group: " << pg << "\n" << "Weak Form: " << wf << "\n"
          << "Velocity BC: " << bc << Uvector << "\n" <<"Pressure BC: " << bc2 << pvector << "\n" << "\n";
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
    this->OpenFoam_Menu->hide();
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
    this->OpenFoam_Menu->hide();
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
    this->OpenFoam_Menu->hide();
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
    this->OpenFoam_Menu->hide();
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
    this->OpenFoam_Menu->hide();
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
    this->OpenFoam_Menu->hide();
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
    this->OpenFoam_Menu->hide();
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
    this->OpenFoam_Menu->hide();
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
    this->OpenFoam_Menu->hide();
}

/* ========================= fOpenFoam_cb ========================= */

void WindowConstructor::fOpenFoam_cb(){

    this->Elasticity2D_Menu->hide();
    this->ElasticityPositional2D_Menu->hide();
    this->Elasticity3D_Menu->hide();
    this->ElasticTruss_Menu->hide();
    this->PositionalTruss_Menu->hide();
    this->Stokes_Menu->hide();
    this->NavierStokes_Menu->hide();
    this->Poisson_Menu->hide();
    this->L2Projection_Menu->hide();
    this->OpenFoam_Menu->show();
}

/* ========================= ffixedValue_cb ========================= */

void WindowConstructor::ffixedValue_cb(){

  string text1 = BC1_OpenFoam->text();
  string text2 = BC2_OpenFoam->text();

   if(text1 != "fixedValue"){

     this->x->deactivate();
     this->y->deactivate();
     this->z->deactivate();
     this->Values1->deactivate();
   }

   if(text1 == "fixedValue"){

     this->x->activate();
     this->y->activate();
     this->z->activate();
     this->Values1->activate();
   }

   if(text2 != "fixedValue"){
   this->x2->deactivate();
   }

   if(text2 == "fixedValue"){
   this->x2->activate();
   }

  //if noSlip or fixedValue --> zeroGradient:
  if(text1 == "noSlip"){
    
    int index1 = BC1_OpenFoam->find_index("noSlip");
    BC1_OpenFoam->value(index1);
    int index2 = BC2_OpenFoam->find_index("zeroGradient");
    BC2_OpenFoam->value(index2);

  }

    if(text1 == "slip"){
    
    int index1 = BC1_OpenFoam->find_index("slip");
    BC1_OpenFoam->value(index1);
    int index2 = BC2_OpenFoam->find_index("zeroGradient");
    BC2_OpenFoam->value(index2);

  }

  if(text1 == "fixedValue"){

    int index1 = BC1_OpenFoam->find_index("fixedValue");
    BC1_OpenFoam->value(index1);
    int index = BC2_OpenFoam->find_index("zeroGradient");
    BC2_OpenFoam->value(index);

  }

  if(text2 == "noSlip"){
    
    int index1 = BC2_OpenFoam->find_index("noSlip");
    BC2_OpenFoam->value(index1);
    int index2 = BC1_OpenFoam->find_index("zeroGradient");
    BC1_OpenFoam->value(index2);

  }

    if(text2 == "slip"){
    
    int index1 = BC2_OpenFoam->find_index("slip");
    BC2_OpenFoam->value(index1);
    int index2 = BC1_OpenFoam->find_index("zeroGradient");
    BC1_OpenFoam->value(index2);

  }

  if(text2 == "fixedValue"){

    int index1 = BC2_OpenFoam->find_index("fixedValue");
    BC2_OpenFoam->value(index1);
    int index = BC1_OpenFoam->find_index("zeroGradient");
    BC1_OpenFoam->value(index);
    const Fl_Menu_Item* items = BC1_OpenFoam->menu();
    const Fl_Menu_Item& item = items[index];

  }

  //if empty --> empty
  if(text2 == "empty"){

    int index1 = BC1_OpenFoam->find_index("empty");
    BC1_OpenFoam->value(index1);
  }

   if(text1 == "empty"){

    int index2 = BC2_OpenFoam->find_index("empty");
    BC2_OpenFoam->value(index2);
  }
}

/* ========================= fParameter_OpenFoam_cb ========================= */

void WindowConstructor::fParameter_OpenFoam_cb(){

  Parameter_OpenFoam_Menu->show();
}

/* ========================= fControlDict_OpenFoam_cb ========================= */

void WindowConstructor::fControlDict_OpenFoam_cb(){

  Parameter_OpenFoam_Menu->hide(); 
}

/* ========================= fU_OpenFoam_cb ========================= */

void WindowConstructor::fU_OpenFoam_cb(){

  U_OpenFoam_Menu->show();
  p_OpenFoam_Menu->hide();
}

/* ========================= fp_OpenFoam_cb ========================= */

void WindowConstructor::fp_OpenFoam_cb(){

  p_OpenFoam_Menu->show();
  U_OpenFoam_Menu->hide();
}



