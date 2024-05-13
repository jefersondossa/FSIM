/*
This Implem_WC.cpp file describes the Container_1 created in WindowConstructor.h. The Window is divided in Containers, being:

-> Container_1: Contains the geo. and msh. files input and output (File Path (Container_3)). Besides, reads and interprets the geo. file,
granting the user, in the Container_2, to choose all analysis type and variant. 

-> Container_2: Allows the user to input the desired variants values. This values can be verified at Script (Container_3).

-> Container_3: It's a text display which allows the user to see the entries that he has placed.

-> Container_4: Allerts erros for the user.

*/

#include "WindowConstructor.h"
#include <string>
#include <iostream>
#include <fstream>
#include <regex>

using namespace std;

Fl_Text_Buffer *buff = new Fl_Text_Buffer; //buff for weakform variables script;

//Define a txt file where will be all Script Memory information
fstream ScriptMemory; //fstream library: out to write mode (overwrites), app (append), in (read)


//_________________________________________________________________________________
//Container_2 and Container_3

//Shows Container_2_1 and hides Containers_2_i;
void WindowConstructor::WF1_f() {

      this->Container_2_1->show();
      this->Container_2_2->hide();
      this->Container_2_3->hide();
      this->Container_2_4->hide();
      this->Container_2_5->hide();
      this->Container_2_6->hide();
      this->Container_2_7->hide();
      this->Container_2_8->hide();
      this->Container_2_9->hide();
}

//Shows Container_2_2 and hides Containers_2_i;
void WindowConstructor::WF2_f() {

      this->Container_2_2->show();
      this->Container_2_1->hide();
      this->Container_2_3->hide();
      this->Container_2_4->hide();
      this->Container_2_5->hide();
      this->Container_2_6->hide();
      this->Container_2_7->hide();
      this->Container_2_8->hide();
      this->Container_2_9->hide();
}

//Shows Container_2_3 and hides Containers_2_i;
void WindowConstructor::WF3_f() {

      this->Container_2_3->show();
      this->Container_2_2->hide();
      this->Container_2_1->hide();
      this->Container_2_4->hide();
      this->Container_2_5->hide();
      this->Container_2_6->hide();
      this->Container_2_7->hide();
      this->Container_2_8->hide();
      this->Container_2_9->hide();
}

//Shows Container_2_4 and hides Containers_2_i;
void WindowConstructor::WF4_f() {

      this->Container_2_4->show();
      this->Container_2_2->hide();
      this->Container_2_3->hide();
      this->Container_2_1->hide();
      this->Container_2_5->hide();
      this->Container_2_6->hide();
      this->Container_2_7->hide();
      this->Container_2_8->hide();
      this->Container_2_9->hide();
}

//Shows Container_2_5 and hides Containers_2_i;
void WindowConstructor::WF5_f() {

      this->Container_2_5->show();
      this->Container_2_2->hide();
      this->Container_2_3->hide();
      this->Container_2_4->hide();
      this->Container_2_1->hide();
      this->Container_2_6->hide();
      this->Container_2_7->hide();
      this->Container_2_8->hide();
      this->Container_2_9->hide();
}

//Shows Container_2_6 and hides Containers_2_i;
void WindowConstructor::WF6_f() {

      this->Container_2_6->show();
      this->Container_2_2->hide();
      this->Container_2_3->hide();
      this->Container_2_4->hide();
      this->Container_2_5->hide();
      this->Container_2_1->hide();
      this->Container_2_7->hide();
      this->Container_2_8->hide();
      this->Container_2_9->hide();
}

//Shows Container_2_7 and hides Containers_2_i;
void WindowConstructor::WF7_f() {

      this->Container_2_7->show();
      this->Container_2_2->hide();
      this->Container_2_3->hide();
      this->Container_2_4->hide();
      this->Container_2_5->hide();
      this->Container_2_6->hide();
      this->Container_2_1->hide();
      this->Container_2_8->hide();
      this->Container_2_9->hide();
}

//Shows Container_2_8 and hides Containers_2_i;
void WindowConstructor::WF8_f() {
      
      this->Container_2_8->show();
      this->Container_2_2->hide();
      this->Container_2_3->hide();
      this->Container_2_4->hide();
      this->Container_2_5->hide();
      this->Container_2_6->hide();
      this->Container_2_7->hide();
      this->Container_2_1->hide();
      this->Container_2_9->hide();
}

//Shows Container_2_9 and hides Containers_2_i;
void WindowConstructor::WF9_f() {

      this->Container_2_9->show();
      this->Container_2_2->hide();
      this->Container_2_3->hide();
      this->Container_2_4->hide();
      this->Container_2_5->hide();
      this->Container_2_6->hide();
      this->Container_2_7->hide();
      this->Container_2_8->hide();
      this->Container_2_1->hide();
}

void WindowConstructor::Apply_2_f(){

  double a; //Object that saves Container_2_i input variable;
  double b; //Object that saves Container_2_i input variable;
  double c; //Object that saves Container_2_i input variable;
  bool d; //Object that saves Container_2_i input variable;

  string wf; //Object that saves Weak Form choice from Container_2;
  string pg; //Object that saves Physical Surface choice from Container_2;
  string bc; //Object that saves boundary COndition choice from Container_2_9;

  scriptdisplay->buffer(buff);
  std::ostringstream oss;

  string str; //string that saves a copie from ScriptMemory.

  //Elasticity2D Variables
  if(this->Container_2_1->visible() != 0){

    //Variables Values
    a = this->E2D_young->value();
    b = this->E2D_poisson->value();
    c = NULL;
    d = this->E2D_PlaneStress->value();

    wf = this->Elasticity2D->label();
    pg = this->PhysicalGroups->text();
    bc = "";

    //Read File and find PhysicalGroups->text() first character position

    ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory.

    //Copies ScriptMemory text to str.
    if(ScriptMemory.is_open()){ //Opens ScriptMemory.
      string line;
      while(getline(ScriptMemory, line)){
        str += line;
        str += "\n";
      }
      ScriptMemory.close();//Closes ScriptMemory.
    }

    size_t pos = str.find(PhysicalGroups->text(),0); //Finds and returns the first Pyisical Groups character position. If was not found, returns npos. pos = position.

    //Reads string str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      str.erase(pos-17, pos_null - pos +19); //Delete the pos and pos_null gap text (-17 because it's necessary to erase the "Physical Group: ").

      /*
      At this moment, str is modified. Since we want the text to be copied to the ScriptMemory (for future modifications in str and header),
      besides the script print, the str is copied to ScriptMemory and oss (oss -> buff -> scriptdisplay).
      */
      
      //Firt clear the ScriptMemory
      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
      if(ScriptMemory.is_open()){
        ScriptMemory.close();
      }

      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
      if(ScriptMemory.is_open()){

      ScriptMemory << str; //Append str to ScriptMemory
      ScriptMemory << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n" << "Plane Stress: " << d << "\n"<< "\n";
      ScriptMemory.close(); //Closes ScriptMemory.
      }
      buff->remove(0,buff->length());
      oss << str;
      oss <<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n" << "Plane Stress: " << d << "\n"<< "\n";
      buff->append(oss.str().c_str());

    }
    else{
      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
      ScriptMemory << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n" << "Plane Stress: " << d << "\n"<< "\n";
      ScriptMemory.close(); //Closes ScriptMemory.
      
      buff->remove(0,buff->length());
      oss << str;
      oss <<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n" << "Plane Stress: " << d << "\n"<< "\n";
      buff->append(oss.str().c_str());
    }
  }

  //ElasticityPositional2D Variables
  if(this->Container_2_2->visible() != 0){

    //Variables Values
    a = this->EP2D_young->value();
    b = this->EP2D_poisson->value();
    c = NULL;
    d = this->EP2D_PlaneStress->value();

    wf = this->ElasticityPositional2D->label();
    pg = this->PhysicalGroups->text();
    bc = "";

    //Read File and find PhysicalGroups->text() first character position

    ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory.

    //Copies ScriptMemory text to str.
    if(ScriptMemory.is_open()){ //Opens ScriptMemory.
      string line;
      while(getline(ScriptMemory, line)){
        str += line;
        str += "\n";
      }
      ScriptMemory.close();//Closes ScriptMemory.
    }

    size_t pos = str.find(PhysicalGroups->text(),0); //Finds and returns the first Pyisical Groups character position. If was not found, returns npos. pos = position.

    //Reads string str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      str.erase(pos-17, pos_null - pos +19); //Delete the pos and pos_null gap text (-17 because it's necessary to erase the "Physical Group: ").

      /*
      At this moment, str is modified. Since we want the text to be copied to the ScriptMemory (for future modifications in str and header),
      besides the script print, the str is copied to ScriptMemory and oss (oss -> buff -> scriptdisplay).
      */
      
      //Firt clear the ScriptMemory
      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
      if(ScriptMemory.is_open()){
        ScriptMemory.close();
      }

      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
      if(ScriptMemory.is_open()){

      ScriptMemory << str; //Append str to ScriptMemory
      ScriptMemory << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n" << "Plane Stress: " << d << "\n"<< "\n";
      ScriptMemory.close(); //Closes ScriptMemory.
      }
      buff->remove(0,buff->length());
      oss << str;
      oss <<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n" << "Plane Stress: " << d << "\n"<< "\n";
      buff->append(oss.str().c_str());

    }
    else{
      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
      ScriptMemory << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n" << "Plane Stress: " << d << "\n"<< "\n";
      ScriptMemory.close(); //Closes ScriptMemory.
      
      buff->remove(0,buff->length());
      oss << str;
      oss <<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n" << "Plane Stress: " << d << "\n"<< "\n";
      buff->append(oss.str().c_str());
    }
  }

  //Elasticity3D Variables
  if(this->Container_2_3->visible() != 0){

    a = this->E3D_young->value();
    b = this->E3D_poisson->value();
    c = NULL;
    d = NULL;

    wf = this->Elasticity3D->label();
    pg = this->PhysicalGroups->text();
    bc = "";

  //Read File and find PhysicalGroups->text() first character position

    ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory.

    //Copies ScriptMemory text to str.
    if(ScriptMemory.is_open()){ //Opens ScriptMemory.
      string line;
      while(getline(ScriptMemory, line)){
        str += line;
        str += "\n";
      }
      ScriptMemory.close();//Closes ScriptMemory.
    }

    size_t pos = str.find(PhysicalGroups->text(),0); //Finds and returns the first Pyisical Groups character position. If was not found, returns npos. pos = position.

    //Reads string str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      str.erase(pos-17, pos_null - pos +19); //Delete the pos and pos_null gap text (-17 because it's necessary to erase the "Physical Group: ").

      /*
      At this moment, str is modified. Since we want the text to be copied to the ScriptMemory (for future modifications in str and header),
      besides the script print, the str is copied to ScriptMemory and oss (oss -> buff -> scriptdisplay).
      */
      
      //Firt clear the ScriptMemory
      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
      if(ScriptMemory.is_open()){
        ScriptMemory.close();
      }

      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
      if(ScriptMemory.is_open()){

      ScriptMemory << str; //Append str to ScriptMemory
      ScriptMemory << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n"<< "\n";
      ScriptMemory.close(); //Closes ScriptMemory.
      }
      buff->remove(0,buff->length());
      oss << str;
      oss <<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n"<< "\n";
      buff->append(oss.str().c_str());

    }
    else{
      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
      ScriptMemory << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n"<< "\n";
      ScriptMemory.close(); //Closes ScriptMemory.
      
      buff->remove(0,buff->length());
      oss << str;
      oss <<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n"<< "\n";
      buff->append(oss.str().c_str());
    }
  }

  //ElasticTruss Variables
  if(this->Container_2_4->visible() != 0){

    a = this->ET_young->value();
    b = this->ET_area->value(); 
    c = NULL;
    d = NULL;

    wf = this->ElasticTruss->label();
    pg = this->PhysicalGroups->text();
    bc = "";

   //Read File and find PhysicalGroups->text() first character position

    ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory.

    //Copies ScriptMemory text to str.
    if(ScriptMemory.is_open()){ //Opens ScriptMemory.
      string line;
      while(getline(ScriptMemory, line)){
        str += line;
        str += "\n";
      }
      ScriptMemory.close();//Closes ScriptMemory.
    }

    size_t pos = str.find(PhysicalGroups->text(),0); //Finds and returns the first Pyisical Groups character position. If was not found, returns npos. pos = position.

    //Reads string str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      str.erase(pos-17, pos_null - pos +19); //Delete the pos and pos_null gap text (-17 because it's necessary to erase the "Physical Group: ").

      /*
      At this moment, str is modified. Since we want the text to be copied to the ScriptMemory (for future modifications in str and header),
      besides the script print, the str is copied to ScriptMemory and oss (oss -> buff -> scriptdisplay).
      */
      
      //Firt clear the ScriptMemory
      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
      if(ScriptMemory.is_open()){
        ScriptMemory.close();
      }

      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
      if(ScriptMemory.is_open()){

      ScriptMemory << str; //Append str to ScriptMemory
      ScriptMemory << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Area: " << b << "\n"<< "\n";
      ScriptMemory.close(); //Closes ScriptMemory.
      }
      buff->remove(0,buff->length());
      oss << str;
      oss <<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Area: " << b << "\n"<< "\n";
      buff->append(oss.str().c_str());

    }
    else{
      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
      ScriptMemory << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Area: " << b << "\n"<< "\n";
      ScriptMemory.close(); //Closes ScriptMemory.
      
      buff->remove(0,buff->length());
      oss << str;
      oss <<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Area: " << b << "\n"<< "\n";
      buff->append(oss.str().c_str());
    }
  }  
  

  //PositionalTruss Variables
  if(this->Container_2_5->visible() != 0){

    a = this->PT_young->value();
    b = this->PT_area->value();
    c = NULL;
    d = NULL;

    wf = this->PositionalTruss->label();
    pg = this->PhysicalGroups->text();
    bc = "";

  //Read File and find PhysicalGroups->text() first character position

    ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory.

    //Copies ScriptMemory text to str.
    if(ScriptMemory.is_open()){ //Opens ScriptMemory.
      string line;
      while(getline(ScriptMemory, line)){
        str += line;
        str += "\n";
      }
      ScriptMemory.close();//Closes ScriptMemory.
    }

    size_t pos = str.find(PhysicalGroups->text(),0); //Finds and returns the first Pyisical Groups character position. If was not found, returns npos. pos = position.

    //Reads string str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      str.erase(pos-17, pos_null - pos +19); //Delete the pos and pos_null gap text (-17 because it's necessary to erase the "Physical Group: ").

      /*
      At this moment, str is modified. Since we want the text to be copied to the ScriptMemory (for future modifications in str and header),
      besides the script print, the str is copied to ScriptMemory and oss (oss -> buff -> scriptdisplay).
      */
      
      //Firt clear the ScriptMemory
      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
      if(ScriptMemory.is_open()){
        ScriptMemory.close();
      }

      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
      if(ScriptMemory.is_open()){

      ScriptMemory << str; //Append str to ScriptMemory
      ScriptMemory << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Area: " << b << "\n"<< "\n";
      ScriptMemory.close(); //Closes ScriptMemory.
      }
      buff->remove(0,buff->length());
      oss << str;
      oss <<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Area: " << b << "\n"<< "\n";
      buff->append(oss.str().c_str());

    }
    else{
      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
      ScriptMemory << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Area: " << b << "\n"<< "\n";
      ScriptMemory.close(); //Closes ScriptMemory.
      
      buff->remove(0,buff->length());
      oss << str;
      oss <<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Young Modulus: " << a << "\n" << "Area: " << b << "\n"<< "\n";
      buff->append(oss.str().c_str());
    }
  }

  //Stokes Variables
  if(this->Container_2_6->visible() != 0){

    a = this->S_density->value();
    b = this->S_viscosity->value();
    c = NULL;
    d = NULL;

    wf = this->Stokes->label();
    pg = this->PhysicalGroups->text();
    bc = "";

  //Read File and find PhysicalGroups->text() first character position

    ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory.

    //Copies ScriptMemory text to str.
    if(ScriptMemory.is_open()){ //Opens ScriptMemory.
      string line;
      while(getline(ScriptMemory, line)){
        str += line;
        str += "\n";
      }
      ScriptMemory.close();//Closes ScriptMemory.
    }

    size_t pos = str.find(PhysicalGroups->text(),0); //Finds and returns the first Pyisical Groups character position. If was not found, returns npos. pos = position.

    //Reads string str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      str.erase(pos-17, pos_null - pos +19); //Delete the pos and pos_null gap text (-17 because it's necessary to erase the "Physical Group: ").

      /*
      At this moment, str is modified. Since we want the text to be copied to the ScriptMemory (for future modifications in str and header),
      besides the script print, the str is copied to ScriptMemory and oss (oss -> buff -> scriptdisplay).
      */
      
      //Firt clear the ScriptMemory
      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
      if(ScriptMemory.is_open()){
        ScriptMemory.close();
      }

      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
      if(ScriptMemory.is_open()){

      ScriptMemory << str; //Append str to ScriptMemory
      ScriptMemory << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Density: " << a << "\n" << "Viscosity: " << b << "\n"<< "\n";
      ScriptMemory.close(); //Closes ScriptMemory.
      }
      buff->remove(0,buff->length());
      oss << str;
      oss <<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Density: " << a << "\n" << "Viscosity: " << b << "\n"<< "\n";
      buff->append(oss.str().c_str());

    }
    else{
      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
      ScriptMemory << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Density: " << a << "\n" << "Viscosity: " << b << "\n"<< "\n";
      ScriptMemory.close(); //Closes ScriptMemory.
      
      buff->remove(0,buff->length());
      oss << str;
      oss <<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Density: " << a << "\n" << "Viscosity: " << b << "\n"<< "\n";
      buff->append(oss.str().c_str());
    }
  }

  //NavierStokes Variables
  if(this->Container_2_7->visible() != 0){

    a = this->NS_density->value();
    b = this->NS_viscosity->value();
    c = NULL;
    d = NULL;

    wf = this->NavierStokes->label();
    pg = this->PhysicalGroups->text();
    bc = "";

  //Read File and find PhysicalGroups->text() first character position

    ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory.

    //Copies ScriptMemory text to str.
    if(ScriptMemory.is_open()){ //Opens ScriptMemory.
      string line;
      while(getline(ScriptMemory, line)){
        str += line;
        str += "\n";
      }
      ScriptMemory.close();//Closes ScriptMemory.
    }

    size_t pos = str.find(PhysicalGroups->text(),0); //Finds and returns the first Pyisical Groups character position. If was not found, returns npos. pos = position.

    //Reads string str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      str.erase(pos-17, pos_null - pos +19); //Delete the pos and pos_null gap text (-17 because it's necessary to erase the "Physical Group: ").

      /*
      At this moment, str is modified. Since we want the text to be copied to the ScriptMemory (for future modifications in str and header),
      besides the script print, the str is copied to ScriptMemory and oss (oss -> buff -> scriptdisplay).
      */
      
      //Firt clear the ScriptMemory
      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
      if(ScriptMemory.is_open()){
        ScriptMemory.close();
      }

      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
      if(ScriptMemory.is_open()){

      ScriptMemory << str; //Append str to ScriptMemory
      ScriptMemory << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Density: " << a << "\n" << "Viscosity: " << b << "\n"<< "\n";
      ScriptMemory.close(); //Closes ScriptMemory.
      }
      buff->remove(0,buff->length());
      oss << str;
      oss <<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Density: " << a << "\n" << "Viscosity: " << b << "\n"<< "\n";
      buff->append(oss.str().c_str());

    }
    else{
      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
      ScriptMemory << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Density: " << a << "\n" << "Viscosity: " << b << "\n"<< "\n";
      ScriptMemory.close(); //Closes ScriptMemory.
      
      buff->remove(0,buff->length());
      oss << str;
      oss <<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Density: " << a << "\n" << "Viscosity: " << b << "\n"<< "\n";
      buff->append(oss.str().c_str());
    }
  }

  //Poisson Variables
  if(this->Container_2_8->visible() != 0){

    a = this->P_nstate->value();
    b = NULL;
    c = NULL;
    d = NULL;

    wf = this->Poisson->label();
    pg = this->PhysicalGroups->text();
    bc = "";

  //Read File and find PhysicalGroups->text() first character position

    ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory.

    //Copies ScriptMemory text to str.
    if(ScriptMemory.is_open()){ //Opens ScriptMemory.
      string line;
      while(getline(ScriptMemory, line)){
        str += line;
        str += "\n";
      }
      ScriptMemory.close();//Closes ScriptMemory.
    }

    size_t pos = str.find(PhysicalGroups->text(),0); //Finds and returns the first Pyisical Groups character position. If was not found, returns npos. pos = position.

    //Reads string str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      str.erase(pos-17, pos_null - pos +19); //Delete the pos and pos_null gap text (-17 because it's necessary to erase the "Physical Group: ").

      /*
      At this moment, str is modified. Since we want the text to be copied to the ScriptMemory (for future modifications in str and header),
      besides the script print, the str is copied to ScriptMemory and oss (oss -> buff -> scriptdisplay).
      */
      
      //Firt clear the ScriptMemory
      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
      if(ScriptMemory.is_open()){
        ScriptMemory.close();
      }

      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
      if(ScriptMemory.is_open()){

      ScriptMemory << str; //Append str to ScriptMemory
      ScriptMemory << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Number of State Variables: " << a << "\n"<< "\n";
      ScriptMemory.close(); //Closes ScriptMemory.
      }
      buff->remove(0,buff->length());
      oss << str;
      oss <<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Number of State Variables: " << a << "\n"<< "\n";
      buff->append(oss.str().c_str());

    }
    else{
      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
      ScriptMemory << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Number of State Variables: " << a << "\n"<< "\n";
      ScriptMemory.close(); //Closes ScriptMemory.
      
      buff->remove(0,buff->length());
      oss << str;
      oss <<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Number of State Variables: " << a << "\n"<< "\n";
      buff->append(oss.str().c_str());
    }
  }

  //L2Projection Variables
  if(this->Container_2_9->visible() != 0){

    a = this->valx->value();
    b = this->valy->value();
    c = this->valz->value();
    d = NULL;

    wf = this->L2Projection->label();
    pg = this->PhysicalGroups->text();
    bc = this->BoundaryCondition->text();

   //Read File and find PhysicalGroups->text() first character position

    ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in); //Allows to read ScriptMemory.

    //Copies ScriptMemory text to str.
    if(ScriptMemory.is_open()){ //Opens ScriptMemory.
      string line;
      while(getline(ScriptMemory, line)){
        str += line;
        str += "\n";
      }
      ScriptMemory.close();//Closes ScriptMemory.
    }

    size_t pos = str.find(PhysicalGroups->text(),0); //Finds and returns the first Pyisical Groups character position. If was not found, returns npos. pos = position.

    //Reads string str and finds the next NULL line position.
    if(pos != string::npos){
      size_t pos_null = str.find("\n\n",pos); //find and return the NULL line position. If not found, returns npos.
      str.erase(pos-17, pos_null - pos +19); //Delete the pos and pos_null gap text (-17 because it's necessary to erase the "Physical Group: ").

      /*
      At this moment, str is modified. Since we want the text to be copied to the ScriptMemory (for future modifications in str and header),
      besides the script print, the str is copied to ScriptMemory and oss (oss -> buff -> scriptdisplay).
      */
      
      //Firt clear the ScriptMemory
      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::out);//Allows to append text
      if(ScriptMemory.is_open()){
        ScriptMemory.close();
      }

      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);//Allows to append text
      if(ScriptMemory.is_open()){

      ScriptMemory << str; //Append str to ScriptMemory
      ScriptMemory << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Boundary Condition: " << bc << "\n" <<"X: " << a << "\n" << "Y: " << b << "\n" << "Z: " << c << "\n"<< "\n";
      ScriptMemory.close(); //Closes ScriptMemory.
      }
      buff->remove(0,buff->length());
      oss << str;
      oss <<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Boundary Condition: " << bc << "\n" <<"X: " << a << "\n" << "Y: " << b << "\n" << "Z: " << c << "\n"<< "\n";
      buff->append(oss.str().c_str());

    }
    else{
      ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::app);
      ScriptMemory << "Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Boundary Condition: " << bc << "\n" <<"X: " << a << "\n" << "Y: " << b << "\n" << "Z: " << c << "\n"<< "\n";
      ScriptMemory.close(); //Closes ScriptMemory.
      
      buff->remove(0,buff->length());
      oss << str;
      oss <<"Physical Groups: " << pg << "\n" << "Weak Form: " << wf << "\n" << "Boundary Condition: " << bc << "\n" <<"X: " << a << "\n" << "Y: " << b << "\n" << "Z: " << c << "\n"<< "\n";
      buff->append(oss.str().c_str());
    }
  } 
}

//______________________________________________________________________________

//Save Script Button
void WindowConstructor::Save_Script_f(){
  this->Container_3_1->show();
}

//Save Script popup (Container_3_1)
void WindowConstructor::Save_txt_f(){
  
  string str;
  string name = Input_txt->value();
  string name_txt = name + ".txt";
  string name_path_txt = "../GraphicInterface/ScriptFiles/" + name_txt;

  fstream new_script;

  new_script.open(name_path_txt,ios::out); //Allows to edit new_script.
  ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in);//Allow to read ScriptMemory.

  if(ScriptMemory.is_open()){ //Opens ScriptMemory.
      string line;
      while(getline(ScriptMemory, line)){
        str += line;
        str += "\n";
      }
      ScriptMemory.close();//Closes ScriptMemory.
      if(new_script.is_open()){
        new_script << str; //Copy text to new_script

        new_script.close();
      }
    }
  //______________________________________________________
  // new .h Creation

  fstream new_header;

  string name_h = name + ".h";
  string name_path_h = "../GraphicInterface/ScriptFiles/" + name_h;

  new_header.open(name_path_h,ios::out);
  
  struct PhysicalGroup {
    string Name;
    string WeakForm;
    string BoundaryCondition;
    map<string, double> Properties;
  };
  ScriptMemory.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in);
  if (ScriptMemory.is_open() && new_header.is_open()) {

    vector<PhysicalGroup> groups;

    new_header << "// Defines the problem dimension" << endl;
    new_header << "const int dimension = " << "d" << endl;
    new_header << endl;
    new_header << "{" << endl;

    new_header << "auto yieldFunction = [](const double &plast, double &sigma_y, double &hardening){" << endl;
    new_header << "hardening = 500.;" << endl;
    new_header << "sigma_y =.5+hardening*plast;" << endl;
    new_header << "};" << endl;
    new_header << endl; 

    new_header << "CompMesh* cmesh = new CompMesh();" << endl;
    new_header << endl;

    string line;
    smatch match;
    PhysicalGroup group;

    regex GroupRegex(R"(Physical Groups: ([^\n]+))");
    regex WeakFormRegex(R"(Weak Form: (.+))");
    regex BoundaryConditionRegex(R"(Boundary Condition: ([^\n]+))");
    regex PropertyRegex(R"(([^:]+): ([\d.]+))");


  while (getline(ScriptMemory, line)) {
    if (regex_search(line, match, GroupRegex)) {
      group.Name = match[1];
} else if (regex_search(line, match, WeakFormRegex)) {
  group.WeakForm = match[1];
} else if (regex_search(line, match, BoundaryConditionRegex)) {
  group.BoundaryCondition = match[1];                
} else if (regex_search(line, match, PropertyRegex)) {
  group.Properties[match[1]] = stod(match[2]);
} else if (line.empty()) {
  groups.push_back(group);
  group = {};
}
}

for (const auto& group : groups) {
new_header << "//Physical Group " << group.Name << endl;
new_header << endl;

if (group.WeakForm == "Elasticity 2D") {
  new_header << "Elasticity2D * " << group.Name << " = new Elasticity2D(\"matid\", "
    << group.Properties.at("Young Modulus") << ", "
    << group.Properties.at("Poisson Modulus") << ", "
    << group.Properties.at("Plane Stress") << ");" << endl;
    new_header << endl;
    new_header << "cmesh->InsertMaterial(" << "(" << group.Name << ");" << endl;
  } 
  else if (group.WeakForm == "Elasticity Positional 2D") {  
    new_header << "ElasticityPositional2D * " << group.Name << " = new ElasticityPositional2D(\"matid\", "
    << group.Properties.at("Young Modulus") << ", "
    << group.Properties.at("Poisson Modulus") << ", "
    << group.Properties.at("Plane Stress") << ");" << endl;    
    new_header << endl;
    new_header << "cmesh->InsertMaterial(" << "(" << group.Name << ");" << endl;
  }
  else if (group.WeakForm == "Elasticity 3D") {  
    new_header << "Elasticity3D * " << group.Name << " = new Elasticity3D(\"matid\", "
    << group.Properties.at("Young Modulus") << ", "
    << group.Properties.at("Poisson Modulus") << ");" << endl;
    new_header << endl;
    new_header << "cmesh->InsertMaterial(" << "(" << group.Name << ");" << endl;
  }  
  else if (group.WeakForm == "Elastic Truss") {
    new_header << "ElasticTruss * " << group.Name << " = new ElasticTruss(\"matid\", \"dimension\", "
    << group.Properties.at("Young Modulus") << ", "
    << group.Properties.at("Area") << ");" << endl;
    new_header << endl;
    new_header << "cmesh->InsertMaterial(" << "(" << group.Name << ");" << endl << endl;
  }
  else if (group.WeakForm == "Positional Truss") {
    new_header << "PositionalTruss * " << group.Name << " = new PositionalTruss(\"matid\", \"dimension\", "
    << group.Properties.at("Young Modulus") << ", "
    << group.Properties.at("Area") << ");" << endl << endl;
    new_header << "cmesh->InsertMaterial(" << "(" << group.Name << ");" << endl << endl;
  }    
  else if (group.WeakForm == "Stokes") {
  new_header << "Stokes * " << group.Name << " = new Stokes(\"matid\", \"dimension\", "
    << group.Properties.at("Density") << ", "
    << group.Properties.at("Viscosity") << ");" << endl << endl;
    new_header << "cmesh->InsertMaterial(" << "(" << group.Name << ");" << endl << endl;
  }
  else if (group.WeakForm == "Navier Stokes") {
  new_header << "NavierStokes * " << group.Name << " = new NavierStokes(\"matid\", \"dimension\", "
    << group.Properties.at("Density") << ", "
    << group.Properties.at("Viscosity") << ");" << endl << endl;
    new_header << "cmesh->InsertMaterial(" << "(" << group.Name << ");" << endl << endl;
  }
  else if (group.WeakForm == "Poisson") {
  new_header << "Poisson * " << group.Name << " = new Poisson(\"matid\", \"dimension\", "
    << group.Properties.at("Number of State Variables") << ");" << endl << endl;
    new_header << "cmesh->InsertMaterial(" << "(" << group.Name << ");" << endl << endl;
  }
  else if (group.WeakForm == "L2Projection") {
  new_header << "L2Projection * " << group.Name << " = new L2Projection(\"matid\", \"dimension\", "
    << group.BoundaryCondition << ", "
    << group.Properties.at("X") << ", "
    << group.Properties.at("Y") << ", "
    << group.Properties.at("Z") << ");" << endl << endl;
    new_header << "cmesh->InsertMaterial(" << "(" << group.Name << ");" << endl << endl;
  }
}

    new_header << "};" << endl;
}
  ScriptMemory.close(); // Fecha ScriptMemory.
  new_header.close();   // Fecha o novo arquivo de script.


this->Container_3_1->hide();
}

void WindowConstructor::Cancel_f(){
  this->Container_3_1->hide();
}

