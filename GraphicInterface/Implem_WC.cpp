#include "WindowConstructor.h"

#ifdef HAS_GRAPHINTERFACE

#include <FL/Fl_File_Chooser.H>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;
//Functions definitions from WindowConstructor
//Get variables values 

//Open Window
void WindowConstructor::Open_f() {
  this->Window->show();
}

//_________________________________________________________________________________
//Container_1

void WindowConstructor::Input_geo_f() {
  const char * gp = fl_file_chooser(".geo file path", "*.geo*", "", 0);
  
  this->geo_path->value(gp);
}

void WindowConstructor::Input_msh_f() {
  const char * mp = fl_file_chooser(".msh file path", "*.msh*", "", 0);
  
  this->msh_path->value(mp);
}

void WindowConstructor::SavePath_f(){

    string gp = this->geo_path->value();
    string mp = this->msh_path->value();

}
    //To implement dimension output

//_________________________________________________________________________________
//Container_2

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

void WindowConstructor::Apply_2_f(int &WF_id){

  double a;
  double b;
  double c;
  bool d;

  string wf;
  string sm;
  string bc;

  Fl_Text_Buffer *buff = new Fl_Text_Buffer;
  scriptdisplay->buffer(buff);

  std::ostringstream oss;
  // oss <<"Texto bla" << variable...
  //buff->text(oss.str().c_str());

  if(this->Container_2_1->visible() != 0){

    a = this->E2D_young->value();
    b = this->E2D_poisson->value();
    c = NULL;
    d = this->E2D_PlaneStress->value();

    wf = this->Elasticity2D->label();
    sm = this->surfaces->label();
    bc = "";

    oss << "Weak Form: " << wf << ", " << "Surface: " << sm << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n" << "Plane Stress: " << d << "\n";
    buff->text(oss.str().c_str());

}

  if(this->Container_2_2->visible() != 0){

    a = this->EP2D_young->value();
    b = this->EP2D_poisson->value();
    c = NULL;
    d = this->EP2D_PlaneStress->value();

    wf = this->ElasticityPositional2D->label();
    sm = this->surfaces->label();
    bc = "";

    oss << "Weak Form: " << wf << ", " << "Surface: " << sm << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n" << "Plane Stress: " << d << "\n";
    buff->text(oss.str().c_str());

}

  if(this->Container_2_3->visible() != 0){

    a = this->E3D_young->value();
    b = this->E3D_poisson->value();
    c = NULL;
    d = NULL;

    wf = this->Elasticity3D->label();
    sm = this->surfaces->label();
    bc = "";

    oss << "Weak Form: " << wf << ", " << "Surface: " << sm << "\n" << "Young Modulus: " << a << "\n" << "Poisson Modulus: " << b << "\n";
    buff->text(oss.str().c_str());

}
        
  if(this->Container_2_4->visible() != 0){

    a = this->ET_young->value();
    b = this->ET_area->value(); 
    c = NULL;
    d = NULL;

    wf = this->ElasticTruss->label();
    sm = this->surfaces->label();
    bc = "";

    oss << "Weak Form: " << wf << ", " << "Surface: " << sm << "\n" << "Young Modulus: " << a << "\n" << "Area: " << b << "\n";
    buff->text(oss.str().c_str());

}

  if(this->Container_2_5->visible() != 0){

    a = this->PT_young->value();
    b = this->PT_area->value();
    c = NULL;
    d = NULL;

    wf = this->PositionalTruss->label();
    sm = this->surfaces->label();
    bc = "";

    oss << "Weak Form: " << wf << ", " << "Surface: " << sm << "\n" << "Young Modulus: " << a << "\n" << "Area: " << b << "\n";
    buff->text(oss.str().c_str());

}

  if(this->Container_2_6->visible() != 0){

    a = this->S_density->value();
    b = this->S_viscosity->value();
    c = NULL;
    d = NULL;

    wf = this->Stokes->label();
    sm = this->surfaces->label();
    bc = "";

    oss << "Weak Form: " << wf << ", " << "Surface: " << sm << "\n" << "Density: " << a << "\n" << "Viscosity: " << b << "\n";
    buff->text(oss.str().c_str());

}

  if(this->Container_2_7->visible() != 0){

    a = this->NS_density->value();
    b = this->NS_viscosity->value();
    c = NULL;
    d = NULL;

    wf = this->NavierStokes->label();
    sm = this->surfaces->label();
    bc = "";

    oss << "Weak Form: " << wf << ", " << "Surface: " << sm << "\n" << "Density: " << a << "\n" << "Viscosity: " << b << "\n";
    buff->text(oss.str().c_str());

}

  if(this->Container_2_8->visible() != 0){

    a = this->P_nstate->value();
    b = NULL;
    c = NULL;
    d = NULL;

    wf = this->Poisson->label();
    sm = this->surfaces->label();
    bc = "";

    oss << "Weak Form: " << wf << ", " << "Surface: " << sm << "\n" << "Number of State Variables: " << a << "\n";
    buff->text(oss.str().c_str());

}

  if(this->Container_2_9->visible() != 0){

    a = this->valx->value();
    b = this->valy->value();
    c = this->valz->value();
    d = NULL;

    wf = this->L2Projection->label();
    sm = this->surfaces->label();
    bc = this->BounderyCondition->label();
    //analiser ainda a saída do bc

    oss << "Weak Form: " << wf << ", " << "Surface: " << sm << "\n" << "Boundery Condition: " << bc << "X: " << a << "\n" << "Y: " << b << "\n" << "Z: " << c << "\n";
    buff->text(oss.str().c_str());
      
}
}

//_________________________________________________________________________________
//Container_3

void WindowConstructor::Update(){

  Fl_Text_Buffer *buffer = new Fl_Text_Buffer;

}



#endif