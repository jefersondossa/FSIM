//This cpp file describes the fElasticTruss function (ElasticTruss Submenu from Weak Form Menu).

#include "WindowConstructor.h"

void WindowConstructor::fElasticTruss(Fl_Group *group){

  ElasticTruss_Menu = new Fl_Group(1175, 175, 360, 175);
  ElasticTruss_Menu->hide();

  young = new Fl_Value_Input(1310, 190, 95, 25, "Elasticity Modulus: ");

  area = new Fl_Value_Input(1310, 235, 95, 25, "Area: ");

  ElasticTruss_Menu->end();

  group->add(ElasticTruss_Menu);

  this->WeakForm->add("Elastic Truss",0,(Fl_Callback*)WindowConstructor::fStatic_ElasticTruss,0,0);

}