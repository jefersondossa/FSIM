/* This cpp file describes the following Script Display functions:

1.  fButton_SaveScript_cb
2.  fButton_TextSave_cb
3.  fButton_TextCancel_cb

*/

#include "WindowConstructor.h"

/* ========================= fButton_SaveScript_cb() ========================= */

void WindowConstructor::fButton_SaveScript_cb(){

    this->Script_Popup_Menu->show();

}

/* ========================= fButton_TextSave_cb() ========================= */

void WindowConstructor::fButton_TextSave_cb(){

//**********************************************************************************************
/*This function is separated in 3 parts:

  1. Reads the txt input file name and creates the name.txt;
  2. Reads the txt input file name and creates the name.h;

//**********************************************************************************************

/* ========================= PART 1. ========================= */

  string name = Script_Text_Name->value();
  string name_txt = name + ".txt";
  string name_path_txt = "../GraphicInterface/ScriptFiles/" + name_txt;

  fstream new_script;

  new_script.open(name_path_txt,ios::out); //Allows to edit new_script.
  ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in);//Allow to read ScriptMemory.

  if(ScriptMemory_txt.is_open()){ //Opens ScriptMemory.
      Str.clear();
      while(getline(ScriptMemory_txt, Line)){
        Str += Line;
        Str += "\n";
      }
      ScriptMemory_txt.close();//Closes ScriptMemory.
      if(new_script.is_open()){
        new_script << Str; //Copy text to new_script

        new_script.close();
      }
    }

/* ========================= PART 2. ========================= */

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
  ScriptMemory_txt.open("../GraphicInterface/ScriptFiles/ScriptMemory.txt",ios::in);
  if (ScriptMemory_txt.is_open() && new_header.is_open()) {

    vector<PhysicalGroup> groups;

    new_header << "// Defines the problem dimension" << endl;
    new_header << "const int dimension = " << this->dimension->value() << ";" << endl;
    new_header << endl;
    new_header << "{" << endl;

    new_header << "auto yieldFunction = [](const double &plast, double &sigma_y, double &hardening){" << endl;
    new_header << "hardening = 500.;" << endl;
    new_header << "sigma_y =.5+hardening*plast;" << endl;
    new_header << "};" << endl;
    new_header << endl; 

    new_header << "CompMesh* cmesh = new CompMesh();" << endl;
    new_header << endl;

    smatch match;
    PhysicalGroup group;

    regex GroupRegex(R"(Physical Groups: ([^\n]+))");
    regex WeakFormRegex(R"(Weak Form: (.+))");
    regex BoundaryConditionRegex(R"(Boundary Condition: ([^\n]+))");
    regex PropertyRegex(R"(([^:]+): (-?[\d.]+))");


  while (getline(ScriptMemory_txt, Line)) {
    if (regex_search(Line, match, GroupRegex)) {
      group.Name = match[1];
} else if (regex_search(Line, match, WeakFormRegex)) {
  group.WeakForm = match[1];
} else if (regex_search(Line, match, BoundaryConditionRegex)) {
  group.BoundaryCondition = match[1];                
} else if (regex_search(Line, match, PropertyRegex)) {
  group.Properties[match[1]] = stod(match[2]);
} else if (Line.empty()) {
  groups.push_back(group);
  group = {};
}
}

for (const auto& group : groups) {

if (group.WeakForm == "Elasticity 2D") {

  bool planeStress = group.Properties.at("Plane Stress") == 1;

  new_header << "//Physical Group " << group.Name<< endl;
  new_header << "Elasticity2D * " << group.Name << " = new Elasticity2D(\"matid\", "
    << group.Properties.at("Young Modulus") << ", "
    << group.Properties.at("Poisson Modulus") << ", "
    << (planeStress ? "true" : "false") << ");" << endl << endl;
  new_header << "cmesh->InsertMaterial(" << group.Name << ");" << endl << endl;
  } 

  else if (group.WeakForm == "Elasticity Positional 2D") { 
  
  bool planeStress = group.Properties.at("Plane Stress") == 1;
    
  new_header << "//Physical Group " << group.Name<< endl;
  new_header << "ElasticityPositional2D * " << group.Name << " = new ElasticityPositional2D(\"matid\", "
    << group.Properties.at("Young Modulus") << ", "
    << group.Properties.at("Poisson Modulus") << ", "
    << (planeStress ? "true" : "false") << ");" << endl << endl;
    new_header << "cmesh->InsertMaterial(" << group.Name << ");" << endl << endl;
  }
  else if (group.WeakForm == "Elasticity 3D") {  
  new_header << "//Physical Group " << group.Name<< endl;
  new_header << "Elasticity3D * " << group.Name << " = new Elasticity3D(\"matid\", "
    << group.Properties.at("Young Modulus") << ", "
    << group.Properties.at("Poisson Modulus") << ");" << endl << endl;
  new_header << "cmesh->InsertMaterial(" << group.Name << ");" << endl << endl;
  }  
  else if (group.WeakForm == "Elastic Truss") {
  new_header << "//Physical Group " << group.Name<< endl;
  new_header << "ElasticTruss * " << group.Name << " = new ElasticTruss(\"matid\", \"dimension\", "
    << group.Properties.at("Young Modulus") << ", "
    << group.Properties.at("Area") << ");" << endl << endl;
  new_header << "cmesh->InsertMaterial(" << group.Name << ");" << endl << endl;
  }
  else if (group.WeakForm == "Positional Truss") {
  new_header << "//Physical Group " << group.Name<< endl;
  new_header << "PositionalTruss * " << group.Name << " = new PositionalTruss(\"matid\", \"dimension\", "
    << group.Properties.at("Young Modulus") << ", "
    << group.Properties.at("Area") << ");" << endl << endl;
  new_header << "cmesh->InsertMaterial(" << group.Name << ");" << endl << endl;
  }    
  else if (group.WeakForm == "Stokes") {
  new_header << "//Physical Group " << group.Name<< endl;
  new_header << "Stokes * " << group.Name << " = new Stokes(\"matid\", \"dimension\", "
    << group.Properties.at("Density") << ", "
    << group.Properties.at("Viscosity") << ");" << endl << endl;
  new_header << "cmesh->InsertMaterial(" << group.Name << ");" << endl << endl;
  }
  else if (group.WeakForm == "Navier Stokes") {
  new_header << "//Physical Group " << group.Name<< endl;
  new_header << "NavierStokes * " << group.Name << " = new NavierStokes(\"matid\", \"dimension\", "
    << group.Properties.at("Density") << ", "
    << group.Properties.at("Viscosity") << ");" << endl << endl;
  new_header << "cmesh->InsertMaterial(" << group.Name << ");" << endl << endl;
  }
  else if (group.WeakForm == "Poisson") {
  new_header << "//Physical Group " << group.Name<< endl;
  new_header << "Poisson * " << group.Name << " = new Poisson(\"matid\", \"dimension\", "
    << group.Properties.at("Number of State Variables") << ");" << endl << endl;
  new_header << "cmesh->InsertMaterial(" << group.Name << ");" << endl << endl;
  }
  else if (group.WeakForm == "L2 Projection") {
  new_header << "//Physical Group " << group.Name<< endl;
  new_header << "L2Projection * " << group.Name << " = new L2Projection(\"matid\", \"dimension\", "
    << group.BoundaryCondition << ", "
    << group.Properties.at("X") << ", "
    << group.Properties.at("Y") << ", "
    << group.Properties.at("Z") << ");" << endl << endl;
  new_header << "cmesh->InsertMaterial(" << group.Name << ");" << endl << endl;
  }
}

    new_header << "};" << endl;
}
  ScriptMemory_txt.close(); // Fecha ScriptMemory.
  new_header.close();   // Fecha o novo arquivo de script.


this->Script_Popup_Menu->hide();

}

/* ========================= fButton_TextCancel_cb() ========================= */

void WindowConstructor::fButton_TextCancel_cb(){
    this->Script_Popup_Menu->hide();
}