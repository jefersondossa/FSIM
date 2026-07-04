#include "WindowConstructor.h"

void WindowConstructor::fHeaderFile(){

    fstream new_header;
    string name = Script_Text_Name->value();
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

    new_header << "auto yieldFunction = [](const REAL &plast, REAL &sigma_y, REAL &hardening){" << endl;
    new_header << "hardening = 500.;" << endl;
    new_header << "sigma_y =.5+hardening*plast;" << endl;
    new_header << "};" << endl;
    new_header << endl; 

    new_header << "CompMesh* cmesh = new CompMesh();" << endl;
    new_header << endl;

    PhysicalGroup group;

    regex GroupRegex(R"(Physical Groups: ([^\n]+))");
    regex WeakFormRegex(R"(Weak Form: (.+))");
    regex BoundaryConditionRegex(R"(Boundary Condition: ([^\n]+))");
    regex PropertyRegex(R"(([^:]+): (-?[\d.]+))");


  while (getline(ScriptMemory_txt, Line)) {
    if (regex_search(Line, Match, GroupRegex)) {
      group.Name = Match[1];
} else if (regex_search(Line, Match, WeakFormRegex)) {
  group.WeakForm = Match[1];
} else if (regex_search(Line, Match, BoundaryConditionRegex)) {
  group.BoundaryCondition = Match[1];                
} else if (regex_search(Line, Match, PropertyRegex)) {
  group.Properties[Match[1]] = stod(Match[2]);
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
  ScriptMemory_txt.close();
  new_header.close();

}