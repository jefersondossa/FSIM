

auto exactSol = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];

    u[0] = .5*(x*x-1);
    gradU(0,0) = x;
    gradU(1,0) = 0.;
};

// Defines the problem dimension
const int dimension = 1;
auto forcing = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    force[0] = -1;
};

{

    CompMesh* cmesh = new CompMesh(ApproxType::EHierarquic); 
    cmesh->SetDefaultOrder(3);
    // CompMesh* cmesh = new CompMesh(ApproxType::EIsoparametric); 

    Poisson * matelas = new Poisson(4,1,1);
    matelas->SetForcingFunction(forcing);
    matelas->SetExactSolution(exactSol);

    //BC
    MatrixDouble val1(1,1);
    val1.setZero();
    VecDouble val2(1);
    val2.setZero();
    //Left
    val2[0] = 0.;
    L2Projection * matbc1 = new L2Projection(2,0,BoundaryConditionType::kDirichlet,val1,val2);

    val2[0] = 0.;
    //Right
    L2Projection * matbc3 = new L2Projection(3,0,BoundaryConditionType::kDirichlet,val1,val2);

    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc3);
 
    cmesh->InsertMaterial(matelas);
      
    GmshTools::Read(*cmesh,"../Beam.msh");


    cmesh->Print("cmesh.txt");
    LinearAnalysis an(cmesh,SolverType::ELU);
   
    std::vector<std::string> ScalarNames, VectorNames;
    // VectorNames = {"Solution"};
    ScalarNames = {"Solution"};
    an.Run();
    cmesh->Print("cmeshup.txt");

    VecDouble errors;
    an.PostProcessError(errors);

    VTUGenerator::PrintResults(cmesh,"resultbeam",ScalarNames,VectorNames);    
        
}           