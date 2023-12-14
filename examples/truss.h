
   // Defines the problem dimension
    const int dimension = 2;
{
//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================
   
	MPI_Barrier(PETSC_COMM_WORLD);   

    CompMesh* truss = new CompMesh();
    truss->SetDimension(dimension);
    PositionalTruss * mattruss = new PositionalTruss(4,dimension,1000.,1.);
    // ElasticTruss * mattruss = new ElasticTruss(4,dimension,1000.,1.);
    
    //BC
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();
    
    L2Projection * matbc1 = new L2Projection(5,2,0,val1,val2);
    val2[1] = -5.;
    L2Projection * matbc2 = new L2Projection(6,2,1,val1,val2);
    
    truss->InsertMaterial(matbc1);
    truss->InsertMaterial(matbc2);
    std::vector<L2Projection *> bcIncrement = {matbc2};

    LinearHardening *plastmodel = new LinearHardening(mattruss,111.,4.); 
    truss->InsertMaterial(plastmodel);
    
    GmshTools::Read(*truss,"../truss.msh");
    
    

    
    IncrementalAnalysis an(truss,SolverType::EUmfpack, 2, bcIncrement);
    // NonLinearAnalysis an(truss,SolverType::EUmfpack);
    // NonLinearAnalysis an(arl.MeshVec(),SolverType::EUmfpack,1.e-6,2);
    // NonLinearAnalysis an(coarseModel,SolverType::EUmfpack);
    // an.Run();
    VecDouble increment(2);
    increment.setZero();
    increment[1] = -1.4;
    an.SetIncrement(increment);

    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"Stress"};
    VectorNames = {"Displacement"};
    // an.Run();
    an.Run("truss",ScalarNames,VectorNames);

}           