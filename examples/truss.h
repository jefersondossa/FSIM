
   // Defines the problem dimension
    const int dimension = 1;
{
//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================
   
	MPI_Barrier(PETSC_COMM_WORLD);   

    CompMesh* truss = new CompMesh();
    truss->SetDimension(dimension);
    Poisson * matpoisson = new Poisson(4,dimension);
    truss->InsertMaterial(matpoisson);
    
    //BC
    MatrixDouble val1(1,1);
    val1.setZero();
    VecDouble val2(1);
    val2.setZero();
    val2[0] = 1.5;
    L2Projection * matbc1 = new L2Projection(2,1,0,val1,val2);
    val2.setZero();
    L2Projection * matbc2 = new L2Projection(3,1,0,val1,val2);
    
    truss->InsertMaterial(matbc1);
    truss->InsertMaterial(matbc2);
    
    GmshTools::Read(*truss,"../truss.msh");
    
    // LinearAnalysis an(coarseModel,SolverType::ESuiteSparse);
    // an.Run();
    LinearAnalysis an(truss,SolverType::ESuiteSparse);
    // NonLinearAnalysis an(arl.MeshVec(),SolverType::ESuiteSparse,1.e-6,2);
    // NonLinearAnalysis an(coarseModel,SolverType::ESuiteSparse);
    an.Run();

    VTUGenerator::PrintResults(truss,"resultCoarse");

    VecDouble errors;
    // an.PostProcessError(errors);

}           