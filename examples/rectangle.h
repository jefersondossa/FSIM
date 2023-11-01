//Analytic Solutions

auto exactSolPoisson = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];

    // u[0] = x * x * (x-1.) * y * y * (y-1.);
    // gradU(0,0) = x*(3*x-2.)*(y-1.)*y*y;
    // gradU(1,0) = x*x*(x-1.)*y*(3.*y-2.);
    u[0] = x;
    gradU(0,0) = 1.;
    gradU(1,0) = 0.;
};

auto forcingFunctionPoisson = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    // force[0] = 1.-x*x;
    // force[0] = -2. * (x*x*(1.-3.*y) - (y-1.)*y*y + 3.*x*(y-1.)*y*y + x*x*x*(3.*y-1.));
};

auto exactSolElasticity2D = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];
    auto pi = M_PI;
    // u[0] = cos(pi*x)*sin(2.*pi*y);
    // u[1] = cos(pi*y)*sin(pi*x);
    // gradU(0,0) = -pi*sin(pi*x)*sin(2.*pi*y);
    // gradU(0,1) = pi*cos(pi*x)*cos(pi*y);
    // gradU(1,0) = 2.*pi*cos(pi*x)*cos(2.*pi*y);
    // gradU(1,1) = -pi*sin(pi*x)*sin(pi*y);
    u[0] = x;
    u[1] = 0.;
    gradU(0,0) = 1.;
    gradU(0,1) = 0.;
    gradU(1,0) = 0.;
    gradU(1,1) = 0.;
};

auto forcingFunctionElasticity2D = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    double E=1.;
    double poisson=0.;
    auto pi = M_PI;
    force[0] = -(E*pi*pi*cos(pi*x)*((1.+poisson)*sin(pi*y) + 2.*(3. - 2.*poisson)*sin(2.*pi*y)))/(2.*(-1. + poisson*poisson));
    force[1] = -(E*pi*pi*(-((-3. + poisson)*cos(pi*y)) + 2.*(1. + poisson)*cos(2*pi*y))*sin(pi*x))/(2.*(-1. + poisson*poisson));
};

auto exactSolStokes = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];

    u[0] = -sin(x)*sin(y);
    u[1] = -cos(x)*cos(y);
    u[2] = cos(x)*cos(y);//-sin(2.*pi*x)*cos(2.*pi*y);
    gradU(0,0) = -cos(x)*sin(y);
    gradU(0,1) = -sin(x)*cos(y);
    gradU(1,0) =  sin(x)*cos(y);
    gradU(1,1) =  cos(x)*sin(y);
    gradU(0,2) = -sin(x)*sin(y);
    gradU(1,2) =  cos(x)*cos(y);
};

auto forcingFunctionStokes = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    double visc = 0.01;

    force[0] = -(1+2.*visc)*sin(x)*sin(y);
    force[1] = (1-2.*visc)*cos(x)*cos(y);
};

auto forcingFunctionNavierStokes = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    double visc = 0.01;

    force[0] = -(sin(x)*(cos(x - 2.*y) + cos(x + 2.*y) - 2.*(1. + 2.*visc)*sin(y)))/2.;
    force[1] = -cos(y)*(cos(x) - 2.*visc*cos(x) + cos(x)*cos(x)*sin(y) - sin(x)*sin(x)*sin(y));
};




   // Defines the problem dimension
    const int dimension = 2;
    const int degree = 1;

for (int k = 3; k < 4; k++)
{
//  Create problem variables 
    
    Arlequin   arlequinProblem; 
   //FSI        coupledProblem;  

//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================
        //Coarse
        Geometry* fluid1 = new Geometry(0);

        Point* p0 = fluid1 -> addPoint({ 0.0, 0.0 },1.5,false);
        Point* p1 = fluid1 -> addPoint({ 1.0, 0.0 },1.5,false);
        Point* p2 = fluid1 -> addPoint({ 1.0, 1.0 },1.5,false);
        Point* p3 = fluid1 -> addPoint({ 0.0, 1.0 },1.5,false);
        
        Line* l0 = fluid1 -> addLine({ p0, p1 });
        
        Line* l2 = fluid1 -> addLine({ p1, p2 });
        Line* l3 = fluid1 -> addLine({ p3, p0 }); 
        Line* l1 = fluid1 -> addLine({ p2, p3 });
        LineLoop* ll0 = fluid1->addLineLoop({ l0, l2, l1, l3});
 
        //std::vector<LineLoop*> lin = {ll0, ll1};
        PlaneSurface* s20 = fluid1 -> addPlaneSurface({ll0});
        
        // double h1 = pow(2,k)+1; double v1 = pow(2,k)+1;
        double h1 = 3; double v1 = 2;
        fluid1 -> transfiniteLine({ l0 }, h1);
        fluid1 -> transfiniteLine({ l2 }, v1);
        fluid1 -> transfiniteLine({ l1 }, h1);
        fluid1 -> transfiniteLine({ l3 }, v1);
        
        fluid1->transfiniteSurface({ s20 }, "Left", {p0,p1,p2,p3});
        
        fluid1 -> addBoundaryCondition("NEUMANN", l0, {0.0}, {0.0}, {},  "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l1, {0.0}, {0.0}, {},  "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l2, {0.1}, {0.0}, {},  "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l3, {0.0}, {0.0}, {},  "GLOBAL");
        

        // Fine
        Geometry* fluid2 = new Geometry(1);

        double radius = 1.0;
        double gluing_zone = 1.0;
        double addit_space = 3.0;

        Point* p10 = fluid2 -> addPoint({0.5,0.},0.1,false);
        Point* p11 = fluid2 -> addPoint({1.5,0.},0.1,false);
        Point* p12 = fluid2 -> addPoint({1.5,1.},0.1,false);
        Point* p13 = fluid2 -> addPoint({0.5,1.},0.1,false);

        Line* l10 = fluid2 -> addLine({ p10, p11 });
        Line* l11 = fluid2 -> addLine({ p11, p12 });
        Line* l12 = fluid2 -> addLine({ p12, p13 });
        Line* l13 = fluid2 -> addLine({ p13, p10 });
        
        LineLoop* ll10 = fluid2->addLineLoop({ l10 , l11, l12, l13});

        PlaneSurface* s1 = fluid2 -> addPlaneSurface({ll10});

        fluid2 -> addBoundaryCondition("NEUMANN", l10, {0.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l11, {1.5}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l12, {0.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l13, {0.0}, {0.0}, {},  "GLOBAL");
        
        

        //Transfinite lines 
        int corn = 20; int side = 7;
        // double h2 = pow(2,k)/2+1;
        // double v2 = pow(2,k)+1;
        double h2 = 3;
        double v2 = 2;
        //corners
        fluid2 -> transfiniteLine({ l10 }, h2);
        fluid2 -> transfiniteLine({ l11 }, v2);
        fluid2 -> transfiniteLine({ l12 }, h2);
        fluid2 -> transfiniteLine({ l13 }, v2);
        
        
        //Surfaces
        fluid2->transfiniteSurface({ s1 }, "Left", {p10,p11,p12,p13});
        // fluid2->transfiniteSurface({ s2 }, "Left", {p11,p12,p17,p16});
        
    if (rank == 0){
        
        PartitionOfUnity elType;
        switch (degree){
        case 1:
            elType=T3;
            break;
        case 2:
            elType = T6;
            break;
        case 3:
            elType = T10;
            break;
        
        default:
            break;
        }

        FluidDomain* problem = new FluidDomain(fluid1);
        // problem -> addSurfaceMaterial({ s20,s21,s22,s23 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
        problem -> generateMesh(elType, DELAUNAY, "coarse", "", false, true);

        FluidDomain* problem2 = new FluidDomain(fluid2);
        // // problem2 -> addSurfaceMaterial({ s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
        problem2 -> generateMesh(elType, DELAUNAY, "fine", "", false, true);

        //problem -> readInput("exemplo.msh",0);
	};

//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================

	MPI_Barrier(PETSC_COMM_WORLD);   
    ProblemParameters pParameters;
    // pParameters.ProbType() = ESolidPositional;
    pParameters.ProbType() = EElastic;
    if (pParameters.ProbType() == EPoisson){
        pParameters.setForcingFunction(forcingFunctionPoisson);
        pParameters.setExactSolution(exactSolPoisson);
    } else if (pParameters.ProbType() == EElastic || pParameters.ProbType() == ESolidPositional){
        pParameters.SetElasticity(1.,0.);
        // pParameters.setForcingFunction(forcingFunctionElasticity2D);
        // pParameters.setExactSolution(exactSolElasticity2D);
    } else if (pParameters.ProbType() == EStokes){
        pParameters.SetIncompressibleFluid(0.01,1.);
        pParameters.setTimeStep(1.);//Needed for stabilization parameter
        pParameters.setForcingFunction(forcingFunctionStokes);
        pParameters.setExactSolution(exactSolStokes);
    } else if (pParameters.ProbType() == ENavierStokes){
        pParameters.SetIncompressibleFluid(0.01,1.);
        pParameters.setTimeStep(1.);//Needed for stabilization parameter
        pParameters.setForcingFunction(forcingFunctionNavierStokes);
        pParameters.setExactSolution(exactSolStokes);
    }
    pParameters.Solver() = SolverType::ESuiteSparse;
    pParameters.ArlequinStab() = ArlequinStabType::ENoStab;
    // pParameters.ArlequinStab() = ArlequinStabType::EOption1;
    pParameters.setSpectralRadius(1.);
    pParameters.setArlequinOperatorConstants(1.,0.0);

    CompMesh* coarseModel = new CompMesh(pParameters,dimension,degree);
    CompMesh* fineModel = new CompMesh(pParameters,dimension,degree);  

    GmshTools::MeshReading(fluid1,"coarse.msh",coarseModel);
    // GmshTools::Read(*coarseModel,"../coarse_test.msh");
    GmshTools::MeshReading(fluid2,"fine.msh",fineModel);

    // CompMeshTools::InitialSolution(coarseModel);

    // coarseModel->meshReading(fluid1,"problem_data.txt","coarse.msh","mirror.txt",0);
    // fineModel->meshReading(fluid2,"problem_data.txt","fine.msh","mirror_fine.txt",0);
   // } 
	MPI_Barrier(PETSC_COMM_WORLD);
    GmshTools::BoundaryConstrains(coarseModel);
    GmshTools::BoundaryConstrains(fineModel);

    std::vector<CompMesh *> meshvector(2);
    // meshvector[0] = coarseModel;
    // meshvector[1] = fineModel;
    // Arlequin arl(meshvector);
    // arl.SetUp();

    // LinearAnalysis an(coarseModel,SolverType::ESuiteSparse);
    // an.Run();
    // LinearAnalysis an(arl.MeshVec(),SolverType::ESuiteSparse);
    // NonLinearAnalysis an(arl.MeshVec(),SolverType::ESuiteSparse);
    NonLinearAnalysis an(coarseModel,SolverType::ESuiteSparse);
    an.Run();

    VTUGenerator::PrintResults(coarseModel,"resultCoarse");
    // VTUGenerator::PrintResults(fineModel,"resultFine");
    // VTUGenerator::PrintResults(arl.MeshVec()[2],"resultCoupling");

    VecDouble errors;
    an.PostProcessError(errors);

}           