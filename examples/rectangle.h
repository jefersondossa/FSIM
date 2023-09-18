
    // Defines the problem dimension
    const int dimension = 2;
    const int degree = 2;
 

    //Type definition
    typedef Fluid<dimension,degree>         FluidModel;
    typedef Arlequin<dimension,degree>      Arlequin;
    typedef FSInteraction<dimension,degree> FSI;

for (int k = 3; k < 4; k++)
{
   
 
//  Create problem variables 
    FluidModel coarseModel(ProblemType::EPoisson), fineModel(ProblemType::EPoisson);  
    Arlequin   arlequinProblem; 
   //FSI        coupledProblem;  

//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================
        //Coarse
        Geometry* fluid1 = new Geometry(0);

        Point* p0 = fluid1 -> addPoint({ 0.0, 0.0 },1.5,false);
        Point* p1 = fluid1 -> addPoint({ 1., 0.0 },1.5,false);
        Point* p2 = fluid1 -> addPoint({ 1., 1. },1.5,false);
        Point* p3 = fluid1 -> addPoint({ 0.0, 1. },1.5,false);
        
        Line* l0 = fluid1 -> addLine({ p0, p1 });
        Line* l1 = fluid1 -> addLine({ p1, p2 });
        Line* l2 = fluid1 -> addLine({ p2, p3 });
        Line* l3 = fluid1 -> addLine({ p3, p0 }); 
 
        LineLoop* ll0 = fluid1->addLineLoop({ l0, l1, l2, l3});
 
        //std::vector<LineLoop*> lin = {ll0, ll1};
        PlaneSurface* s20 = fluid1 -> addPlaneSurface({ll0});
        
        double h1 = pow(2,k)+1; double v1 = pow(2,k)+1;
        fluid1 -> transfiniteLine({ l0 }, h1);
        fluid1 -> transfiniteLine({ l1 }, v1);
        fluid1 -> transfiniteLine({ l2 }, h1);
        fluid1 -> transfiniteLine({ l3 }, v1);
        
        fluid1->transfiniteSurface({ s20 }, "Left", {p0,p1,p2,p3});
        
        fluid1 -> addBoundaryCondition("DIRICHLET", l0, {0.0}, {0.0}, {},  "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l1, {0.0}, {0.0}, {},  "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l2, {0.0}, {0.0}, {},  "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l3, {0.0}, {0.0}, {},  "GLOBAL");
        

        // Fine
        Geometry* fluid2 = new Geometry(1);

        double radius = 1.0;
        double gluing_zone = 1.0;
        double addit_space = 3.0;

        Point* p10 = fluid2 -> addPoint({0.5,0.},0.1,false);
        Point* p11 = fluid2 -> addPoint({1.,0.},0.1,false);
        Point* p12 = fluid2 -> addPoint({1.,1.},0.1,false);
        Point* p13 = fluid2 -> addPoint({0.5,1.},0.1,false);

        Line* l10 = fluid2 -> addLine({ p10, p11 });
        Line* l11 = fluid2 -> addLine({ p11, p12 });
        Line* l12 = fluid2 -> addLine({ p12, p13 });
        Line* l13 = fluid2 -> addLine({ p13, p10 });
        
        LineLoop* ll10 = fluid2->addLineLoop({ l10 , l11, l12, l13});

        PlaneSurface* s1 = fluid2 -> addPlaneSurface({ll10});

        fluid2 -> addBoundaryCondition("DIRICHLET", l10, {0.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l11, {0.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l12, {0.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l13, {0.0}, {0.0}, {},  "GLOBAL");
        
        

        //Transfinite lines 
        int corn = 20; int side = 7;
        double h2 = pow(2,k)/2+1;
        double v2 = pow(2,k)+1;
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
    
    coarseModel.meshReading(fluid1,"problem_data.txt","coarse.msh","mirror.txt",0);
    fineModel.meshReading(fluid2,"problem_data.txt","fine.msh","mirror_fine.txt",0);
   // } 
	MPI_Barrier(PETSC_COMM_WORLD);


auto exactSol = [](const VecDouble &coord, double &u, VecDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];
    
    //
    u = x * x * (x-1.) * y * y * (y-1.);
    gradU[0] = x*(3*x-2.)*(y-1.)*y*y;
    gradU[1] = x*x*(x-1.)*y*(3.*y-2.);
};

auto forcingFunction = [](const VecDouble &coord, double &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    
    //
    force = -2. * (x*x*(1.-3.*y) - (y-1.)*y*y + 3.*x*(y-1.)*y*y + x*x*x*(3.*y-1.));

};

    coarseModel.getFluidParameters().setForcingFunctionPoisson(forcingFunction);
    coarseModel.getFluidParameters().setExactSolutionPoisson(exactSol);
    coarseModel.getFluidParameters().setSolver(SolverType::ESuiteSparse);
    fineModel.getFluidParameters().setForcingFunctionPoisson(forcingFunction);
    fineModel.getFluidParameters().setExactSolutionPoisson(exactSol);
    fineModel.getFluidParameters().setSolver(SolverType::ESuiteSparse);

    // coarseModel.solvePoisson();

    // arlequinProblem.setArlequinStabilization(ArlequinStabType::EOption2);
    arlequinProblem.setArlequinStabilization(ArlequinStabType::ENoStab);
    arlequinProblem.setFluidModels(coarseModel, fineModel) ; 

    arlequinProblem.solveArlequinProblem(1, 1.e-7, 2, 0); 
}           