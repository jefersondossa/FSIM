
    // Defines the problem dimension
    const int dimension = 2;
    const int degree = 1;

    //Type definition
    typedef Fluid<dimension,degree>         FluidModel;
    typedef Arlequin<dimension,degree>      Arlequin;
    typedef FSInteraction<dimension,degree> FSI;


 
//  Create problem variables 
    FluidModel coarseModel, fineModel, control;  
    Arlequin   arlequinProblem; 
   //FSI        coupledProblem;  
    
//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================
        //Coarse
        Geometry* fluid1 = new Geometry(0);

        Point* p0 = fluid1 -> addPoint({ 0.0, 0.0 },0.5,false);
        Point* p1 = fluid1 -> addPoint({ 61., 0.0 },0.5,false);
        Point* p2 = fluid1 -> addPoint({ 61., 32. },0.5,false);
        Point* p3 = fluid1 -> addPoint({ 0.0, 32. },0.5,false);
        
        Line* l0 = fluid1 -> addLine({ p0, p1 });
        Line* l1 = fluid1 -> addLine({ p1, p2 });
        Line* l2 = fluid1 -> addLine({ p2, p3 });
        Line* l3 = fluid1 -> addLine({ p3, p0 }); 
 
        LineLoop* ll0 = fluid1->addLineLoop({ l0, l1, l2, l3});
 
        //std::vector<LineLoop*> lin = {ll0, ll1};
        PlaneSurface* s20 = fluid1 -> addPlaneSurface({ll0});
        
        double h1 = 62; double v1 = 33;
        fluid1 -> transfiniteLine({ l0 }, h1);
        fluid1 -> transfiniteLine({ l1 }, v1);
        fluid1 -> transfiniteLine({ l2 }, h1);
        fluid1 -> transfiniteLine({ l3 }, v1);
        
        fluid1->transfiniteSurface({ s20 }, "Left", {p0,p1,p2,p3});
        
        fluid1 -> addBoundaryCondition("DIRICHLET", l0, {}, {0.0}, {},  "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l1, {}, {}, {},  "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l2, {}, {0.0}, {},  "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l3, {0.125}, {0.0}, {},  "GLOBAL");
        

        // Fine
        Geometry* fluid2 = new Geometry(1);

        double radius = 1.0;
        double gluing_zone = 1.0;
        double addit_space = 3.0;

        Point* p10 = fluid2 -> addPoint({16.,16.},0.1,false);
        Point* p11 = fluid2 -> addPoint({16.,16.+radius},0.1,false);
        Point* p12 = fluid2 -> addPoint({16.,16.+radius+addit_space},0.1,false);
        Point* p13 = fluid2 -> addPoint({16.,16.+radius+addit_space+gluing_zone},0.1,false);
        Point* p14 = fluid2 -> addPoint({16.,16.-radius},0.1,false);
        Point* p15 = fluid2 -> addPoint({16.,16.-radius-addit_space},0.1,false);
        Point* p16 = fluid2 -> addPoint({16.,16.-radius-addit_space-gluing_zone},0.1,false);
        
        Line* l10 = fluid2 -> addCircle({p11,p10,p14});
        Line* l11 = fluid2 -> addCircle({p12,p10,p15});
        Line* l12 = fluid2 -> addCircle({p13,p10,p16});
        Line* l13 = fluid2 -> addCircle({p14,p10,p11});
        Line* l14 = fluid2 -> addCircle({p15,p10,p12});
        Line* l15 = fluid2 -> addCircle({p16,p10,p13});

        LineLoop* ll10 = fluid2->addLineLoop({ l10 -> operator-(), l13 -> operator-(), l11, l14});
        LineLoop* ll11 = fluid2->addLineLoop({ l11 -> operator-(), l14 -> operator-(), l12, l15});

        PlaneSurface* s1 = fluid2 -> addPlaneSurface({ll10});
        PlaneSurface* s2 = fluid2 -> addPlaneSurface({ll11});

        fluid2 -> addBoundaryCondition("DIRICHLET", l10, {0.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l13, {0.0}, {0.0}, {},  "GLOBAL");
        
        fluid2 -> addBoundaryCondition("NEUMANN", l11, {}, {}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l14, {}, {}, {},  "GLOBAL");

        fluid2 -> addBoundaryCondition("GLUE", l12, {0}, {0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l15, {0}, {0}, {},  "GLOBAL");
        

        //Transfinite lines 
        int corn = 20; int side = 7;
        //corners
        fluid2 -> transfiniteLine({ l10 }, corn);
        fluid2 -> transfiniteLine({ l11 }, corn);
        fluid2 -> transfiniteLine({ l12 }, corn);
        fluid2 -> transfiniteLine({ l13 }, corn);
        fluid2 -> transfiniteLine({ l14 }, corn);
        fluid2 -> transfiniteLine({ l15 }, corn);
        
        
        //Surfaces
        // fluid2->transfiniteSurface({ s1 }, "Left", {p10,p11,p16,p15});
        // fluid2->transfiniteSurface({ s2 }, "Left", {p11,p12,p17,p16});
        
    if (rank == 0){
  
        FluidDomain* problem = new FluidDomain(fluid1);
        // problem -> addSurfaceMaterial({ s20,s21,s22,s23 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
        problem -> generateMesh(T3, DELAUNAY, "coarse", "", false, true);

        FluidDomain* problem2 = new FluidDomain(fluid2);
        // problem2 -> addSurfaceMaterial({ s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
        problem2 -> generateMesh(T3, DELAUNAY, "fine", "", false, true);

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

    coarseModel.getFluidParameters().setSolver(SolverType::EMumps);
    fineModel.getFluidParameters().setSolver(SolverType::EMumps);


    coarseModel.setProblemType(ProblemType::ENavierStokes);
    fineModel.setProblemType(ProblemType::ENavierStokes);

    arlequinProblem.setFluidModels(coarseModel, fineModel) ; 

    arlequinProblem.solveArlequinProblem(3, 1.e-7, 2, 1); 
           