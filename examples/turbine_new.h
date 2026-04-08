
    // Defines the problem dimension
    const int dimension = 2;
    const int degree = 2;

    //Type definition
    typedef Fluid<dimension,degree>         FluidModel;
    typedef Arlequin<dimension,degree>      Arlequin;
    typedef FSInteraction<dimension,degree> FSI;


 
//  Create problem variables 
    FluidModel coarseModel, fineModel, control;  
    Arlequin   arlequinProblem; 
    FSI        coupledProblem;  
    
//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================
        //Coarse
        Geometry* fluid1 = new Geometry(0);
        double e = 5;

        double a1 = 2.5;
        double z = 0.;
        double a = 50.;
        double b = 150.;
        double as = 6.;

        Point* p000 = fluid1 -> addPoint({-a,-a},a1,false);
        Point* p001 = fluid1 -> addPoint({-2-e,-a},a1,false);
        Point* p002 = fluid1 -> addPoint({ 1+4*e,-a},a1,false);
        Point* p003 = fluid1 -> addPoint({b, -a},a1,false);
        Point* p004 = fluid1 -> addPoint({  -a,-5},a1,false);
        Point* p005 = fluid1 -> addPoint({-2-e,-5},a1/as,false);
        Point* p006 = fluid1 -> addPoint({ 1+4*e,-5},a1/as,false);
        Point* p007 = fluid1 -> addPoint({   b,-5},a1,false);
        Point* p008 = fluid1 -> addPoint({  -a, 5},a1,false);
        Point* p009 = fluid1 -> addPoint({-2-e, 5},a1/as,false);
        Point* p010 = fluid1 -> addPoint({ 1+4*e, 5},a1/as,false);
        Point* p011 = fluid1 -> addPoint({   b, 5},a1,false);
        Point* p012 = fluid1 -> addPoint({  -a,a},a1,false);
        Point* p013 = fluid1 -> addPoint({-2-e,a},a1,false);
        Point* p014 = fluid1 -> addPoint({ 1+4*e,a},a1,false);
        Point* p015 = fluid1 -> addPoint({   b,a},a1,false);
        

        Line* l001 = fluid1 -> addLine({p000,p001});
        Line* l002 = fluid1 -> addLine({p001,p002});
        Line* l003 = fluid1 -> addLine({p002,p003});
        Line* l004 = fluid1 -> addLine({p000,p004});
        Line* l005 = fluid1 -> addLine({p001,p005});        
        Line* l006 = fluid1 -> addLine({p002,p006});
        Line* l007 = fluid1 -> addLine({p003,p007});
        Line* l008 = fluid1 -> addLine({p004,p005});
        Line* l009 = fluid1 -> addLine({p005,p006});
        Line* l010 = fluid1 -> addLine({p006,p007});
        Line* l011 = fluid1 -> addLine({p004,p008});
        Line* l012 = fluid1 -> addLine({p005,p009});
        Line* l013 = fluid1 -> addLine({p006,p010});
        Line* l014 = fluid1 -> addLine({p007,p011});
        Line* l015 = fluid1 -> addLine({p008,p009});
        Line* l016 = fluid1 -> addLine({p009,p010});
        Line* l017 = fluid1 -> addLine({p010,p011});        
        Line* l018 = fluid1 -> addLine({p008,p012});
        Line* l019 = fluid1 -> addLine({p009,p013});
        Line* l020 = fluid1 -> addLine({p010,p014});
        Line* l021 = fluid1 -> addLine({p011,p015});
        Line* l022 = fluid1 -> addLine({p012,p013});
        Line* l023 = fluid1 -> addLine({p013,p014});
        Line* l024 = fluid1 -> addLine({p014,p015});
       
        LineLoop* ll001 = fluid1 -> addLineLoop({ l001, l005, l008 -> operator-(), l004 -> operator-() });
        LineLoop* ll002 = fluid1 -> addLineLoop({ l002, l006, l009 -> operator-(), l005 -> operator-() });
        LineLoop* ll003 = fluid1 -> addLineLoop({ l003, l007, l010 -> operator-(), l006 -> operator-() });
        LineLoop* ll004 = fluid1 -> addLineLoop({ l008, l012, l015 -> operator-(), l011 -> operator-() });
        LineLoop* ll005 = fluid1 -> addLineLoop({ l009, l013, l016 -> operator-(), l012 -> operator-() });
        LineLoop* ll006 = fluid1 -> addLineLoop({ l010, l014, l017 -> operator-(), l013 -> operator-() });
        LineLoop* ll007 = fluid1 -> addLineLoop({ l015, l019, l022 -> operator-(), l018 -> operator-() });
        LineLoop* ll008 = fluid1 -> addLineLoop({ l016, l020, l023 -> operator-(), l019 -> operator-() });
        LineLoop* ll009 = fluid1 -> addLineLoop({ l017, l021, l024 -> operator-(), l020 -> operator-() });
        
        //std::vector<LineLoop*> lin = {ll0, ll1};
        PlaneSurface* s001 = fluid1 -> addPlaneSurface({ll001});
        PlaneSurface* s002 = fluid1 -> addPlaneSurface({ll002});
        PlaneSurface* s003 = fluid1 -> addPlaneSurface({ll003});
        PlaneSurface* s004 = fluid1 -> addPlaneSurface({ll004});
        PlaneSurface* s005 = fluid1 -> addPlaneSurface({ll005});
        PlaneSurface* s006 = fluid1 -> addPlaneSurface({ll006});
        PlaneSurface* s007 = fluid1 -> addPlaneSurface({ll007});
        PlaneSurface* s008 = fluid1 -> addPlaneSurface({ll008});
        PlaneSurface* s009 = fluid1 -> addPlaneSurface({ll009});

        
        fluid1 -> addBoundaryCondition("NEUMANN", l005, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l006, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l007, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l012, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l013, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l014, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l019, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l020, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l021, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l008, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l009, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l010, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l015, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l016, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l017, {}, {}, {}, "GLOBAL");
        
        fluid1 -> addBoundaryCondition("DIRICHLET", l001, {}, {0.0}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l002, {}, {0.0}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l003, {}, {0.0}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l022, {}, {0.0}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l023, {}, {0.0}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l024, {}, {0.0}, {}, "GLOBAL");

        fluid1 -> addBoundaryCondition("DIRICHLET", l004, {20.0}, {0.0}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l011, {20.0}, {0.0}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l018, {20.0}, {0.0}, {}, "GLOBAL");
        
        
        // Fine
        Geometry* fluid2 = new Geometry(1);

        double a2 = 0.1;
        double a4 = 0.1;
        double a3 = 0.5;
        double dxToGlueZone = 7.;
        double dyToGlueZone = 5.;
        double gzsize = 1;
        double auxradius = 0.01;

        {
        double centerX = 00.;
        double centerY = 0.;
        double auxcenter = (4.-auxradius)/2.;
		
        Point* p100 = fluid2 -> addPoint({ centerX, centerY+0.000},a2,false);
        Point* p101 = fluid2 -> addPoint({ centerX+1.850, centerY+0.000},a2,false);
        Point* p102 = fluid2 -> addPoint({ centerX-1.850, centerY+0.000},a2,false);
        Point* p103 = fluid2 -> addPoint({ centerX-3.850, centerY+0.000},a2,false);
        Point* p104 = fluid2 -> addPoint({ centerX+0.150, centerY+0.000},a2,false);
        Point* p105 = fluid2 -> addPoint({ centerX-0.150, centerY+0.000},a2,false);
        Point* p106 = fluid2 -> addPoint({ centerX+3.850, centerY+0.000},a2,false);
        Point* p107 = fluid2 -> addPoint({ centerX+auxcenter-0.15,centerY},a2,false);
        Point* p108 = fluid2 -> addPoint({ centerX-auxcenter+0.15,centerY},a2,false);
        Point* p109 = fluid2 -> addPoint({ centerX-3.850+auxradius,centerY},a2,false);
        Point* p110 = fluid2 -> addPoint({ centerX+3.850-auxradius,centerY},a2,false);
        
        Line* l100 = fluid2 -> addCircle({p104,p102,p103}); // l0
        Line* l101 = fluid2 -> addCircle({p105,p101,p106}); // l1
        Line* l106 = fluid2 -> addLine({p106,p110}); // l2
        Line* l107 = fluid2 -> addLine({p103,p109}); // l3
        Line* l110 = fluid2 -> addCircle({p104,p108,p109}); // l4
        Line* l111 = fluid2 -> addCircle({p105,p107,p110}); // l5
        Line* l115 = fluid2 -> addCircle({p104,p100,p105}); // l6
        Line* l116 = fluid2 -> addCircle({p105,p100,p104}); // l7
        //Glue zone
        double rx = 5.5;
        double ry = 4.0;
        double angle = 10.;

        Line * l125 = fluid2 -> addEllipse({p100}, rx, ry, 0, 360, angle);
        Line * l129 = fluid2 -> addEllipse({p100}, rx+gzsize, ry+gzsize, 0, 360, angle);
        
        fluid2 -> transfiniteLine({l125}, 101, 1);
        fluid2 -> transfiniteLine({l129}, 101, 1);

        LineLoop* ll101 = fluid2->addLineLoop({ l125});
        LineLoop* ll105 = fluid2->addLineLoop({ l129->operator-()});
        LineLoop* ll102 = fluid2->addLineLoop({ l100, l107, l110->operator-()});
        LineLoop* ll103 = fluid2->addLineLoop({ l116, l115});
        LineLoop* ll104 = fluid2->addLineLoop({ l101, l106, l111->operator-()});

        std::vector<LineLoop*> lin = {ll101,ll102,ll103,ll104};
        PlaneSurface* s101 = fluid2 -> addPlaneSurface(lin);
        std::vector<LineLoop*> lin2 = {ll101,ll105};
        PlaneSurface* s102 = fluid2 -> addPlaneSurface(lin2);
        
        
        
        fluid2 -> addBoundaryCondition("NEUMANN", l125, {}, {}, {}, "GLOBAL");
        
        fluid2 -> addBoundaryCondition("GLUE", l129, {0}, {0}, {}, "GLOBAL");
        
        fluid2 -> addBoundaryCondition("FSINTERFACE", l100, {0}, {0}, {}, "GLOBAL"); // FSI
        fluid2 -> addBoundaryCondition("FSINTERFACE", l101, {0}, {0}, {}, "GLOBAL"); // FSI
        fluid2 -> addBoundaryCondition("FSINTERFACE", l106, {0}, {0}, {}, "GLOBAL"); // FSI
        fluid2 -> addBoundaryCondition("FSINTERFACE", l107, {0}, {0}, {}, "GLOBAL"); // FSI
        fluid2 -> addBoundaryCondition("FSINTERFACE", l110, {0}, {0}, {}, "GLOBAL"); // FSI
        fluid2 -> addBoundaryCondition("FSINTERFACE", l111, {0}, {0}, {}, "GLOBAL"); // FSI
        fluid2 -> addBoundaryCondition("DIRICHLET", l115, {0}, {0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l116, {0}, {0}, {}, "GLOBAL");    
        }


        {
        double centerX = 14.;
        double centerY = 0.;
        double auxcenter = (4.-auxradius)/2.;
		
        Point* p100 = fluid2 -> addPoint({ centerX, centerY+0.000},a2,false);
        Point* p101 = fluid2 -> addPoint({ centerX+1.850, centerY+0.000},a2,false);
        Point* p102 = fluid2 -> addPoint({ centerX-1.850, centerY+0.000},a2,false);
        Point* p103 = fluid2 -> addPoint({ centerX-3.850, centerY+0.000},a2,false);
        Point* p104 = fluid2 -> addPoint({ centerX+0.150, centerY+0.000},a2,false);
        Point* p105 = fluid2 -> addPoint({ centerX-0.150, centerY+0.000},a2,false);
        Point* p106 = fluid2 -> addPoint({ centerX+3.850, centerY+0.000},a2,false);
        Point* p107 = fluid2 -> addPoint({ centerX+auxcenter-0.15,centerY},a2,false);
        Point* p108 = fluid2 -> addPoint({ centerX-auxcenter+0.15,centerY},a2,false);
        Point* p109 = fluid2 -> addPoint({ centerX-3.850+auxradius,centerY},a2,false);
        Point* p110 = fluid2 -> addPoint({ centerX+3.850-auxradius,centerY},a2,false);
        
        Line* l100 = fluid2 -> addCircle({p104,p102,p103}); // l0
        Line* l101 = fluid2 -> addCircle({p105,p101,p106}); // l1
        Line* l106 = fluid2 -> addLine({p106,p110}); // l2
        Line* l107 = fluid2 -> addLine({p103,p109}); // l3
        Line* l110 = fluid2 -> addCircle({p104,p108,p109}); // l4
        Line* l111 = fluid2 -> addCircle({p105,p107,p110}); // l5
        Line* l115 = fluid2 -> addCircle({p104,p100,p105}); // l6
        Line* l116 = fluid2 -> addCircle({p105,p100,p104}); // l7
        //Glue zone
        double rx = 5.5;
        double ry = 4.0;
        double angle = 10.;

        Line * l125 = fluid2 -> addEllipse({p100}, rx, ry, 0, 360, angle);
        Line * l129 = fluid2 -> addEllipse({p100}, rx+gzsize, ry+gzsize, 0, 360, angle);
        
        fluid2 -> transfiniteLine({l125}, 101, 1);
        fluid2 -> transfiniteLine({l129}, 101, 1);

        LineLoop* ll101 = fluid2->addLineLoop({ l125});
        LineLoop* ll105 = fluid2->addLineLoop({ l129->operator-()});
        LineLoop* ll102 = fluid2->addLineLoop({ l100, l107, l110->operator-()});
        LineLoop* ll103 = fluid2->addLineLoop({ l116, l115});
        LineLoop* ll104 = fluid2->addLineLoop({ l101, l106, l111->operator-()});

        std::vector<LineLoop*> lin = {ll101,ll102,ll103,ll104};
        PlaneSurface* s101 = fluid2 -> addPlaneSurface(lin);
        std::vector<LineLoop*> lin2 = {ll101,ll105};
        PlaneSurface* s102 = fluid2 -> addPlaneSurface(lin2);
        
        
        
        fluid2 -> addBoundaryCondition("NEUMANN", l125, {}, {}, {}, "GLOBAL");
        
        fluid2 -> addBoundaryCondition("GLUE", l129, {0}, {0}, {}, "GLOBAL");
        
        fluid2 -> addBoundaryCondition("FSINTERFACE", l100, {0}, {0}, {}, "GLOBAL"); // FSI
        fluid2 -> addBoundaryCondition("FSINTERFACE", l101, {0}, {0}, {}, "GLOBAL"); // FSI
        fluid2 -> addBoundaryCondition("FSINTERFACE", l106, {0}, {0}, {}, "GLOBAL"); // FSI
        fluid2 -> addBoundaryCondition("FSINTERFACE", l107, {0}, {0}, {}, "GLOBAL"); // FSI
        fluid2 -> addBoundaryCondition("FSINTERFACE", l110, {0}, {0}, {}, "GLOBAL"); // FSI
        fluid2 -> addBoundaryCondition("FSINTERFACE", l111, {0}, {0}, {}, "GLOBAL"); // FSI
        fluid2 -> addBoundaryCondition("DIRICHLET", l115, {0}, {0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l116, {0}, {0}, {}, "GLOBAL");    
        }


    if (rank == 0){
  
        FluidDomain* problem = new FluidDomain(fluid1);
        // problem -> addSurfaceMaterial({ s20,s21,s22,s23 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
        problem -> generateMesh(T6, FRONT, "coarse", "", false, true);

        FluidDomain* problem2 = new FluidDomain(fluid2);
        // problem2 -> addSurfaceMaterial({ s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
        problem2 -> generateMesh(T6, FRONT, "fine", "", false, true);

        // FluidDomain* problem = new FluidDomain(fluid1);
        // problem -> addSurfaceMaterial({ s001,s002,s003,s004,s005,s006,s007,s008,s009 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
        // problem -> generateMesh("T6", "DELAUNAY", "coarse", "", false, true);

        // FluidDomain* problem2 = new FluidDomain(fluid2);
        // problem2 -> addSurfaceMaterial({ s101,s102,s103 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
        // problem2 -> generateMesh("T6", "DELAUNAY", "fine", "", false, true);

        //problem -> readInput("exemplo.msh",0);
	};

//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================

	MPI_Barrier(PETSC_COMM_WORLD);   
    
    coarseModel.meshReading(fluid1,"../turbine_data.txt","coarse.msh","mirror.txt",0);
    fineModel.meshReading(fluid2,"../turbine_data.txt","fine.msh","mirror_fine.txt",0);
   // } 
	MPI_Barrier(PETSC_COMM_WORLD);

    char in_solid[32] = "cantilever.txt";
    coarseModel.getFluidParameters().setSolver(SolverType::EMumps);
    fineModel.getFluidParameters().setSolver(SolverType::EMumps);


    coarseModel.setProblemType(ProblemType::ENavierStokes);
    fineModel.setProblemType(ProblemType::ENavierStokes);
        
  
    arlequinProblem.setFluidModels(coarseModel, fineModel);
    coupledProblem.setArlequinAndSolidModels(arlequinProblem,in_solid);

    arlequinProblem.UnsetUseSNES();
    arlequinProblem.solveArlequinProblem(4, 1.e-7, 2, 1);
    
    // coupledProblem.solveFSIProblemArlequin(100000);

        // coupledProblem.solveFSIProblemGaussSeidelArlequin(1000);

    //control.dataReading(".txt","mirror_coarse.txt");
  
    //Problem solving  
         
    //Solve Problem function needs three parameters:  
    //1- The maximum number of iterations in the Newton-Raphson process
    //2- The maximum relative error in the Newton-Raphson process (DU)
    //3- The type of problem to be solved:  
    //   1- Steady Stokes problem;   
    //   2- Steady Navier-Stokes problem (solves the steady Stokes
    //      in the first step to perform the initial guess);  
     
    //coarseModel.solveTransientProblem(2, 1.e-7, 2);
 
    //fluidCoarseMesh.solveTransientProblem(4, 1.e-16, 2); 
  
    
    //control.solveTransientProblem(2, 1.e-6, 2); 
     
    //fineModel.solveTransientProblem(2, 1.e-16, 2);  
 
    // fineModel.solveSteadyLaplaceProblem(8, 1.e-16, 2);
         
  
    //1- The maximum number of iterations in the Newton-Raphson process
    //2- The maximum relative error in the Newton-Raphson process (DU)
    //3- The type of problem to be solved:
    //   1- Stokes problem;   
    //   2- Navier-Stokes problem (solves the steady Stokes
    //      in the first step to perform the initial guess);
    //4- 0 - Steady problem 
    //   1 - Transient problem  
      
    
	// coarseModel.readInitialValues("initialVelocityCoarse7500.txt","initialPressureCoarse7500.txt");
	// fineModel.readInitialValues("initialVelocityFine7500.txt","initialPressureFine7500.txt");

    // arlequinProblem.setFluidModels(coarseModel, fineModel) ; 

    // arlequinProblem.solveArlequinProblem(3, 1.e-7, 2, 1); 
           