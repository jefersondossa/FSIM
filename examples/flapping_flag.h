
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
    FSI        coupledProblem;  
    
//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================
        //Coarse
        Geometry* fluid1 = new Geometry(0);

        double s = 0.5;

        Point* p000 = fluid1 -> addPoint({0.00,0.00},s,false);
        Point* p001 = fluid1 -> addPoint({21.0,0.00},s,false);
        Point* p002 = fluid1 -> addPoint({21.0,12.0},s,false);
        Point* p003 = fluid1 -> addPoint({0.00,12.0},s,false);

        Line* l000 = fluid1 -> addLine({p000,p001});
        Line* l001 = fluid1 -> addLine({p001,p002});
        Line* l002 = fluid1 -> addLine({p002,p003});
        Line* l003 = fluid1 -> addLine({p003,p000});
        
        LineLoop* ll001 = fluid1 -> addLineLoop({ l000, l001, l002, l003 });

        PlaneSurface* s001 = fluid1 -> addPlaneSurface({ll001});
  
        fluid1 -> addBoundaryCondition("NEUMANN", l001, {}, {}, {}, "GLOBAL");

        fluid1 -> addBoundaryCondition("DIRICHLET", l000, {}, {0.0}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l002, {}, {0.0}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l003, {31.5}, {0.0}, {}, "GLOBAL");
        

        // Fine
        Geometry* fluid2 = new Geometry(1);
        double s2 = 0.1;
        double s3 = 0.3;

        Point* p100 = fluid2 -> addPoint({5.5,5.5},s2,false);
        Point* p101 = fluid2 -> addPoint({6.5,5.5},s2,false);
        Point* p102 = fluid2 -> addPoint({5.5,6.5},s2,false);
        Point* p103 = fluid2 -> addPoint({6.5,6.5},s2,false);
        Point* p104 = fluid2 -> addPoint({6.5,6.0001},s2,false);
        Point* p105 = fluid2 -> addPoint({6.5,5.9999},s2,false);
        Point* p106 = fluid2 -> addPoint({10.5,6.0},s2,false);
        Point* p107 = fluid2 -> addPoint({4.0,4.0},s3,false);
        Point* p108 = fluid2 -> addPoint({12.0,4.0},s3,false);
        Point* p109 = fluid2 -> addPoint({4.0,8.0},s3,false);
        Point* p110 = fluid2 -> addPoint({12.0,8.0},s3,false);
        Point* p111 = fluid2 -> addPoint({4.0,3.0},s3,false);
        Point* p112 = fluid2 -> addPoint({12.0,3.0},s3,false);
        Point* p113 = fluid2 -> addPoint({13.0,4.0},s3,false);
        Point* p114 = fluid2 -> addPoint({13.0,8.0},s3,false);
        Point* p115 = fluid2 -> addPoint({12.0,9.0},s3,false);
        Point* p116 = fluid2 -> addPoint({4.0,9.0},s3,false);
        Point* p117 = fluid2 -> addPoint({3.0,8.0},s3,false);
        Point* p118 = fluid2 -> addPoint({3.0,4.0},s3,false);
        
        Line* l100 = fluid2 -> addLine({p100,p101});
        Line* l101 = fluid2 -> addLine({p100,p102});
        Line* l102 = fluid2 -> addLine({p102,p103});
        Line* l103 = fluid2 -> addLine({p104,p103});
        Line* l104 = fluid2 -> addLine({p101,p105});
        Line* l105 = fluid2 -> addLine({p105,p106});
        Line* l106 = fluid2 -> addLine({p104,p106});
        Line* l107 = fluid2 -> addLine({p107,p108});
        Line* l108 = fluid2 -> addLine({p108,p110});
        Line* l109 = fluid2 -> addLine({p109,p110});
        Line* l110 = fluid2 -> addLine({p107,p109});
        Line* l111 = fluid2 -> addLine({p118,p107});
        Line* l112 = fluid2 -> addLine({p111,p107});
        Line* l113 = fluid2 -> addLine({p112,p108});
        Line* l114 = fluid2 -> addLine({p108,p113});
        Line* l115 = fluid2 -> addLine({p110,p114});
        Line* l116 = fluid2 -> addLine({p110,p115});
        Line* l117 = fluid2 -> addLine({p109,p116});
        Line* l118 = fluid2 -> addLine({p117,p109});
        Line* l119 = fluid2 -> addLine({p111,p112});
        Line* l120 = fluid2 -> addCircle({p112,p108,p113});
        Line* l121 = fluid2 -> addLine({p113,p114});
        Line* l122 = fluid2 -> addCircle({p114,p110,p115});
        Line* l123 = fluid2 -> addLine({p116,p115});
        Line* l124 = fluid2 -> addCircle({p117,p109,p116});
        Line* l125 = fluid2 -> addLine({p118,p117});
        Line* l126 = fluid2 -> addCircle({p111,p107,p118});

        LineLoop* ll100 = fluid2->addLineLoop({ l107, l108, l109 -> operator-(), l110 -> operator-(), l105 -> operator-(), l104 -> operator-(), l100 -> operator-(), l101, l102, l103 -> operator-(), l106 });
        LineLoop* ll101 = fluid2->addLineLoop({ l118, l117, l124 -> operator-() });
        LineLoop* ll102 = fluid2->addLineLoop({ l109, l116, l123 -> operator-(), l117 -> operator-() });
        LineLoop* ll103 = fluid2->addLineLoop({ l115, l122, l116 -> operator-() });
        LineLoop* ll104 = fluid2->addLineLoop({ l111, l110, l118 -> operator-(), l125 -> operator-() });
        LineLoop* ll105 = fluid2->addLineLoop({ l114, l121, l115 -> operator-(), l108 -> operator-() });
        LineLoop* ll106 = fluid2->addLineLoop({ l112, l111 -> operator-(), l126 -> operator-() });
        LineLoop* ll107 = fluid2->addLineLoop({ l119, l113, l107 -> operator-(), l112 -> operator-() });
        LineLoop* ll108 = fluid2->addLineLoop({ l120, l114 -> operator-(), l113 -> operator-() });


        //std::vector<LineLoop*> lin = {ll2, ll1};
        //std::vector<LineLoop*> lin = {ll2};
        //PlaneSurface* s101 = fluid2 -> addPlaneSurface(lin);
        PlaneSurface* s100 = fluid2 -> addPlaneSurface({ll100});
        PlaneSurface* s101 = fluid2 -> addPlaneSurface({ll101});
        PlaneSurface* s102 = fluid2 -> addPlaneSurface({ll102});
        PlaneSurface* s103 = fluid2 -> addPlaneSurface({ll103});
        PlaneSurface* s104 = fluid2 -> addPlaneSurface({ll104});
        PlaneSurface* s105 = fluid2 -> addPlaneSurface({ll105});
        PlaneSurface* s106 = fluid2 -> addPlaneSurface({ll106});
        PlaneSurface* s107 = fluid2 -> addPlaneSurface({ll107});
        PlaneSurface* s108 = fluid2 -> addPlaneSurface({ll108});       

        fluid2 -> addBoundaryCondition("NEUMANN", l111, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l112, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l113, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l114, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l115, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l116, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l117, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l118, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l107, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l108, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l109, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l110, {}, {}, {}, "GLOBAL");
        
        fluid2 -> addBoundaryCondition("GLUE", l119, {0}, {0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l120, {0}, {0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l121, {0}, {0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l122, {0}, {0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l123, {0}, {0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l124, {0}, {0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l125, {0}, {0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l126, {0}, {0}, {}, "GLOBAL");
        
        fluid2 -> addBoundaryCondition("FSINTERFACE", l105, {0}, {0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("FSINTERFACE", l106, {0}, {0}, {}, "GLOBA {},L");
         
        fluid2 -> addBoundaryCondition("DIRICHLET", l100, {0}, {0}, {}, "GLOBAL"); 
        fluid2 -> addBoundaryCondition("DIRICHLET", l101, {0}, {0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l102, {0}, {0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l103, {0}, {0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l104, {0}, {0}, {}, "GLOBAL");
        
        
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
    
    
    coarseModel.meshReading(fluid1,"flapping_flag_data.txt","coarse.msh","mirror.txt",0);
    fineModel.meshReading(fluid2,"flapping_flag_data.txt","fine.msh","mirror_fine.txt",0);
   // } 
	MPI_Barrier(PETSC_COMM_WORLD);

    coarseModel.getFluidParameters().setSolver(SolverType::EMumps);
    fineModel.getFluidParameters().setSolver(SolverType::EMumps);


    coarseModel.setProblemType(ProblemType::ENavierStokes);
    fineModel.setProblemType(ProblemType::ENavierStokes);

    arlequinProblem.setArlequinStabilization(ArlequinStabType::EOption2);
    arlequinProblem.setFluidModels(coarseModel, fineModel) ; 


    char in_solid[32] = "cantilever.txt";

    coupledProblem.setArlequinAndSolidModels(arlequinProblem,in_solid);
    coupledProblem.solveFSIProblemArlequin(100000);

    // arlequinProblem.solveArlequinProblem(3, 1.e-7, 2, 1); 
           