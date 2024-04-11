
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
   //FSI        coupledProblem;  
    
//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================
        //Coarse
        Geometry* fluid1 = new Geometry(0);

        Point* p0 = fluid1 -> addPoint({ 0.0, 0.0 },0.5,false);
        Point* p1 = fluid1 -> addPoint({ 0.5, 0.0 },0.5,false);
        Point* p2 = fluid1 -> addPoint({ 1.0, 0.0 },0.5,false);
        Point* p3 = fluid1 -> addPoint({ 0.0, 0.5 },0.5,false);
        Point* p4 = fluid1 -> addPoint({ 0.5, 0.5 },0.5,false);
        Point* p5 = fluid1 -> addPoint({ 1.0, 0.5 },0.5,false);
        Point* p6 = fluid1 -> addPoint({ 0.0, 1.0 },0.5,false);
        Point* p7 = fluid1 -> addPoint({ 0.5, 1.0 },0.5,false);
        Point* p8 = fluid1 -> addPoint({ 1.0, 1.0 },0.5,false);

        Line* l0 = fluid1 -> addLine({ p0, p1 });
        Line* l1 = fluid1 -> addLine({ p1, p2 });
        Line* l2 = fluid1 -> addLine({ p0, p3 });
        Line* l3 = fluid1 -> addLine({ p1, p4 }); 
        Line* l4 = fluid1 -> addLine({ p2, p5 });
        Line* l5 = fluid1 -> addLine({ p3, p4 });
        Line* l6 = fluid1 -> addLine({ p4, p5 });
        Line* l7 = fluid1 -> addLine({ p3, p6 }); 
        Line* l8 = fluid1 -> addLine({ p4, p7 });
        Line* l9 = fluid1 -> addLine({ p5, p8 });
        Line* l110 = fluid1 -> addLine({ p6, p7 });
        Line* l111 = fluid1 -> addLine({ p7, p8 });
 
        LineLoop* ll0 = fluid1->addLineLoop({ l0, l3, l5 -> operator-(), l2 -> operator-()});
        LineLoop* ll1 = fluid1->addLineLoop({ l1, l4, l6 -> operator-(), l3 -> operator-()});
        LineLoop* ll2 = fluid1->addLineLoop({ l5, l8, l110 -> operator-(), l7 -> operator-()});
        LineLoop* ll3 = fluid1->addLineLoop({ l6, l9, l111 -> operator-(), l8 -> operator-()});

        //std::vector<LineLoop*> lin = {ll0, ll1};
        PlaneSurface* s20 = fluid1 -> addPlaneSurface({ll0});
        PlaneSurface* s21 = fluid1 -> addPlaneSurface({ll1});
        PlaneSurface* s22 = fluid1 -> addPlaneSurface({ll2});
        PlaneSurface* s23 = fluid1 -> addPlaneSurface({ll3});
        double h1 = 11; double v1 = 11;
        fluid1 -> transfiniteLine({ l0 }, h1);
        fluid1 -> transfiniteLine({ l1 }, v1);
        fluid1 -> transfiniteLine({ l2 }, h1);
        fluid1 -> transfiniteLine({ l3 }, v1);
        fluid1 -> transfiniteLine({ l4 }, v1);
        fluid1 -> transfiniteLine({ l5 }, v1);
        fluid1 -> transfiniteLine({ l6 }, v1);
        fluid1 -> transfiniteLine({ l7 }, v1);
        fluid1 -> transfiniteLine({ l8 }, v1);
        fluid1 -> transfiniteLine({ l9 }, v1);
        fluid1 -> transfiniteLine({ l110 }, v1);
        fluid1 -> transfiniteLine({ l111 }, v1);
        fluid1->transfiniteSurface({ s20 }, "Left", {p0,p1,p4,p3});
        fluid1->transfiniteSurface({ s21 }, "Right", {p1,p2,p5,p4});
        fluid1->transfiniteSurface({ s22 }, "Right", {p3,p4,p7,p6});
        fluid1->transfiniteSurface({ s23 }, "Left", {p4,p5,p8,p7});

        fluid1 -> addBoundaryCondition("DIRICHLET", l110, {1.0}, {0.0}, {},  "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l111, {1.0}, {0.0}, {},  "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l0, {0.0}, {0.0}, {},  "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l1, {0.0}, {0.0}, {},  "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l2, {0.0}, {0.0}, {},  "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l3, {}, {}, {},  "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l4, {0.0}, {0.0}, {},  "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l5, {}, {}, {},  "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l6, {}, {}, {},  "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l7, {0.0}, {0.0}, {},  "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l8, {}, {}, {},  "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l9, {0.0}, {0.0}, {},  "GLOBAL");



        // Fine
        Geometry* fluid2 = new Geometry(1);

        Point* p10 = fluid2 -> addPoint({0.0,0.0},0.1,false);
        Point* p11 = fluid2 -> addPoint({0.2,0.0},0.1,false);
        Point* p12 = fluid2 -> addPoint({0.5,0.0},0.1,false);
        Point* p13 = fluid2 -> addPoint({0.8,0.0},0.1,false);
        Point* p14 = fluid2 -> addPoint({1.0,0.0},0.1,false);
        Point* p15 = fluid2 -> addPoint({0.0,0.2},0.1,false);
        Point* p16 = fluid2 -> addPoint({0.2,0.2},0.1,false);
        Point* p17 = fluid2 -> addPoint({0.5,0.2},0.1,false);
        Point* p18 = fluid2 -> addPoint({0.8,0.2},0.1,false);
        Point* p19 = fluid2 -> addPoint({1.0,0.2},0.1,false);
        Point* p20 = fluid2 -> addPoint({0.0,0.5},0.1,false);
        Point* p21 = fluid2 -> addPoint({0.2,0.5},0.1,false);
        Point* p22 = fluid2 -> addPoint({0.8,0.5},0.1,false);
        Point* p23 = fluid2 -> addPoint({1.0,0.5},0.1,false);
        Point* p24 = fluid2 -> addPoint({0.0,0.8},0.1,false);
        Point* p25 = fluid2 -> addPoint({0.2,0.8},0.1,false);
        Point* p26 = fluid2 -> addPoint({0.5,0.8},0.1,false);
        Point* p27 = fluid2 -> addPoint({0.8,0.8},0.1,false);
        Point* p28 = fluid2 -> addPoint({1.0,0.8},0.1,false);
        Point* p29 = fluid2 -> addPoint({0.0,1.0},0.1,false);
        Point* p30 = fluid2 -> addPoint({0.2,1.0},0.1,false);
        Point* p31 = fluid2 -> addPoint({0.5,1.0},0.1,false);
        Point* p32 = fluid2 -> addPoint({0.8,1.0},0.1,false);
        Point* p33 = fluid2 -> addPoint({1.0,1.0},0.1,false);

        Line* l10 = fluid2 -> addLine({p10,p11});
        Line* l11 = fluid2 -> addLine({p11,p12});
        Line* l12 = fluid2 -> addLine({p12,p13});
        Line* l13 = fluid2 -> addLine({p13,p14});
        Line* l14 = fluid2 -> addLine({p10,p15});
        Line* l15 = fluid2 -> addLine({p11,p16});
        Line* l16 = fluid2 -> addLine({p12,p17});
        Line* l17 = fluid2 -> addLine({p13,p18});
        Line* l18 = fluid2 -> addLine({p14,p19});
        Line* l19 = fluid2 -> addLine({p15,p16});
        Line* l20 = fluid2 -> addLine({p16,p17});
        Line* l21 = fluid2 -> addLine({p17,p18});
        Line* l22 = fluid2 -> addLine({p18,p19});
        Line* l23 = fluid2 -> addLine({p15,p20});
        Line* l24 = fluid2 -> addLine({p16,p21});
        Line* l25 = fluid2 -> addLine({p18,p22});
        Line* l26 = fluid2 -> addLine({p19,p23});
        Line* l27 = fluid2 -> addLine({p20,p21});
        Line* l28 = fluid2 -> addLine({p22,p23});
        Line* l29 = fluid2 -> addLine({p20,p24});
        Line* l30 = fluid2 -> addLine({p21,p25});
        Line* l31 = fluid2 -> addLine({p22,p27});
        Line* l32 = fluid2 -> addLine({p23,p28});
        Line* l33 = fluid2 -> addLine({p24,p25});
        Line* l34 = fluid2 -> addLine({p25,p26});
        Line* l35 = fluid2 -> addLine({p26,p27});
        Line* l36 = fluid2 -> addLine({p27,p28});
        Line* l37 = fluid2 -> addLine({p24,p29});
        Line* l38 = fluid2 -> addLine({p25,p30});
        Line* l39 = fluid2 -> addLine({p26,p31});
        Line* l40 = fluid2 -> addLine({p27,p32});
        Line* l41 = fluid2 -> addLine({p28,p33});
        Line* l42 = fluid2 -> addLine({p29,p30});
        Line* l43 = fluid2 -> addLine({p30,p31});
        Line* l44 = fluid2 -> addLine({p31,p32});
        Line* l45 = fluid2 -> addLine({p32,p33});

        LineLoop* ll10 = fluid2->addLineLoop({ l10, l15, l19 -> operator-(), l14 -> operator-() });
        LineLoop* ll11 = fluid2->addLineLoop({ l11, l16, l20 -> operator-(), l15 -> operator-() });
        LineLoop* ll12 = fluid2->addLineLoop({ l12, l17, l21 -> operator-(), l16 -> operator-() });
        LineLoop* ll13 = fluid2->addLineLoop({ l13, l18, l22 -> operator-(), l17 -> operator-() });
        LineLoop* ll14 = fluid2->addLineLoop({ l19, l24, l27 -> operator-(), l23 -> operator-() });
        LineLoop* ll15 = fluid2->addLineLoop({ l22, l26, l28 -> operator-(), l25 -> operator-() });
        LineLoop* ll16 = fluid2->addLineLoop({ l27, l30, l33 -> operator-(), l29 -> operator-() });
        LineLoop* ll17 = fluid2->addLineLoop({ l28, l32, l36 -> operator-(), l31 -> operator-() });
        LineLoop* ll18 = fluid2->addLineLoop({ l33, l38, l42 -> operator-(), l37 -> operator-() });
        LineLoop* ll19 = fluid2->addLineLoop({ l34, l39, l43 -> operator-(), l38 -> operator-() });
        LineLoop* ll20 = fluid2->addLineLoop({ l35, l40, l44 -> operator-(), l39 -> operator-() });
        LineLoop* ll21 = fluid2->addLineLoop({ l36, l41, l45 -> operator-(), l40 -> operator-() });


        PlaneSurface* s1 = fluid2 -> addPlaneSurface({ll10});
        PlaneSurface* s2 = fluid2 -> addPlaneSurface({ll11});
        PlaneSurface* s3 = fluid2 -> addPlaneSurface({ll12});
        PlaneSurface* s4 = fluid2 -> addPlaneSurface({ll13});
        PlaneSurface* s5 = fluid2 -> addPlaneSurface({ll14});
        PlaneSurface* s6 = fluid2 -> addPlaneSurface({ll15});
        PlaneSurface* s7 = fluid2 -> addPlaneSurface({ll16});
        PlaneSurface* s8 = fluid2 -> addPlaneSurface({ll17});
        PlaneSurface* s9 = fluid2 -> addPlaneSurface({ll18});
        PlaneSurface* s10 = fluid2 -> addPlaneSurface({ll19});
        PlaneSurface* s11 = fluid2 -> addPlaneSurface({ll20});
        PlaneSurface* s12 = fluid2 -> addPlaneSurface({ll21});

        fluid2 -> addBoundaryCondition("DIRICHLET", l42, {1.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l43, {1.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l44, {1.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l45, {1.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l10, {0.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l11, {0.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l12, {0.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l13, {0.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l14, {0.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l18, {0.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l23, {0.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l26, {0.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l29, {0.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l32, {0.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l37, {0.0}, {0.0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l41, {0.0}, {0.0}, {},  "GLOBAL");
        

        fluid2 -> addBoundaryCondition("NEUMANN", l15, {}, {}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l16, {}, {}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l17, {}, {}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l19, {}, {}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l22, {}, {}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l27, {}, {}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l28, {}, {}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l33, {}, {}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l36, {}, {}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l38, {}, {}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l39, {}, {}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l40, {}, {}, {},  "GLOBAL");
        
        fluid2 -> addBoundaryCondition("GLUE", l20, {0}, {0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l21, {0}, {0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l24, {0}, {0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l25, {0}, {0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l30, {0}, {0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l31, {0}, {0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l34, {0}, {0}, {},  "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l35, {0}, {0}, {},  "GLOBAL");

        //Transfinite lines 
        int corn = 5; int side = 7;
        //corners
        fluid2 -> transfiniteLine({ l10 }, corn);
        fluid2 -> transfiniteLine({ l13 }, corn);
        fluid2 -> transfiniteLine({ l14 }, corn);
        fluid2 -> transfiniteLine({ l15 }, corn);
        fluid2 -> transfiniteLine({ l16 }, corn);
        fluid2 -> transfiniteLine({ l17 }, corn);
        fluid2 -> transfiniteLine({ l18 }, corn);
        fluid2 -> transfiniteLine({ l19 }, corn);
        fluid2 -> transfiniteLine({ l22 }, corn);
        fluid2 -> transfiniteLine({ l27 }, corn);
        fluid2 -> transfiniteLine({ l28 }, corn);
        fluid2 -> transfiniteLine({ l33 }, corn);
        fluid2 -> transfiniteLine({ l36 }, corn);
        fluid2 -> transfiniteLine({ l37 }, corn);
        fluid2 -> transfiniteLine({ l38 }, corn);
        fluid2 -> transfiniteLine({ l39 }, corn);
        fluid2 -> transfiniteLine({ l40 }, corn);
        fluid2 -> transfiniteLine({ l41 }, corn);
        fluid2 -> transfiniteLine({ l42 }, corn);
        fluid2 -> transfiniteLine({ l45 }, corn);

        //sides
        fluid2 -> transfiniteLine({ l11 }, side);
        fluid2 -> transfiniteLine({ l12 }, side);
        fluid2 -> transfiniteLine({ l20 }, side);
        fluid2 -> transfiniteLine({ l21 }, side);
        fluid2 -> transfiniteLine({ l23 }, side);
        fluid2 -> transfiniteLine({ l24 }, side);
        fluid2 -> transfiniteLine({ l25 }, side);
        fluid2 -> transfiniteLine({ l26 }, side);
        fluid2 -> transfiniteLine({ l29 }, side);
        fluid2 -> transfiniteLine({ l30 }, side);
        fluid2 -> transfiniteLine({ l31 }, side);
        fluid2 -> transfiniteLine({ l32 }, side);
        fluid2 -> transfiniteLine({ l34 }, side);
        fluid2 -> transfiniteLine({ l35 }, side);
        fluid2 -> transfiniteLine({ l43 }, side);
        fluid2 -> transfiniteLine({ l44 }, side);

        //Surfaces
        fluid2->transfiniteSurface({ s1 }, "Left", {p10,p11,p16,p15});
        fluid2->transfiniteSurface({ s2 }, "Left", {p11,p12,p17,p16});
        fluid2->transfiniteSurface({ s3 }, "Right", {p12,p13,p18,p17});
        fluid2->transfiniteSurface({ s4 }, "Right", {p13,p14,p19,p18});
        fluid2->transfiniteSurface({ s5 }, "Left", {p15,p16,p21,p20});
        fluid2->transfiniteSurface({ s6 }, "Right", {p18,p19,p23,p22});
        fluid2->transfiniteSurface({ s7 }, "Right", {p20,p21,p25,p24});
        fluid2->transfiniteSurface({ s8 }, "Left", {p22,p23,p28,p27});
        fluid2->transfiniteSurface({ s9 }, "Right", {p24,p25,p30,p29});
        fluid2->transfiniteSurface({ s10 }, "Right", {p25,p26,p31,p30});
        fluid2->transfiniteSurface({ s11 }, "Left", {p26,p27,p32,p31});
        fluid2->transfiniteSurface({ s12 }, "Left", {p27,p28,p33,p32});


    if (rank == 0){
  
        FluidDomain* problem = new FluidDomain(fluid1);
        // problem -> addSurfaceMaterial({ s20,s21,s22,s23 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
        problem -> generateMesh(T6, DELAUNAY, "coarse", "", false, true);

        FluidDomain* problem2 = new FluidDomain(fluid2);
        // problem2 -> addSurfaceMaterial({ s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
        problem2 -> generateMesh(T6, DELAUNAY, "fine", "", false, true);

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


 	//control.dataReading("cylinder_coarse_ref.txt","mirror.txt");
    //control.solveTransientProblem(2, 1.e-6, 2); 

    //Data reading     
    // Cylinder
    // fineModel.dataReading("cyl_fine_str3.txt","mirror_fine.txt");  
    // coarseModel.dataReading("cyl_coarse2.txt","mirror_coarse.txt");
    // fineModel.dataReading("cylinder_fine_ref2.txt","mirror_fine.txt");  
    // coarseModel.dataReading("cylinder_coarse_ref.txt","mirror_coarse.txt");
    // fineModel.dataReading("cyl_fine_novo3.txt","mirror_fine.txt");  
    // coarseModel.dataReading("cyl_coarse_novo2.txt","mirror_coarse.txt");
  
    // Cavity   
    // fineModel.dataReading("fine4.txt","mirror_fine.txt"); 
    // coarseModel.dataReading("coarse.txt","mirror_coarse.txt");
  
    // Helice 
    // fineModel.dataReading("cavpropfcir.txt","mirror_fine.txt"); 
    // coarseModel.dataReading("cavity1.txt","mirror_coarse.txt");

    // Flutter
    // fineModel.dataReading("flutter_fine2.txt","mirror_fine.txt"); 
    // coarseModel.dataReading("flutter_coarse.txt","mirror_coarse.txt");
     
    //control.dataReading("32x32.txt","mirror_control.txt");
 
    //char in_solid[32] = "cantilever.txt";
        
  
    // arlequinProblem.setFluidModels(coarseModel, fineModel);
    // coupledProblem.setArlequinAndSolidModels(arlequinProblem,in_solid);
 
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

    arlequinProblem.setFluidModels(coarseModel, fineModel) ; 

    arlequinProblem.solveArlequinProblem(3, 1.e-7, 2, 1); 
           