//------------------------------------------------------------------------------
//--------------------------Universidade de Sao Paulo---------------------------
//----------------------Escola de Engenharia de Sao Carlos----------------------
//----------------Departamento de Engenharia de Estruturas - SET----------------
//------------------------------Sao Carlos - 2017-------------------------------
//------------------------------------------------------------------------------
 
///-----------------------------------------------------------------------------
///----Software developed for analysis of Fluid-Structure Interaction (FSI)-----
///---problems with shell structures and incompressible fluid. The fluid is-----
///---modeled in Arlequin framework and Arbitrary Lagrangian-Eulerian (ALE)-----
///---description in mixed formulation, with triangular/tetrahedral elements----
///with quadratic (six/ten nodes) and linear (four/three nodes) approximations--
///--for velocity and pressure fields, respectively. The structure is modeled---
///-with shell elements free rotations as degrees of freedom in the nonlinear---
///-Positional Finite Element Method framework. Triangular elements with cubic--
///approximation (ten nodes) and seven degrees of freedom per node (three nodal-
///--positions, three components of unconstrained vector and strain rate. The---
///-fluid-shell coupling is performed by a partitioned loosely coupled scheme.--
///-----------------------------------------------------------------------------
   
//------------------------------------------------------------------------------
//---------------------------------Developed by---------------------------------
//-------Jeferson Wilian Dossa Fernandes and Rodolfo Andre Kuche Sanches--------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------STARTS MAIN PROGRAM-----------------------------
//------------------------------------------------------------------------------
static char help[] = "Solves the Incompressible flow problem";

// C++ standard libraries
#include <fstream> 
  
// Developed Header Files
#include "src/FSInteraction.hpp"  
#include "src/fluidDomain.h"


int main(int argc, char **args) {

    // Starts main program invoking PETSc
    PetscInitialize(&argc, &args, (char*)0, help);

    int rank, size;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
 
    if (rank == 0){
        std::cout << "2D Incompressible Flows Numerical Analysis" << std::endl;
        std::cout << "Starting.." << std::endl;
        std::cout << "Type the input file name:" << std::endl;
    };

    // Defines the problem dimension
    const int dimension = 2;

    //Type definition
    typedef Fluid<dimension>         FluidModel;
    typedef Arlequin<dimension>      Arlequin;
    typedef FSInteraction<dimension> FSI;
 
//  Create problem variables 
    FluidModel coarseModel, fineModel, control;  
    Arlequin   arlequinProblem; 
   //FSI        coupledProblem;  
    
//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================
        //Coarse
    Geometry* fluid1 = new Geometry(0);

    double x1 = -9.5;
    double x2 = -0.6;
    double x3 = 1.6;
    double x4 = 20.5;

	double y1 = -10;
    double y2 = -0.85;
    double y3 = 0.85;
    double y4 = 10;
 
    double h1 = 13; double h2 = 25; double h3 = 25;
    double v1 = 13; double v2 = 15;

    Point* p001 = fluid1 -> addPoint({x1,y1},0.5,false);
    Point* p002 = fluid1 -> addPoint({x2,y1},0.5,false);
    Point* p003 = fluid1 -> addPoint({x3,y1},0.5,false); 
    Point* p004 = fluid1 -> addPoint({x4,y1},0.5,false);
    Point* p005 = fluid1 -> addPoint({x1,y2},0.5,false);
    Point* p006 = fluid1 -> addPoint({x2,y2},0.5,false);
    Point* p007 = fluid1 -> addPoint({x3,y2},0.5,false);
    Point* p008 = fluid1 -> addPoint({x4,y2},0.5,false);
    Point* p009 = fluid1 -> addPoint({x1,y3},0.5,false);
    Point* p010 = fluid1 -> addPoint({x2,y3},0.5,false);
    Point* p011 = fluid1 -> addPoint({x3,y3},0.5,false);
    Point* p012 = fluid1 -> addPoint({x4,y3},0.5,false);
    Point* p013 = fluid1 -> addPoint({x1,y4},0.5,false);
    Point* p014 = fluid1 -> addPoint({x2,y4},0.5,false);
    Point* p015 = fluid1 -> addPoint({x3,y4},0.5,false);
    Point* p016 = fluid1 -> addPoint({x4,y4},0.5,false);

    Line* l001 = fluid1 -> addLine({p001,p002});
    Line* l002 = fluid1 -> addLine({p002,p003});
    Line* l003 = fluid1 -> addLine({p003,p004});
    Line* l004 = fluid1 -> addLine({p001,p005}); 
    Line* l005 = fluid1 -> addLine({p002,p006});
    Line* l006 = fluid1 -> addLine({p003,p007});
    Line* l007 = fluid1 -> addLine({p004,p008});
    Line* l008 = fluid1 -> addLine({p005,p006}); 
    Line* l009 = fluid1 -> addLine({p006,p007});
    Line* l010 = fluid1 -> addLine({p007,p008});
    Line* l011 = fluid1 -> addLine({p005,p009});
    Line* l012 = fluid1 -> addLine({p006,p010}); 
    Line* l013 = fluid1 -> addLine({p007,p011});
    Line* l014 = fluid1 -> addLine({p008,p012});
    Line* l015 = fluid1 -> addLine({p009,p010});
    Line* l016 = fluid1 -> addLine({p010,p011}); 
    Line* l017 = fluid1 -> addLine({p011,p012});
    Line* l018 = fluid1 -> addLine({p009,p013});
    Line* l019 = fluid1 -> addLine({p010,p014});
    Line* l020 = fluid1 -> addLine({p011,p015}); 
    Line* l021 = fluid1 -> addLine({p012,p016});
    Line* l022 = fluid1 -> addLine({p013,p014});
    Line* l023 = fluid1 -> addLine({p014,p015});
    Line* l024 = fluid1 -> addLine({p015,p016}); 
     
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

    double prog1 = 1.07;
    double prog2 = 1.1;

    fluid1 -> transfiniteLine({l001}, h1, 1/prog2);
    fluid1 -> transfiniteLine({l002}, 4);
    fluid1 -> transfiniteLine({l003}, h3, prog1);
    fluid1 -> transfiniteLine({l004}, v1, 1/prog1);
    fluid1 -> transfiniteLine({l005}, v1, 1/prog1);
    fluid1 -> transfiniteLine({l006}, v1, 1/prog1);
    fluid1 -> transfiniteLine({l007}, v1, 1/prog1);
    fluid1 -> transfiniteLine({l008}, h1, 1/prog2);
    fluid1 -> transfiniteLine({l009}, h2);
    fluid1 -> transfiniteLine({l010}, 50, 1.04);
    fluid1 -> transfiniteLine({l011}, 4);
    fluid1 -> transfiniteLine({l012}, v2);
    fluid1 -> transfiniteLine({l013}, v2);
    fluid1 -> transfiniteLine({l014}, 4);
    fluid1 -> transfiniteLine({l015}, h1, 1/prog2);
    fluid1 -> transfiniteLine({l016}, h2);
    fluid1 -> transfiniteLine({l017}, 50, 1.04);
    fluid1 -> transfiniteLine({l018}, v1, prog1);
    fluid1 -> transfiniteLine({l019}, v1, prog1);
    fluid1 -> transfiniteLine({l020}, v1, prog1);
    fluid1 -> transfiniteLine({l021}, v1, prog1);
    fluid1 -> transfiniteLine({l022}, h1, 1/prog2); 
    fluid1 -> transfiniteLine({l023}, 4);
    fluid1 -> transfiniteLine({l024}, h3, prog1);
  
    // fluid1 -> transfiniteSurface({s001}, "Alternated", {p001,p002,p006,p005});
    // fluid1 -> transfiniteSurface({s002}, "Alternated", {p002,p003,p007,p006});
    // fluid1 -> transfiniteSurface({s003}, "Alternated", {p003,p004,p008,p007});
    // fluid1 -> transfiniteSurface({s004}, "Alternated", {p005,p006,p010,p009});
    // fluid1 -> transfiniteSurface({s005}, "Alternated", {p006,p007,p011,p010});
    // fluid1 -> transfiniteSurface({s006}, "Alternated", {p007,p008,p012,p011});
    // fluid1 -> transfiniteSurface({s007}, "Alternated", {p009,p010,p014,p013});
    // fluid1 -> transfiniteSurface({s008}, "Alternated", {p010,p011,p015,p014});
    // fluid1 -> transfiniteSurface({s009}, "Alternated", {p011,p012,p016,p015});

    fluid1 -> addBoundaryCondition("NEUMANN", l005, {}, {}, "GLOBAL");
    fluid1 -> addBoundaryCondition("NEUMANN", l006, {}, {}, "GLOBAL");
    fluid1 -> addBoundaryCondition("NEUMANN", l008, {}, {}, "GLOBAL");
    fluid1 -> addBoundaryCondition("NEUMANN", l009, {}, {}, "GLOBAL");
    fluid1 -> addBoundaryCondition("NEUMANN", l010, {}, {}, "GLOBAL");
    fluid1 -> addBoundaryCondition("NEUMANN", l012, {}, {}, "GLOBAL");
    fluid1 -> addBoundaryCondition("NEUMANN", l013, {}, {}, "GLOBAL");
    fluid1 -> addBoundaryCondition("NEUMANN", l015, {}, {}, "GLOBAL");
    fluid1 -> addBoundaryCondition("NEUMANN", l016, {}, {}, "GLOBAL");
    fluid1 -> addBoundaryCondition("NEUMANN", l017, {}, {}, "GLOBAL");
    fluid1 -> addBoundaryCondition("NEUMANN", l019, {}, {}, "GLOBAL");
    fluid1 -> addBoundaryCondition("NEUMANN", l020, {}, {}, "GLOBAL");

    fluid1 -> addBoundaryCondition("DIRICHLET", l001, {1.0}, {0.0}, "GLOBAL");
    fluid1 -> addBoundaryCondition("DIRICHLET", l002, {1.0}, {0.0}, "GLOBAL");
    fluid1 -> addBoundaryCondition("DIRICHLET", l003, {1.0}, {0.0}, "GLOBAL");
    fluid1 -> addBoundaryCondition("DIRICHLET", l022, {1.0}, {0.0}, "GLOBAL");
    fluid1 -> addBoundaryCondition("DIRICHLET", l023, {1.0}, {0.0}, "GLOBAL");
    fluid1 -> addBoundaryCondition("DIRICHLET", l024, {1.0}, {0.0}, "GLOBAL");

    fluid1 -> addBoundaryCondition("NEUMANN", l007, {}, {}, "GLOBAL");
    fluid1 -> addBoundaryCondition("NEUMANN", l014, {}, {}, "GLOBAL");
    fluid1 -> addBoundaryCondition("NEUMANN", l021, {}, {}, "GLOBAL");

    fluid1 -> addBoundaryCondition("DIRICHLET", l004, {1.0}, {0.0}, "GLOBAL");
    fluid1 -> addBoundaryCondition("DIRICHLET", l011, {1.0}, {0.0}, "GLOBAL");
    fluid1 -> addBoundaryCondition("DIRICHLET", l018, {1.0}, {0.0}, "GLOBAL");
    


    // Fine
    Geometry* fluid2 = new Geometry(1);
    double elSize = 2.e-0;

    Point* p1001 = fluid2 -> addPoint({1.000000000000000E+00,0.000000000000000E+00},elSize,false);
    Point* p1002 = fluid2 -> addPoint({9.970300197601318E-01,-4.207300080452114E-04},elSize,false);
    Point* p1003 = fluid2 -> addPoint({9.937499761581421E-01,-8.835400221869349E-04},elSize,false);
    Point* p1004 = fluid2 -> addPoint({9.901300072669983E-01,-1.392400008626282E-03},elSize,false);
    Point* p1005 = fluid2 -> addPoint({9.861299991607666E-01,-1.951699960045516E-03},elSize,false);
    Point* p1006 = fluid2 -> addPoint({9.817100167274475E-01,-2.566199982538819E-03},elSize,false);
    Point* p1007 = fluid2 -> addPoint({9.768300056457520E-01,-3.240799997001886E-03},elSize,false);
    Point* p1008 = fluid2 -> addPoint({9.714400172233582E-01,-3.981099929660559E-03},elSize,false);
    Point* p1009 = fluid2 -> addPoint({9.654899835586548E-01,-4.792999941855669E-03},elSize,false);
    Point* p1010 = fluid2 -> addPoint({9.589200019836426E-01,-5.682699847966433E-03},elSize,false);
    Point* p1011 = fluid2 -> addPoint({9.516599774360657E-01,-6.657199934124947E-03},elSize,false);
    Point* p1012 = fluid2 -> addPoint({9.436399936676025E-01,-7.723500020802021E-03},elSize,false);
    Point* p1013 = fluid2 -> addPoint({9.347800016403198E-01,-8.889400400221348E-03},elSize,false);
    Point* p1014 = fluid2 -> addPoint({9.250000119209290E-01,-1.016299985349178E-02},elSize,false);
    Point* p1015 = fluid2 -> addPoint({9.142000079154968E-01,-1.155200041830540E-02},elSize,false);
    Point* p1016 = fluid2 -> addPoint({9.022700190544128E-01,-1.306699961423874E-02},elSize,false);
    Point* p1017 = fluid2 -> addPoint({8.890900015830994E-01,-1.471399981528521E-02},elSize,false);
    Point* p1018 = fluid2 -> addPoint({8.745399713516235E-01,-1.650499925017357E-02},elSize,false);
    Point* p1019 = fluid2 -> addPoint({8.584600090980530E-01,-1.844800077378750E-02},elSize,false);
    Point* p1020 = fluid2 -> addPoint({8.416299819946289E-01,-2.044299989938736E-02},elSize,false);
    Point* p1021 = fluid2 -> addPoint({8.249099850654602E-01,-2.238800004124641E-02},elSize,false);
    Point* p1022 = fluid2 -> addPoint({8.082900047302246E-01,-2.428100071847439E-02},elSize,false);
    Point* p1023 = fluid2 -> addPoint({7.917699813842773E-01,-2.612599916756153E-02},elSize,false);
    Point* p1024 = fluid2 -> addPoint({7.753499746322632E-01,-2.792399935424328E-02},elSize,false);
    Point* p1025 = fluid2 -> addPoint({7.590299844741821E-01,-2.967499941587448E-02},elSize,false);
    Point* p1026 = fluid2 -> addPoint({7.428100109100342E-01,-3.138099983334541E-02},elSize,false);
    Point* p1027 = fluid2 -> addPoint({7.266700267791748E-01,-3.304199874401093E-02},elSize,false);
    Point* p1028 = fluid2 -> addPoint({7.106199860572815E-01,-3.465899825096130E-02},elSize,false);
    Point* p1029 = fluid2 -> addPoint({6.946499943733215E-01,-3.623199835419655E-02},elSize,false);
    Point* p1030 = fluid2 -> addPoint({6.787599921226501E-01,-3.776200115680695E-02},elSize,false);
    Point* p1031 = fluid2 -> addPoint({6.629499793052673E-01,-3.924899920821190E-02},elSize,false);
    Point* p1032 = fluid2 -> addPoint({6.472100019454956E-01,-4.069200158119202E-02},elSize,false);
    Point* p1033 = fluid2 -> addPoint({6.315400004386902E-01,-4.209100082516670E-02},elSize,false);
    Point* p1034 = fluid2 -> addPoint({6.159300208091736E-01,-4.344699904322624E-02},elSize,false);
    Point* p1035 = fluid2 -> addPoint({6.003900170326233E-01,-4.475900158286095E-02},elSize,false);
    Point* p1036 = fluid2 -> addPoint({5.849099755287170E-01,-4.602500051259995E-02},elSize,false);
    Point* p1037 = fluid2 -> addPoint({5.694800019264221E-01,-4.724600166082382E-02},elSize,false);
    Point* p1038 = fluid2 -> addPoint({5.541099905967712E-01,-4.842000082135201E-02},elSize,false);
    Point* p1039 = fluid2 -> addPoint({5.387799739837646E-01,-4.954700171947479E-02},elSize,false);
    Point* p1040 = fluid2 -> addPoint({5.235000252723694E-01,-5.062500014901161E-02},elSize,false);
    Point* p1041 = fluid2 -> addPoint({5.082700252532959E-01,-5.165300145745277E-02},elSize,false);
    Point* p1042 = fluid2 -> addPoint({4.930700063705444E-01,-5.262900143861771E-02},elSize,false);
    Point* p1043 = fluid2 -> addPoint({4.778999984264374E-01,-5.355200171470642E-02},elSize,false);
    Point* p1044 = fluid2 -> addPoint({4.627699851989746E-01,-5.441899970173836E-02},elSize,false);
    Point* p1045 = fluid2 -> addPoint({4.476700127124786E-01,-5.522900074720383E-02},elSize,false);
    Point* p1046 = fluid2 -> addPoint({4.325900077819824E-01,-5.597800016403198E-02},elSize,false);
    Point* p1047 = fluid2 -> addPoint({4.175300002098083E-01,-5.666499957442284E-02},elSize,false);
    Point* p1048 = fluid2 -> addPoint({4.024899899959564E-01,-5.728600174188614E-02},elSize,false);
    Point* p1049 = fluid2 -> addPoint({3.874700069427490E-01,-5.783800035715103E-02},elSize,false);
    Point* p1050 = fluid2 -> addPoint({3.724600076675415E-01,-5.831800028681755E-02},elSize,false);
    Point* p1051 = fluid2 -> addPoint({3.574599921703339E-01,-5.872299894690514E-02},elSize,false);
    Point* p1052 = fluid2 -> addPoint({3.424600064754486E-01,-5.904699862003326E-02},elSize,false);
    Point* p1053 = fluid2 -> addPoint({3.274700045585632E-01,-5.928599834442139E-02},elSize,false);
    Point* p1054 = fluid2 -> addPoint({3.124699890613556E-01,-5.943400040268898E-02},elSize,false);
    Point* p1055 = fluid2 -> addPoint({2.974700033664703E-01,-5.948600172996521E-02},elSize,false);
    Point* p1056 = fluid2 -> addPoint({2.824699878692627E-01,-5.943600088357925E-02},elSize,false);
    Point* p1057 = fluid2 -> addPoint({2.674500048160553E-01,-5.927500128746033E-02},elSize,false);
    Point* p1058 = fluid2 -> addPoint({2.524299919605255E-01,-5.899700149893761E-02},elSize,false);
    Point* p1059 = fluid2 -> addPoint({2.373899966478348E-01,-5.859300121665001E-02},elSize,false);
    Point* p1060 = fluid2 -> addPoint({2.223400026559830E-01,-5.804999917745590E-02},elSize,false);
    Point* p1061 = fluid2 -> addPoint({2.072599977254868E-01,-5.735800042748451E-02},elSize,false);
    Point* p1062 = fluid2 -> addPoint({1.921699941158295E-01,-5.650499835610390E-02},elSize,false);
    Point* p1063 = fluid2 -> addPoint({1.770599931478500E-01,-5.547099933028221E-02},elSize,false);
    Point* p1064 = fluid2 -> addPoint({1.619299948215485E-01,-5.423900112509727E-02},elSize,false);
    Point* p1065 = fluid2 -> addPoint({1.467899978160858E-01,-5.278500169515610E-02},elSize,false);
    Point* p1066 = fluid2 -> addPoint({1.322900056838989E-01,-5.116200074553490E-02},elSize,false);
    Point* p1067 = fluid2 -> addPoint({1.190399974584579E-01,-4.945399984717369E-02},elSize,false);
    Point* p1068 = fluid2 -> addPoint({1.069300025701523E-01,-4.768000170588493E-02},elSize,false);
    Point* p1069 = fluid2 -> addPoint({9.586799889802933E-02,-4.585599899291992E-02},elSize,false);
    Point* p1070 = fluid2 -> addPoint({8.577000349760056E-02,-4.399599879980087E-02},elSize,false);
    Point* p1071 = fluid2 -> addPoint({7.655599713325500E-02,-4.211099818348885E-02},elSize,false);
    Point* p1072 = fluid2 -> addPoint({6.815300136804581E-02,-4.021200165152550E-02},elSize,false);
    Point* p1073 = fluid2 -> addPoint({6.049599871039391E-02,-3.830600157380104E-02},elSize,false);
    Point* p1074 = fluid2 -> addPoint({5.352399870753288E-02,-3.640000149607658E-02},elSize,false);
    Point* p1075 = fluid2 -> addPoint({4.718200117349625E-02,-3.449999913573265E-02},elSize,false);
    Point* p1076 = fluid2 -> addPoint({4.142000153660774E-02,-3.260999917984009E-02},elSize,false);
    Point* p1077 = fluid2 -> addPoint({3.619199991226196E-02,-3.073300048708916E-02},elSize,false);
    Point* p1078 = fluid2 -> addPoint({3.145699948072433E-02,-2.887200005352497E-02},elSize,false);
    Point* p1079 = fluid2 -> addPoint({2.717700041830540E-02,-2.702900022268295E-02},elSize,false);
    Point* p1080 = fluid2 -> addPoint({2.331699989736080E-02,-2.520499937236309E-02},elSize,false);
    Point* p1081 = fluid2 -> addPoint({1.984800025820732E-02,-2.340300008654594E-02},elSize,false);
    Point* p1082 = fluid2 -> addPoint({1.673999987542629E-02,-2.162200026214123E-02},elSize,false);
    Point* p1083 = fluid2 -> addPoint({1.396999973803759E-02,-1.986500062048435E-02},elSize,false);
    Point* p1084 = fluid2 -> addPoint({1.151499990373850E-02,-1.813299953937531E-02},elSize,false);
    Point* p1085 = fluid2 -> addPoint({9.355300106108189E-03,-1.642799936234951E-02},elSize,false);
    Point* p1086 = fluid2 -> addPoint({7.471399847418070E-03,-1.475299987941980E-02},elSize,false);
    Point* p1087 = fluid2 -> addPoint({5.847000051289797E-03,-1.311200018972158E-02},elSize,false);
    Point* p1088 = fluid2 -> addPoint({4.466100130230188E-03,-1.151099987328053E-02},elSize,false);
    Point* p1089 = fluid2 -> addPoint({3.312600078061223E-03,-9.954700246453285E-03},elSize,false);
    Point* p1090 = fluid2 -> addPoint({2.369200112298131E-03,-8.452000096440315E-03},elSize,false);
    Point* p1091 = fluid2 -> addPoint({1.617899979464710E-03,-7.011000066995621E-03},elSize,false);
    Point* p1092 = fluid2 -> addPoint({1.039199996739626E-03,-5.639600101858377E-03},elSize,false);
    Point* p1093 = fluid2 -> addPoint({6.129799876362085E-04,-4.344400018453598E-03},elSize,false);
    Point* p1094 = fluid2 -> addPoint({3.168400144204497E-04,-3.131099976599216E-03},elSize,false);
    Point* p1095 = fluid2 -> addPoint({1.291899970965460E-04,-2.002600114792585E-03},elSize,false);
    Point* p1096 = fluid2 -> addPoint({2.965000021504238E-05,-9.594199946150184E-04},elSize,false);
    Point* p1097 = fluid2 -> addPoint({0.000000000000000E+00,0.000000000000000E+00},elSize,false);
    Point* p1098 = fluid2 -> addPoint({2.965000021504238E-05,9.594199946150184E-04},elSize,false);
    Point* p1099 = fluid2 -> addPoint({1.291899970965460E-04,2.002600114792585E-03},elSize,false);
    Point* p1100 = fluid2 -> addPoint({3.168400144204497E-04,3.131099976599216E-03},elSize,false);
    Point* p1101 = fluid2 -> addPoint({6.129799876362085E-04,4.344400018453598E-03},elSize,false);
    Point* p1102 = fluid2 -> addPoint({1.039199996739626E-03,5.639600101858377E-03},elSize,false);
    Point* p1103 = fluid2 -> addPoint({1.617899979464710E-03,7.011000066995621E-03},elSize,false);
    Point* p1104 = fluid2 -> addPoint({2.369200112298131E-03,8.452000096440315E-03},elSize,false);
    Point* p1105 = fluid2 -> addPoint({3.312600078061223E-03,9.954700246453285E-03},elSize,false);
    Point* p1106 = fluid2 -> addPoint({4.466100130230188E-03,1.151099987328053E-02},elSize,false);
    Point* p1107 = fluid2 -> addPoint({5.847000051289797E-03,1.311200018972158E-02},elSize,false);
    Point* p1108 = fluid2 -> addPoint({7.471399847418070E-03,1.475299987941980E-02},elSize,false);
    Point* p1109 = fluid2 -> addPoint({9.355300106108189E-03,1.642799936234951E-02},elSize,false);
    Point* p1110 = fluid2 -> addPoint({1.151499990373850E-02,1.813299953937531E-02},elSize,false);
    Point* p1111 = fluid2 -> addPoint({1.396999973803759E-02,1.986500062048435E-02},elSize,false);
    Point* p1112 = fluid2 -> addPoint({1.673999987542629E-02,2.162200026214123E-02},elSize,false);
    Point* p1113 = fluid2 -> addPoint({1.984800025820732E-02,2.340300008654594E-02},elSize,false);
    Point* p1114 = fluid2 -> addPoint({2.331699989736080E-02,2.520499937236309E-02},elSize,false);
    Point* p1115 = fluid2 -> addPoint({2.717700041830540E-02,2.702900022268295E-02},elSize,false);
    Point* p1116 = fluid2 -> addPoint({3.145699948072433E-02,2.887200005352497E-02},elSize,false);
    Point* p1117 = fluid2 -> addPoint({3.619199991226196E-02,3.073300048708916E-02},elSize,false);
    Point* p1118 = fluid2 -> addPoint({4.142000153660774E-02,3.260999917984009E-02},elSize,false);
    Point* p1119 = fluid2 -> addPoint({4.718200117349625E-02,3.449999913573265E-02},elSize,false);
    Point* p1120 = fluid2 -> addPoint({5.352399870753288E-02,3.640000149607658E-02},elSize,false);
    Point* p1121 = fluid2 -> addPoint({6.049599871039391E-02,3.830600157380104E-02},elSize,false);
    Point* p1122 = fluid2 -> addPoint({6.815300136804581E-02,4.021200165152550E-02},elSize,false);
    Point* p1123 = fluid2 -> addPoint({7.655599713325500E-02,4.211099818348885E-02},elSize,false);
    Point* p1124 = fluid2 -> addPoint({8.577000349760056E-02,4.399599879980087E-02},elSize,false);
    Point* p1125 = fluid2 -> addPoint({9.586799889802933E-02,4.585599899291992E-02},elSize,false);
    Point* p1126 = fluid2 -> addPoint({1.069300025701523E-01,4.768000170588493E-02},elSize,false);
    Point* p1127 = fluid2 -> addPoint({1.190399974584579E-01,4.945399984717369E-02},elSize,false);
    Point* p1128 = fluid2 -> addPoint({1.322900056838989E-01,5.116200074553490E-02},elSize,false);
    Point* p1129 = fluid2 -> addPoint({1.467899978160858E-01,5.278500169515610E-02},elSize,false);
    Point* p1130 = fluid2 -> addPoint({1.619299948215485E-01,5.423900112509727E-02},elSize,false);
    Point* p1131 = fluid2 -> addPoint({1.770599931478500E-01,5.547099933028221E-02},elSize,false);
    Point* p1132 = fluid2 -> addPoint({1.921699941158295E-01,5.650499835610390E-02},elSize,false);
    Point* p1133 = fluid2 -> addPoint({2.072599977254868E-01,5.735800042748451E-02},elSize,false);
    Point* p1134 = fluid2 -> addPoint({2.223400026559830E-01,5.804999917745590E-02},elSize,false);
    Point* p1135 = fluid2 -> addPoint({2.373899966478348E-01,5.859300121665001E-02},elSize,false);
    Point* p1136 = fluid2 -> addPoint({2.524299919605255E-01,5.899700149893761E-02},elSize,false);
    Point* p1137 = fluid2 -> addPoint({2.674500048160553E-01,5.927500128746033E-02},elSize,false);
    Point* p1138 = fluid2 -> addPoint({2.824699878692627E-01,5.943600088357925E-02},elSize,false);
    Point* p1139 = fluid2 -> addPoint({2.974700033664703E-01,5.948600172996521E-02},elSize,false);
    Point* p1140 = fluid2 -> addPoint({3.124699890613556E-01,5.943400040268898E-02},elSize,false);
    Point* p1141 = fluid2 -> addPoint({3.274700045585632E-01,5.928599834442139E-02},elSize,false);
    Point* p1142 = fluid2 -> addPoint({3.424600064754486E-01,5.904699862003326E-02},elSize,false);
    Point* p1143 = fluid2 -> addPoint({3.574599921703339E-01,5.872299894690514E-02},elSize,false);
    Point* p1144 = fluid2 -> addPoint({3.724600076675415E-01,5.831800028681755E-02},elSize,false);
    Point* p1145 = fluid2 -> addPoint({3.874700069427490E-01,5.783800035715103E-02},elSize,false);
    Point* p1146 = fluid2 -> addPoint({4.024899899959564E-01,5.728600174188614E-02},elSize,false);
    Point* p1147 = fluid2 -> addPoint({4.175300002098083E-01,5.666499957442284E-02},elSize,false);
    Point* p1148 = fluid2 -> addPoint({4.325900077819824E-01,5.597800016403198E-02},elSize,false);
    Point* p1149 = fluid2 -> addPoint({4.476700127124786E-01,5.522900074720383E-02},elSize,false);
    Point* p1150 = fluid2 -> addPoint({4.627699851989746E-01,5.441899970173836E-02},elSize,false);
    Point* p1151 = fluid2 -> addPoint({4.778999984264374E-01,5.355200171470642E-02},elSize,false);
    Point* p1152 = fluid2 -> addPoint({4.930700063705444E-01,5.262900143861771E-02},elSize,false);
    Point* p1153 = fluid2 -> addPoint({5.082700252532959E-01,5.165300145745277E-02},elSize,false);
    Point* p1154 = fluid2 -> addPoint({5.235000252723694E-01,5.062500014901161E-02},elSize,false);
    Point* p1155 = fluid2 -> addPoint({5.387799739837646E-01,4.954700171947479E-02},elSize,false);
    Point* p1156 = fluid2 -> addPoint({5.541099905967712E-01,4.842000082135201E-02},elSize,false);
    Point* p1157 = fluid2 -> addPoint({5.694800019264221E-01,4.724600166082382E-02},elSize,false);
    Point* p1158 = fluid2 -> addPoint({5.849099755287170E-01,4.602500051259995E-02},elSize,false);
    Point* p1159 = fluid2 -> addPoint({6.003900170326233E-01,4.475900158286095E-02},elSize,false);
    Point* p1160 = fluid2 -> addPoint({6.159300208091736E-01,4.344699904322624E-02},elSize,false);
    Point* p1161 = fluid2 -> addPoint({6.315400004386902E-01,4.209100082516670E-02},elSize,false);
    Point* p1162 = fluid2 -> addPoint({6.472100019454956E-01,4.069200158119202E-02},elSize,false);
    Point* p1163 = fluid2 -> addPoint({6.629499793052673E-01,3.924899920821190E-02},elSize,false);
    Point* p1164 = fluid2 -> addPoint({6.787599921226501E-01,3.776200115680695E-02},elSize,false);
    Point* p1165 = fluid2 -> addPoint({6.946499943733215E-01,3.623199835419655E-02},elSize,false);
    Point* p1166 = fluid2 -> addPoint({7.106199860572815E-01,3.465899825096130E-02},elSize,false);
    Point* p1167 = fluid2 -> addPoint({7.266700267791748E-01,3.304199874401093E-02},elSize,false);
    Point* p1168 = fluid2 -> addPoint({7.428100109100342E-01,3.138099983334541E-02},elSize,false);
    Point* p1169 = fluid2 -> addPoint({7.590299844741821E-01,2.967499941587448E-02},elSize,false);
    Point* p1170 = fluid2 -> addPoint({7.753499746322632E-01,2.792399935424328E-02},elSize,false);
    Point* p1171 = fluid2 -> addPoint({7.917699813842773E-01,2.612599916756153E-02},elSize,false);
    Point* p1172 = fluid2 -> addPoint({8.082900047302246E-01,2.428100071847439E-02},elSize,false);
    Point* p1173 = fluid2 -> addPoint({8.249099850654602E-01,2.238800004124641E-02},elSize,false);
    Point* p1174 = fluid2 -> addPoint({8.416299819946289E-01,2.044299989938736E-02},elSize,false);
    Point* p1175 = fluid2 -> addPoint({8.584600090980530E-01,1.844800077378750E-02},elSize,false);
    Point* p1176 = fluid2 -> addPoint({8.745399713516235E-01,1.650499925017357E-02},elSize,false);
    Point* p1177 = fluid2 -> addPoint({8.890900015830994E-01,1.471399981528521E-02},elSize,false);
    Point* p1178 = fluid2 -> addPoint({9.022700190544128E-01,1.306699961423874E-02},elSize,false);
    Point* p1179 = fluid2 -> addPoint({9.142000079154968E-01,1.155200041830540E-02},elSize,false);
    Point* p1180 = fluid2 -> addPoint({9.250000119209290E-01,1.016299985349178E-02},elSize,false);
    Point* p1181 = fluid2 -> addPoint({9.347800016403198E-01,8.889400400221348E-03},elSize,false);
    Point* p1182 = fluid2 -> addPoint({9.436399936676025E-01,7.723500020802021E-03},elSize,false);
    Point* p1183 = fluid2 -> addPoint({9.516599774360657E-01,6.657199934124947E-03},elSize,false);
    Point* p1184 = fluid2 -> addPoint({9.589200019836426E-01,5.682699847966433E-03},elSize,false);
    Point* p1185 = fluid2 -> addPoint({9.654899835586548E-01,4.792999941855669E-03},elSize,false);
    Point* p1186 = fluid2 -> addPoint({9.714400172233582E-01,3.981099929660559E-03},elSize,false);
    Point* p1187 = fluid2 -> addPoint({9.768300056457520E-01,3.240799997001886E-03},elSize,false);
    Point* p1188 = fluid2 -> addPoint({9.817100167274475E-01,2.566199982538819E-03},elSize,false);
    Point* p1189 = fluid2 -> addPoint({9.861299991607666E-01,1.951699960045516E-03},elSize,false);
    Point* p1190 = fluid2 -> addPoint({9.901300072669983E-01,1.392400008626282E-03},elSize,false);
    Point* p1191 = fluid2 -> addPoint({9.937499761581421E-01,8.835400221869349E-04},elSize,false);
    Point* p1192 = fluid2 -> addPoint({9.970300197601318E-01,4.207300080452114E-04},elSize,false);
    
    Line* l1001 = fluid2 -> addLine({p1001,p1002,p1003,p1004,p1005,p1006,p1007,p1008,p1009,p1010,
                                     p1011,p1012,p1013,p1014,p1015,p1016,p1017,p1018,p1019,p1020,
                                     p1021,p1022,p1023,p1024,p1025,p1026,p1027,p1028,p1029,p1030,
                                     p1031,p1032,p1033,p1034,p1035,p1036,p1037,p1038,p1039,p1040,
                                     p1041,p1042,p1043,p1044,p1045,p1046,p1047,p1048,p1049,p1050,
                                     p1051,p1052,p1053,p1054,p1055,p1056,p1057});

    Line* l1002 = fluid2 -> addLine({p1057,p1058,p1059,p1060,p1061,p1062,p1063,p1064,p1065,p1066,
                                    p1067,p1068,p1069,p1070,p1071,p1072,p1073,p1074,p1075,p1076,
                                    p1077,p1078,p1079,p1080,p1081,p1082,p1083,p1084,p1085,p1086,
                                    p1087,p1088,p1089,p1090,p1091,p1092,p1093,p1094,p1095,p1096,
                                    p1097});

    Line* l1003 = fluid2 -> addLine({p1097,p1098,p1099,p1100,p1101,p1102,p1103,p1104,p1105,p1106,
                                    p1107,p1108,p1109,p1110,p1111,p1112,p1113,p1114,p1115,p1116,
                                    p1117,p1118,p1119,p1120,p1121,p1122,p1123,p1124,p1125,p1126,
                                    p1127,p1128,p1129,p1130,p1131,p1132,p1133,p1134,p1135,p1136,
                                    p1137});

    Line* l1004 = fluid2 -> addLine({p1137,p1138,p1139,p1140,p1141,p1142,p1143,p1144,p1145,p1146,
                                    p1147,p1148,p1149,p1150,p1151,p1152,p1153,p1154,p1155,p1156,
                                    p1157,p1158,p1159,p1160,p1161,p1162,p1163,p1164,p1165,p1166,
                                    p1167,p1168,p1169,p1170,p1171,p1172,p1173,p1174,p1175,p1176,
                                    p1177,p1178,p1179,p1180,p1181,p1182,p1183,p1184,p1185,p1186,
                                    p1187,p1188,p1189,p1190,p1191,p1192,p1001});


    double elSize2 = 0.075;
    double elSize3 = 0.03;
    double elSize4 = 0.075;

    Point* p1386 = fluid2 -> addPoint({1.25,0.0},elSize2,false);
    Point* p1387 = fluid2 -> addPoint({.4,0.},elSize2,false);
    Point* p1388 = fluid2 -> addPoint({-.4,0.},elSize2,false);
    Point* p1389 = fluid2 -> addPoint({0,.4},elSize2,false);
    Point* p1390 = fluid2 -> addPoint({1,.4},elSize4,false);
    Point* p1391 = fluid2 -> addPoint({1.4,0.},elSize3,false);
    Point* p1392 = fluid2 -> addPoint({1.,-.4},elSize4,false);
    Point* p1393 = fluid2 -> addPoint({0.,-.4},elSize2,false);
    
    Point* p1394 = fluid2 -> addPoint({0,-0.6},elSize2,false);
    Point* p1395 = fluid2 -> addPoint({-0.6,0.},elSize2,false);
    Point* p1396 = fluid2 -> addPoint({0.,0.6},elSize2,false);
    Point* p1397 = fluid2 -> addPoint({1.,.6},elSize4,false);
    Point* p1398 = fluid2 -> addPoint({1.6,.0},elSize4,false);
    Point* p1399 = fluid2 -> addPoint({1.,-0.6},elSize4,false);
    
    Line* l1014 = fluid2 -> addLine({p1393,p1097,p1388});
    Line* l1015 = fluid2 -> addLine({p1388,p1097,p1389});
    Line* l1016 = fluid2 -> addLine({p1390,p1001,p1391});
    Line* l1017 = fluid2 -> addLine({p1391,p1001,p1392});
    Line* l1018 = fluid2 -> addLine({p1389,p1390});
    Line* l1019 = fluid2 -> addLine({p1392,p1393});
    Line* l1020 = fluid2 -> addLine({p1001,p1386});
    Line* l1021 = fluid2 -> addLine({p1386,p1391});

    Line* l1022 = fluid2 -> addLine({p1394,p1097,p1395});
    Line* l1023 = fluid2 -> addLine({p1395,p1097,p1396});
    Line* l1024 = fluid2 -> addLine({p1396,p1397});
    Line* l1025 = fluid2 -> addLine({p1397,p1001,p1398});
    Line* l1026 = fluid2 -> addLine({p1398,p1001,p1399});
    Line* l1027 = fluid2 -> addLine({p1399,p1394});

    
    LineLoop* ll1005 = fluid2->addLineLoop({ l1017 -> operator-(), l1019 -> operator-(), 
                                             l1014 -> operator-(), l1015 -> operator-(), 
                                             l1018 -> operator-(), l1016 -> operator-(), 
                                             l001, l002, l003, l004 });
   
    LineLoop* ll1006 = fluid2->addLineLoop({ l1027 -> operator-(), l1026 -> operator-(), 
                                             l1025 -> operator-(), l1024 -> operator-(),
                                             l1023 -> operator-(), l1022 -> operator-(),
                                             l1019, l1017, l1016, l1018, l1015, l1014}); 
    

    //Transfinite lines 
    int t1 = 20; int t2 = 50; int t3 = 25; int t4 = 4;
    //corners
    fluid2 -> transfiniteLine({ l1001 }, t2, 1.02);
    fluid2 -> transfiniteLine({ l1002 -> operator-()}, t3, 1.05);
    fluid2 -> transfiniteLine({ l1003 }, t3, 1.05);
    fluid2 -> transfiniteLine({ l1004 -> operator-()}, t2, 1.02);
    
    PlaneSurface* s1005 = fluid2 -> addPlaneSurface({ll1005});
    PlaneSurface* s1006 = fluid2 -> addPlaneSurface({ll1006});
   
    fluid2 -> addBoundaryCondition("NEUMANN", l1020, {}, {}, "GLOBAL");
    fluid2 -> addBoundaryCondition("NEUMANN", l1021, {}, {}, "GLOBAL");


    fluid2 -> addBoundaryCondition("NEUMANN", l1014, {}, {}, "GLOBAL");
    fluid2 -> addBoundaryCondition("NEUMANN", l1015, {}, {}, "GLOBAL");
    fluid2 -> addBoundaryCondition("NEUMANN", l1016, {}, {}, "GLOBAL");
    fluid2 -> addBoundaryCondition("NEUMANN", l1017, {}, {}, "GLOBAL");
    fluid2 -> addBoundaryCondition("NEUMANN", l1018, {}, {}, "GLOBAL");
    fluid2 -> addBoundaryCondition("NEUMANN", l1019, {}, {}, "GLOBAL");

    fluid2 -> addBoundaryCondition("GLUE", l1022, {0}, {0}, "GLOBAL");
    fluid2 -> addBoundaryCondition("GLUE", l1023, {0}, {0}, "GLOBAL");
    fluid2 -> addBoundaryCondition("GLUE", l1024, {0}, {0}, "GLOBAL");
    fluid2 -> addBoundaryCondition("GLUE", l1025, {0}, {0}, "GLOBAL");
    fluid2 -> addBoundaryCondition("GLUE", l1026, {0}, {0}, "GLOBAL");
    fluid2 -> addBoundaryCondition("GLUE", l1027, {0}, {0}, "GLOBAL");
    
    fluid2 -> addBoundaryCondition("DIRICHLET", l1001, {0.0}, {0.0}, "GLOBAL");
    fluid2 -> addBoundaryCondition("DIRICHLET", l1002, {0}, {0.0}, "GLOBAL");
    fluid2 -> addBoundaryCondition("DIRICHLET", l1003, {0}, {0}, "GLOBAL");
    fluid2 -> addBoundaryCondition("DIRICHLET", l1004, {0}, {0}, "GLOBAL");
   

    if (rank == 0){
  
        FluidDomain* problem = new FluidDomain(fluid1);
        problem -> addSurfaceMaterial({ s001,s002,s003,s004,s005,s006,s007,s008,s009 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
        problem -> generateMesh("T6", "DELAUNAY", "coarse", "", false, false);

        FluidDomain* problem2 = new FluidDomain(fluid2);
        problem2 -> addSurfaceMaterial({ s1005, s1006 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
        problem2 -> generateMesh("T6", "DELAUNAY", "fine", "", false, false);

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
     
    //coarseModel.solveTransientProblem(1, 1.e-15, 2);
 
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
      
    coarseModel.readInitialValues("saidaVelCoarse106155.txt","saidaVelCoarse106156.txt");
    fineModel.readInitialValues("saidaVelFine106155.txt","saidaVelFine106156.txt");


       
     arlequinProblem.setFluidModels(coarseModel, fineModel) ; 
     arlequinProblem.solveArlequinProblem(5, 1.e-7, 2, 1); 
     arlequinProblem.solveArlequinProblemMoving(5, 1.e-7, 2, 1); 
          
     // 
    //Finalize main program   
    PetscFinalize();
 
    return 0; 
}
 
 
  




 
