
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

        double x1 = -9.5;
        double x2 = -0.6;
        double x3 = 1.6;
        double x4 = 20.5;

        double y1 = -10;
        double y2 = -0.85;
        double y3 = 0.85;
        double y4 = 10;
    
        double h1 = 25; double h2 = 51; double h3 = 41;
        double v1 = 25; double v2 = 29;

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

        double prog1 = 1.03;
        double prog2 = 1.05;

        fluid1 -> transfiniteLine({l001}, h1, 1/prog2);
        fluid1 -> transfiniteLine({l002}, 7);
        fluid1 -> transfiniteLine({l003}, h3, prog1);
        fluid1 -> transfiniteLine({l004}, v1, 1/prog1);
        fluid1 -> transfiniteLine({l005}, v1, 1/prog1);
        fluid1 -> transfiniteLine({l006}, v1, 1/prog1);
        fluid1 -> transfiniteLine({l007}, v1, 1/prog1);
        fluid1 -> transfiniteLine({l008}, h1, 1/prog2);
        fluid1 -> transfiniteLine({l009}, h2);
        fluid1 -> transfiniteLine({l010}, 80, 1.02);
        fluid1 -> transfiniteLine({l011}, 7);
        fluid1 -> transfiniteLine({l012}, v2);
        fluid1 -> transfiniteLine({l013}, v2);
        fluid1 -> transfiniteLine({l014}, 7);
        fluid1 -> transfiniteLine({l015}, h1, 1/prog2);
        fluid1 -> transfiniteLine({l016}, h2);
        fluid1 -> transfiniteLine({l017}, 80, 1.02);
        fluid1 -> transfiniteLine({l018}, v1, prog1);
        fluid1 -> transfiniteLine({l019}, v1, prog1);
        fluid1 -> transfiniteLine({l020}, v1, prog1);
        fluid1 -> transfiniteLine({l021}, v1, prog1);
        fluid1 -> transfiniteLine({l022}, h1, 1/prog2); 
        fluid1 -> transfiniteLine({l023}, 7);
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

        fluid1 -> addBoundaryCondition("NEUMANN", l005, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l006, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l008, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l009, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l010, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l012, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l013, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l015, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l016, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l017, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l019, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l020, {}, {}, {}, "GLOBAL");

        fluid1 -> addBoundaryCondition("DIRICHLET", l001, {}, {0.0}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l002, {}, {0.0}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l003, {}, {0.0}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l022, {}, {0.0}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l023, {}, {0.0}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l024, {}, {0.0}, {}, "GLOBAL");

        fluid1 -> addBoundaryCondition("NEUMANN", l007, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l014, {}, {}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("NEUMANN", l021, {}, {}, {}, "GLOBAL");

        fluid1 -> addBoundaryCondition("DIRICHLET", l004, {1.0}, {0.0}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l011, {1.0}, {0.0}, {}, "GLOBAL");
        fluid1 -> addBoundaryCondition("DIRICHLET", l018, {1.0}, {0.0}, {}, "GLOBAL");

        
        // Fine
        Geometry* fluid2 = new Geometry(1);
        double elSize = 1.e-2;

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


        Point* p1193 = fluid2 -> addPoint({1.0021039076 ,   -0.0148517195},elSize,false);
        Point* p1194 = fluid2 -> addPoint({0.9991296319 ,   -0.0152730574},elSize,false);
        Point* p1195 = fluid2 -> addPoint({0.9958416605 ,   -0.015736986},elSize,false);
        Point* p1196 = fluid2 -> addPoint({0.9922123312 ,   -0.0162471611},elSize,false);
        Point* p1197 = fluid2 -> addPoint({0.9882010663 ,   -0.0168080347},elSize,false);
        Point* p1198 = fluid2 -> addPoint({0.9837695207 ,   -0.0174241421},elSize,false);
        Point* p1199 = fluid2 -> addPoint({0.9788772213 ,   -0.0181004402},elSize,false);
        Point* p1200 = fluid2 -> addPoint({0.9734742202 ,   -0.0188425271},elSize,false);
        Point* p1201 = fluid2 -> addPoint({0.9675100674 ,   -0.0196563529},elSize,false);
        Point* p1202 = fluid2 -> addPoint({0.9609237881 ,   -0.0205482587},elSize,false);
        Point* p1203 = fluid2 -> addPoint({0.9536457466 ,   -0.0215251763},elSize,false);
        Point* p1204 = fluid2 -> addPoint({0.9456064673 ,   -0.0225940407},elSize,false);
        Point* p1205 = fluid2 -> addPoint({0.936726521  ,   -0.0237625662},elSize,false);
        Point* p1206 = fluid2 -> addPoint({0.926924642  ,   -0.0250390141},elSize,false);
        Point* p1207 = fluid2 -> addPoint({0.9161009667 ,   -0.026431058},elSize,false);
        Point* p1208 = fluid2 -> addPoint({0.9041441065 ,   -0.0279494657},elSize,false);
        Point* p1209 = fluid2 -> addPoint({0.8909355865 ,   -0.0296000274},elSize,false);
        Point* p1210 = fluid2 -> addPoint({0.8763551336 ,   -0.0313947671},elSize,false);        
        Point* p1211 = fluid2 -> addPoint({0.8602421914 ,   -0.033341752},elSize,false);
        Point* p1212 = fluid2 -> addPoint({0.8433795075 ,   -0.0353406225},elSize,false);
        Point* p1213 = fluid2 -> addPoint({0.8266254123 ,   -0.0372895875},elSize,false);
        Point* p1214 = fluid2 -> addPoint({0.8099712562 ,   -0.0391864827},elSize,false);
        Point* p1215 = fluid2 -> addPoint({0.7934188514 ,   -0.041035098},elSize,false);
        Point* p1216 = fluid2 -> addPoint({0.7769664979 ,   -0.0428366401},elSize,false);
        Point* p1217 = fluid2 -> addPoint({0.7606146459 ,   -0.0445910595},elSize,false);
        Point* p1218 = fluid2 -> addPoint({0.7443623555 ,   -0.0463004578},elSize,false);
        Point* p1219 = fluid2 -> addPoint({0.728189659  ,   -0.0479648241},elSize,false);
        Point* p1220 = fluid2 -> addPoint({0.7121070017 ,   -0.0495851091},elSize,false);
        Point* p1221 = fluid2 -> addPoint({0.6961040387 ,   -0.0511613572},elSize,false);
        Point* p1222 = fluid2 -> addPoint({0.6801811724 ,   -0.0526945243},elSize,false);
        Point* p1223 = fluid2 -> addPoint({0.6643370383 ,   -0.0541847304},elSize,false);
        Point* p1224 = fluid2 -> addPoint({0.648561693  ,   -0.0556309752},elSize,false);
        Point* p1225 = fluid2 -> addPoint({0.6328560372 ,   -0.0570331575},elSize,false);
        Point* p1226 = fluid2 -> addPoint({0.6172100864 ,   -0.0583922803},elSize,false);
        Point* p1227 = fluid2 -> addPoint({0.6016323481 ,   -0.0597074668},elSize,false);
        Point* p1228 = fluid2 -> addPoint({0.5861129783 ,   -0.0609766822},elSize,false);
        Point* p1229 = fluid2 -> addPoint({0.5706428883 ,   -0.0622008569},elSize,false);
        Point* p1230 = fluid2 -> addPoint({0.5552311121 ,   -0.063378045},elSize,false);
        Point* p1231 = fluid2 -> addPoint({0.5398577115 ,   -0.0645082344},elSize,false);
        Point* p1232 = fluid2 -> addPoint({0.5245329672 ,   -0.0655893922},elSize,false);
        Point* p1233 = fluid2 -> addPoint({0.509255731  ,   -0.0666205792},elSize,false);
        Point* p1234 = fluid2 -> addPoint({0.4940061096 ,   -0.0675997633},elSize,false);
        Point* p1235 = fluid2 -> addPoint({0.4787845933 ,   -0.0685258954},elSize,false);
        Point* p1236 = fluid2 -> addPoint({0.4636008281 ,   -0.069395972},elSize,false);
        Point* p1237 = fluid2 -> addPoint({0.4484438327 ,   -0.0702090275},elSize,false);
        Point* p1238 = fluid2 -> addPoint({0.433303862  ,   -0.0709610043},elSize,false);
        Point* p1239 = fluid2 -> addPoint({0.4181812133 ,   -0.071650857},elSize,false);
        Point* p1240 = fluid2 -> addPoint({0.4030748747 ,   -0.0722745944},elSize,false);
        Point* p1241 = fluid2 -> addPoint({0.3879851869 ,   -0.0728291507},elSize,false);
        Point* p1242 = fluid2 -> addPoint({0.3729021671 ,   -0.073311482},elSize,false);
        Point* p1243 = fluid2 -> addPoint({0.3578243838 ,   -0.0737185723},elSize,false);
        Point* p1244 = fluid2 -> addPoint({0.3427415505 ,   -0.0740443562},elSize,false);
        Point* p1245 = fluid2 -> addPoint({0.3276635537 ,   -0.0742847496},elSize,false);
        Point* p1246 = fluid2 -> addPoint({0.3125699885 ,   -0.0744336671},elSize,false);
        Point* p1247 = fluid2 -> addPoint({0.2974710036 ,   -0.0744860017},elSize,false);
        Point* p1248 = fluid2 -> addPoint({0.2823645605 ,   -0.0744356304},elSize,false);
        Point* p1249 = fluid2 -> addPoint({0.2672308208 ,   -0.0742733998},elSize,false);
        Point* p1250 = fluid2 -> addPoint({0.2520897602 ,   -0.0739931424},elSize,false);
        Point* p1251 = fluid2 -> addPoint({0.2369181452 ,   -0.0735855779},elSize,false);
        Point* p1252 = fluid2 -> addPoint({0.2217256824 ,   -0.0730374143},elSize,false);
        Point* p1253 = fluid2 -> addPoint({0.2064928556 ,   -0.0723383707},elSize,false);
        Point* p1254 = fluid2 -> addPoint({0.1912345663 ,   -0.0714758024},elSize,false);
        Point* p1255 = fluid2 -> addPoint({0.175939129  ,   -0.0704290628},elSize,false);
        Point* p1256 = fluid2 -> addPoint({0.1606041844 ,   -0.0691802938},elSize,false);
        Point* p1257 = fluid2 -> addPoint({0.1452411345 ,   -0.0677048215},elSize,false);
        Point* p1258 = fluid2 -> addPoint({0.1305022974 ,   -0.0660550897},elSize,false);
        Point* p1259 = fluid2 -> addPoint({0.1169995987 ,   -0.0643145777},elSize,false);
        Point* p1260 = fluid2 -> addPoint({0.1046284773 ,   -0.0625023828},elSize,false);
        Point* p1261 = fluid2 -> addPoint({0.0932951677 ,   -0.0606337031},elSize,false);
        Point* p1262 = fluid2 -> addPoint({0.0829143876 ,   -0.0587216722},elSize,false);
        Point* p1263 = fluid2 -> addPoint({0.0734059574 ,   -0.0567765104},elSize,false);
        Point* p1264 = fluid2 -> addPoint({0.0646948759 ,   -0.0548079384},elSize,false);
        Point* p1265 = fluid2 -> addPoint({0.0567136293 ,   -0.0528212929},elSize,false);
        Point* p1266 = fluid2 -> addPoint({0.0494011757 ,   -0.0508222874},elSize,false);
        Point* p1267 = fluid2 -> addPoint({0.0426997929 ,   -0.0488146702},elSize,false);
        Point* p1268 = fluid2 -> addPoint({0.0365562941 ,   -0.0467995851},elSize,false);
        Point* p1269 = fluid2 -> addPoint({0.0309228332 ,   -0.0447770698},elSize,false);
        Point* p1270 = fluid2 -> addPoint({0.0257563583 ,   -0.0427465339},elSize,false);
        Point* p1271 = fluid2 -> addPoint({0.0210159407 ,   -0.0407053061},elSize,false);
        Point* p1272 = fluid2 -> addPoint({0.0166653558 ,   -0.0386495385},elSize,false);
        Point* p1273 = fluid2 -> addPoint({0.0126720964 ,   -0.036575183},elSize,false);
        Point* p1274 = fluid2 -> addPoint({0.0090045669 ,   -0.034473579},elSize,false);
        Point* p1275 = fluid2 -> addPoint({0.0056401317 ,   -0.0323395066},elSize,false);
        Point* p1276 = fluid2 -> addPoint({0.002555131  ,   -0.0301629933},elSize,false);
        Point* p1277 = fluid2 -> addPoint({-0.0002648114    ,   -0.0279368418},elSize,false);
        Point* p1278 = fluid2 -> addPoint({-0.0028322301    ,   -0.0256541562},elSize,false);
        Point* p1279 = fluid2 -> addPoint({-0.005153576 ,   -0.0233094178},elSize,false);
        Point* p1280 = fluid2 -> addPoint({-0.0072314457    ,   -0.0209007507},elSize,false);
        Point* p1281 = fluid2 -> addPoint({-0.0090596527    ,   -0.0184357002},elSize,false);
        Point* p1282 = fluid2 -> addPoint({-0.0106304376    ,   -0.0159359441},elSize,false);
        Point* p1283 = fluid2 -> addPoint({-0.0119422488    ,   -0.0134236723},elSize,false);
        Point* p1284 = fluid2 -> addPoint({-0.0129971616    ,   -0.0109292649},elSize,false);
        Point* p1285 = fluid2 -> addPoint({-0.0138012787    ,   -0.0084952006},elSize,false);
        Point* p1286 = fluid2 -> addPoint({-0.0143729621    ,   -0.0061658508},elSize,false);
        Point* p1287 = fluid2 -> addPoint({-0.0147413436    ,   -0.0039691278},elSize,false);
        Point* p1288 = fluid2 -> addPoint({-0.0149392344    ,   -0.0019250797},elSize,false);
        Point* p1289 = fluid2 -> addPoint({-0.015   ,  0},elSize,false);
        Point* p1290 = fluid2 -> addPoint({-0.0149392344    ,   0.0019250797},elSize,false);
        Point* p1291 = fluid2 -> addPoint({-0.0147413436    ,   0.0039691278},elSize,false);
        Point* p1292 = fluid2 -> addPoint({-0.0143729621    ,   0.0061658508},elSize,false);
        Point* p1293 = fluid2 -> addPoint({-0.0138012787    ,   0.0084952006},elSize,false);
        Point* p1294 = fluid2 -> addPoint({-0.0129971616    ,   0.0109292649},elSize,false);
        Point* p1295 = fluid2 -> addPoint({-0.0119422488    ,   0.0134236723},elSize,false);
        Point* p1296 = fluid2 -> addPoint({-0.0106304376    ,   0.0159359441},elSize,false);
        Point* p1297 = fluid2 -> addPoint({-0.0090596527    ,   0.0184357002},elSize,false);
        Point* p1298 = fluid2 -> addPoint({-0.0072314457    ,   0.0209007507},elSize,false);
        Point* p1299 = fluid2 -> addPoint({-0.005153576 ,   0.0233094178},elSize,false);
        Point* p1300 = fluid2 -> addPoint({-0.0028322301    ,   0.0256541562},elSize,false);
        Point* p1301 = fluid2 -> addPoint({-0.0002648114    ,   0.0279368418},elSize,false);
        Point* p1302 = fluid2 -> addPoint({0.002555131  ,   0.0301629933},elSize,false);
        Point* p1303 = fluid2 -> addPoint({0.0056401317 ,   0.0323395066},elSize,false);
        Point* p1304 = fluid2 -> addPoint({0.0090045669 ,   0.034473579},elSize,false);
        Point* p1305 = fluid2 -> addPoint({0.0126720964 ,   0.036575183},elSize,false);
        Point* p1306 = fluid2 -> addPoint({0.0166653558 ,   0.0386495385},elSize,false);
        Point* p1307 = fluid2 -> addPoint({0.0210159407 ,   0.0407053061},elSize,false);
        Point* p1308 = fluid2 -> addPoint({0.0257563583 ,   0.0427465339},elSize,false);
        Point* p1309 = fluid2 -> addPoint({0.0309228332 ,   0.0447770698},elSize,false);
        Point* p1310 = fluid2 -> addPoint({0.0365562941 ,   0.0467995851},elSize,false);
        Point* p1311 = fluid2 -> addPoint({0.0426997929 ,   0.0488146702},elSize,false);
        Point* p1312 = fluid2 -> addPoint({0.0494011757 ,   0.0508222874},elSize,false);
        Point* p1313 = fluid2 -> addPoint({0.0567136293 ,   0.0528212929},elSize,false);
        Point* p1314 = fluid2 -> addPoint({0.0646948759 ,   0.0548079384},elSize,false);
        Point* p1315 = fluid2 -> addPoint({0.0734059574 ,   0.0567765104},elSize,false);
        Point* p1316 = fluid2 -> addPoint({0.0829143876 ,   0.0587216722},elSize,false);
        Point* p1317 = fluid2 -> addPoint({0.0932951677 ,   0.0606337031},elSize,false);
        Point* p1318 = fluid2 -> addPoint({0.1046284773 ,   0.0625023828},elSize,false);
        Point* p1319 = fluid2 -> addPoint({0.1169995987 ,   0.0643145777},elSize,false);
        Point* p1320 = fluid2 -> addPoint({0.1305022974 ,   0.0660550897},elSize,false);
        Point* p1321 = fluid2 -> addPoint({0.1452411345 ,   0.0677048215},elSize,false);
        Point* p1322 = fluid2 -> addPoint({0.1606041844 ,   0.0691802938},elSize,false);
        Point* p1323 = fluid2 -> addPoint({0.175939129  ,   0.0704290628},elSize,false);
        Point* p1324 = fluid2 -> addPoint({0.1912345663 ,   0.0714758024},elSize,false);
        Point* p1325 = fluid2 -> addPoint({0.2064928556 ,   0.0723383707},elSize,false);
        Point* p1326 = fluid2 -> addPoint({0.2217256824 ,   0.0730374143},elSize,false);
        Point* p1327 = fluid2 -> addPoint({0.2369181452 ,   0.0735855779},elSize,false);
        Point* p1328 = fluid2 -> addPoint({0.2520897602 ,   0.0739931424},elSize,false);
        Point* p1329 = fluid2 -> addPoint({0.2672308208 ,   0.0742733998},elSize,false);
        Point* p1330 = fluid2 -> addPoint({0.2823645605 ,   0.0744356304},elSize,false);
        Point* p1331 = fluid2 -> addPoint({0.2974710036 ,   0.0744860017},elSize,false);
        Point* p1332 = fluid2 -> addPoint({0.3125699885 ,   0.0744336671},elSize,false);
        Point* p1333 = fluid2 -> addPoint({0.3276635537 ,   0.0742847496},elSize,false);
        Point* p1334 = fluid2 -> addPoint({0.3427415505 ,   0.0740443562},elSize,false);
        Point* p1335 = fluid2 -> addPoint({0.3578243838 ,   0.0737185723},elSize,false);
        Point* p1336 = fluid2 -> addPoint({0.3729021671 ,   0.073311482},elSize,false);
        Point* p1337 = fluid2 -> addPoint({0.3879851869 ,   0.0728291507},elSize,false);
        Point* p1338 = fluid2 -> addPoint({0.4030748747 ,   0.0722745944},elSize,false);
        Point* p1339 = fluid2 -> addPoint({0.4181812133 ,   0.071650857},elSize,false);
        Point* p1340 = fluid2 -> addPoint({0.433303862  ,   0.0709610043},elSize,false);
        Point* p1341 = fluid2 -> addPoint({0.4484438327 ,   0.0702090275},elSize,false);
        Point* p1342 = fluid2 -> addPoint({0.4636008281 ,   0.069395972},elSize,false);
        Point* p1343 = fluid2 -> addPoint({0.4787845933 ,   0.0685258954},elSize,false);
        Point* p1344 = fluid2 -> addPoint({0.4940061096 ,   0.0675997633},elSize,false);
        Point* p1345 = fluid2 -> addPoint({0.509255731  ,   0.0666205792},elSize,false);
        Point* p1346 = fluid2 -> addPoint({0.5245329672 ,   0.0655893922},elSize,false);
        Point* p1347 = fluid2 -> addPoint({0.5398577115 ,   0.0645082344},elSize,false);
        Point* p1348 = fluid2 -> addPoint({0.5552311121 ,   0.063378045},elSize,false);
        Point* p1349 = fluid2 -> addPoint({0.5706428883 ,   0.0622008569},elSize,false);
        Point* p1350 = fluid2 -> addPoint({0.5861129783 ,   0.0609766822},elSize,false);
        Point* p1351 = fluid2 -> addPoint({0.6016323481 ,   0.0597074668},elSize,false);
        Point* p1352 = fluid2 -> addPoint({0.6172100864 ,   0.0583922803},elSize,false);
        Point* p1353 = fluid2 -> addPoint({0.6328560372 ,   0.0570331575},elSize,false);
        Point* p1354 = fluid2 -> addPoint({0.648561693  ,   0.0556309752},elSize,false);
        Point* p1355 = fluid2 -> addPoint({0.6643370383 ,   0.0541847304},elSize,false);
        Point* p1356 = fluid2 -> addPoint({0.6801811724 ,   0.0526945243},elSize,false);
        Point* p1357 = fluid2 -> addPoint({0.6961040387 ,   0.0511613572},elSize,false);
        Point* p1358 = fluid2 -> addPoint({0.7121070017 ,   0.0495851091},elSize,false);
        Point* p1359 = fluid2 -> addPoint({0.728189659  ,   0.0479648241},elSize,false);
        Point* p1360 = fluid2 -> addPoint({0.7443623555 ,   0.0463004578},elSize,false);
        Point* p1361 = fluid2 -> addPoint({0.7606146459 ,   0.0445910595},elSize,false);
        Point* p1362 = fluid2 -> addPoint({0.7769664979 ,   0.0428366401},elSize,false);
        Point* p1363 = fluid2 -> addPoint({0.7934188514 ,   0.041035098},elSize,false);
        Point* p1364 = fluid2 -> addPoint({0.8099712562 ,   0.0391864827},elSize,false);
        Point* p1365 = fluid2 -> addPoint({0.8266254123 ,   0.0372895875},elSize,false);
        Point* p1366 = fluid2 -> addPoint({0.8433795075 ,   0.0353406225},elSize,false);
        Point* p1367 = fluid2 -> addPoint({0.8602421914 ,   0.033341752},elSize,false);
        Point* p1368 = fluid2 -> addPoint({0.8763551336 ,   0.0313947671},elSize,false);
        Point* p1369 = fluid2 -> addPoint({0.8909355865 ,   0.0296000274},elSize,false);
        Point* p1370 = fluid2 -> addPoint({0.9041441065 ,   0.0279494657},elSize,false);
        Point* p1371 = fluid2 -> addPoint({0.9161009667 ,   0.026431058},elSize,false);
        Point* p1372 = fluid2 -> addPoint({0.926924642  ,   0.0250390141},elSize,false);
        Point* p1373 = fluid2 -> addPoint({0.936726521  ,   0.0237625662},elSize,false);
        Point* p1374 = fluid2 -> addPoint({0.9456064673 ,   0.0225940407},elSize,false);
        Point* p1375 = fluid2 -> addPoint({0.9536457466 ,   0.0215251763},elSize,false);
        Point* p1376 = fluid2 -> addPoint({0.9609237881 ,   0.0205482587},elSize,false);
        Point* p1377 = fluid2 -> addPoint({0.9675100674 ,   0.0196563529},elSize,false);
        Point* p1378 = fluid2 -> addPoint({0.9734742202 ,   0.0188425271},elSize,false);
        Point* p1379 = fluid2 -> addPoint({0.9788772213 ,   0.0181004402},elSize,false);
        Point* p1380 = fluid2 -> addPoint({0.9837695207 ,   0.0174241421},elSize,false);
        Point* p1381 = fluid2 -> addPoint({0.9882010663 ,   0.0168080347},elSize,false);
        Point* p1382 = fluid2 -> addPoint({0.9922123312 ,   0.0162471611},elSize,false);
        Point* p1383 = fluid2 -> addPoint({0.9958416605 ,   0.015736986},elSize,false);
        Point* p1384 = fluid2 -> addPoint({0.9991257419 ,   0.0152736068},elSize,false);
        Point* p1385 = fluid2 -> addPoint({1.0021039076 ,   0.0148517195},elSize,false);
    
        Line* l1005 = fluid2 -> addLine({p1193,p1194,p1195,p1196,p1197,p1198,p1199,p1200,p1201,p1202,
                                        p1203,p1204,p1205,p1206,p1207,p1208,p1209,p1210,p1211,p1212,
                                        p1213,p1214,p1215,p1216,p1217,p1218,p1219,p1220,p1221,p1222,
                                        p1223,p1224,p1225,p1226,p1227,p1228,p1229,p1230,p1231,p1232,
                                        p1233,p1234,p1235,p1236,p1237,p1238,p1239,p1240,p1241,p1242,
                                        p1243,p1244,p1245,p1246,p1247,p1248,p1249});

        Line* l1006 = fluid2 -> addLine({p1249,p1250,p1251,p1252,p1253,p1254,p1255,p1256,p1257,p1258,
                                        p1259,p1260,p1261,p1262,p1263,p1264,p1265,p1266,p1267,p1268,
                                        p1269,p1270,p1271,p1272,p1273,p1274,p1275,p1276,p1277,p1278,
                                        p1279,p1280,p1281,p1282,p1283,p1284,p1285,p1286,p1287,p1288,
                                        p1289});
        
        Line* l1007 = fluid2 -> addLine({p1289,p1290,p1291,p1292,p1293,p1294,p1295,p1296,p1297,p1298,
                                        p1299,p1300,p1301,p1302,p1303,p1304,p1305,p1306,p1307,p1308,
                                        p1309,p1310,p1311,p1312,p1313,p1314,p1315,p1316,p1317,p1318,
                                        p1319,p1320,p1321,p1322,p1323,p1324,p1325,p1326,p1327,p1328,
                                        p1329});

        Line* l1008 = fluid2 -> addLine({p1329,p1330,p1331,p1332,p1333,p1334,p1335,p1336,p1337,p1338,
                                        p1339,p1340,p1341,p1342,p1343,p1344,p1345,p1346,p1347,p1348,
                                        p1349,p1350,p1351,p1352,p1353,p1354,p1355,p1356,p1357,p1358,
                                        p1359,p1360,p1361,p1362,p1363,p1364,p1365,p1366,p1367,p1368,
                                        p1369,p1370,p1371,p1372,p1373,p1374,p1375,p1376,p1377,p1378,
                                        p1379,p1380,p1381,p1382,p1383,p1384,p1385});


        Line* l1009 = fluid2 -> addLine({p1001,p1193});
        Line* l1010 = fluid2 -> addLine({p1057,p1249});
        Line* l1011 = fluid2 -> addLine({p1097,p1289});
        Line* l1012 = fluid2 -> addLine({p1137,p1329});
        Line* l1013 = fluid2 -> addLine({p1001,p1385});

        LineLoop* ll1001 = fluid2 -> addLineLoop({ l1009 -> operator-(), l1005 -> operator-(), l1010, l1001 });
        LineLoop* ll1002 = fluid2 -> addLineLoop({ l1010 -> operator-(), l1006 -> operator-(), l1011, l1002 });
        LineLoop* ll1003 = fluid2 -> addLineLoop({ l1011 -> operator-(), l1007 -> operator-(), l1012, l1003 });
        LineLoop* ll1004 = fluid2 -> addLineLoop({ l1012 -> operator-(), l1008 -> operator-(), l1013, l1004 });

        double elSize2 = 0.075/2.;
        double elSize3 = 0.03/2.;
        double elSize4 = 0.075/2.;

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
        
        Line* l1014 = fluid2 -> addCircle({p1393,p1097,p1388});
        Line* l1015 = fluid2 -> addCircle({p1388,p1097,p1389});
        Line* l1016 = fluid2 -> addCircle({p1390,p1001,p1391});
        Line* l1017 = fluid2 -> addCircle({p1391,p1001,p1392});
        Line* l1018 = fluid2 -> addLine({p1389,p1390});
        Line* l1019 = fluid2 -> addLine({p1392,p1393});
        Line* l1020 = fluid2 -> addLine({p1001,p1386});
        Line* l1021 = fluid2 -> addLine({p1386,p1391});

        Line* l1022 = fluid2 -> addCircle({p1394,p1097,p1395});
        Line* l1023 = fluid2 -> addCircle({p1395,p1097,p1396});
        Line* l1024 = fluid2 -> addLine({p1396,p1397});
        Line* l1025 = fluid2 -> addCircle({p1397,p1001,p1398});
        Line* l1026 = fluid2 -> addCircle({p1398,p1001,p1399});
        Line* l1027 = fluid2 -> addLine({p1399,p1394});

        
        LineLoop* ll1005 = fluid2->addLineLoop({ l1017 -> operator-(), l1019 -> operator-(), 
                                                l1014 -> operator-(), l1015 -> operator-(), 
                                                l1018 -> operator-(), l1016 -> operator-(), 
                                                l1020, l1021, l1013 -> operator-(), l1008, 
                                                l1007, l1006, l1005, l1009, l1020 -> operator-(), 
                                                l1021 -> operator-() });
    
        LineLoop* ll1006 = fluid2->addLineLoop({ l1027 -> operator-(), l1026 -> operator-(), 
                                                l1025 -> operator-(), l1024 -> operator-(),
                                                l1023 -> operator-(), l1022 -> operator-(),
                                                l1019, l1017, l1016, l1018, l1015, l1014}); 
        

        //Transfinite lines 
        int t1 = 40; int t2 = 130; int t3 = 70; int t4 = 7;
        //corners
        fluid2 -> transfiniteLine({ l1001 }, t2, 1.02);
        fluid2 -> transfiniteLine({ l1002 -> operator-()}, t3, 1.05);
        fluid2 -> transfiniteLine({ l1003 }, t3, 1.05);
        fluid2 -> transfiniteLine({ l1004 -> operator-()}, t2, 1.02);
        fluid2 -> transfiniteLine({ l1005 }, t2, 1.02);
        fluid2 -> transfiniteLine({ l1006 -> operator-()}, t3, 1.04);
        fluid2 -> transfiniteLine({ l1007 }, t3, 1.04);
        fluid2 -> transfiniteLine({ l1008 -> operator-()}, t2, 1.02);
        fluid2 -> transfiniteLine({ l1009 }, t4, 1.25);
        fluid2 -> transfiniteLine({ l1010 }, t4, 1.25);
        fluid2 -> transfiniteLine({ l1011 }, t4, 1.25);
        fluid2 -> transfiniteLine({ l1012 }, t4, 1.05);

        fluid2 -> transfiniteLine({ l1013 }, t4, 1.25);
        // fluid2 -> transfiniteLine({ l1014 }, r1);
        // fluid2 -> transfiniteLine({ l1015 }, r1);
        // fluid2 -> transfiniteLine({ l1016 }, r1);
        // //2r
        // fluid2 -> transfiniteLine({ l1017 }, r2);
        // fluid2 -> transfiniteLine({ l1018 }, r2);
        // fluid2 -> transfiniteLine({ l1019 }, r2);
        fluid2 -> transfiniteLine({ l1020 }, 50,1.03);
        fluid2 -> transfiniteLine({ l1021 }, 30,1.03);

        PlaneSurface* s1001 = fluid2 -> addPlaneSurface({ll1001});
        PlaneSurface* s1002 = fluid2 -> addPlaneSurface({ll1002});
        PlaneSurface* s1003 = fluid2 -> addPlaneSurface({ll1003});
        PlaneSurface* s1004 = fluid2 -> addPlaneSurface({ll1004});
        PlaneSurface* s1005 = fluid2 -> addPlaneSurface({ll1005});
        PlaneSurface* s1006 = fluid2 -> addPlaneSurface({ll1006});
    

        fluid2 -> addBoundaryCondition("NEUMANN", l1005, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l1006, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l1007, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l1008, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l1009, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l1010, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l1011, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l1012, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l1013, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l1020, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l1021, {}, {}, {}, "GLOBAL");


        fluid2 -> addBoundaryCondition("NEUMANN", l1014, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l1015, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l1016, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l1017, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l1018, {}, {}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("NEUMANN", l1019, {}, {}, {}, "GLOBAL");

        fluid2 -> addBoundaryCondition("GLUE", l1022, {0}, {0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l1023, {0}, {0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l1024, {0}, {0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l1025, {0}, {0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l1026, {0}, {0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("GLUE", l1027, {0}, {0}, {}, "GLOBAL");

        fluid2 -> addBoundaryCondition("DIRICHLET", l1001, {0.0}, {0.0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l1002, {0}, {0.0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l1003, {0}, {0}, {}, "GLOBAL");
        fluid2 -> addBoundaryCondition("DIRICHLET", l1004, {0}, {0}, {}, "GLOBAL");

        //Surfaces
        fluid2 -> transfiniteSurface({s1001}, "Right", {});
        fluid2 -> transfiniteSurface({s1002}, "Right", {});
        fluid2 -> transfiniteSurface({s1003}, "Right", {});
        fluid2 -> transfiniteSurface({s1004}, "Right", {});
        // fluid2 -> transfiniteSurface({s105}, "Alternated", {p106,p110,p111,p107});
        // fluid2 -> transfiniteSurface({s106}, "Alternated", {p107,p111,p112,p108});
        // fluid2 -> transfiniteSurface({s107}, "Alternated", {p108,p112,p113,p109});
        // fluid2 -> transfiniteSurface({s108}, "Alternated", {p109,p113,p110,p106});
        
        
        //Surfaces
        // fluid2->transfiniteSurface({ s1 }, "Left", {p10,p11,p16,p15});
        // fluid2->transfiniteSurface({ s2 }, "Left", {p11,p12,p17,p16});
        
    if (rank == 0){
  
        FluidDomain* problem = new FluidDomain(fluid1);
        // problem -> addSurfaceMaterial({ s20,s21,s22,s23 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
        problem -> generateMesh(T3, DELAUNAY, "coarse", "", true, true);

        FluidDomain* problem2 = new FluidDomain(fluid2);
        // problem2 -> addSurfaceMaterial({ s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12 }, 1.0, 1.0, 1.0, "PLANE_STRESS");
        problem2 -> generateMesh(T3, DELAUNAY, "fine", "", false, true);

        //problem -> readInput("exemplo.msh",0);
	};

//==========================================================================
//==============================PROBLEM MESHES==============================
//==========================================================================

	MPI_Barrier(PETSC_COMM_WORLD);
    
    
    coarseModel.meshReading(fluid1,"../airfoil_data.txt","coarse.msh","mirror.txt",0);
    fineModel.meshReading(fluid2,"../airfoil_data.txt","fine.msh","mirror_fine.txt",0);
   // } 
	MPI_Barrier(PETSC_COMM_WORLD);

    coarseModel.getFluidParameters().setSolver(SolverType::EMumps);
    fineModel.getFluidParameters().setSolver(SolverType::EMumps);


    coarseModel.setProblemType(ProblemType::ENavierStokes);
    fineModel.setProblemType(ProblemType::ENavierStokes);

    arlequinProblem.setArlequinStabilization(ArlequinStabType::EOption2);
    arlequinProblem.setFluidModels(coarseModel, fineModel) ; 


    // char in_solid[32] = "cantilever.txt";

    // coupledProblem.setArlequinAndSolidModels(arlequinProblem,in_solid);
    // coupledProblem.solveFSIProblemArlequin(100000);

    bool noPitch=true;
    if (noPitch) {
        for (int i = 0; i < arlequinProblem.numNodesFine; i++){
            VecDouble xn(2);
            VecDouble xi = arlequinProblem.fineModel.nodes_[i] -> getInitialCoordinates();       
            VecDouble x = arlequinProblem.fineModel.nodes_[i] -> getCoordinates();       

            double a = -20 * M_PI / 180 + 10 * M_PI / 180 * std::cos(2.*M_PI*0);// + 10 * M_PI / 180;

            // std::cout << " AAA " << a << std::endl;

            xn[0] = 0.5 + (xi[0]-0.5) * std::cos(a) - (xi[1]-0.0) * std::sin(a);
            xn[1] = 0.0 + (xi[0]-0.5) * std::sin(a) + (xi[1]-0.0) * std::cos(a);

            arlequinProblem.fineModel.nodes_[i] -> setPreviousCoordinates(0,x[0]);
            arlequinProblem.fineModel.nodes_[i] -> setPreviousCoordinates(1,x[1]);

            arlequinProblem.fineModel.nodes_[i] -> setCoordinates(xn);
        };


        arlequinProblem.solveArlequinProblem(3, 1.e-7, 2, 1);
    } else {
        arlequinProblem.solveArlequinProblemMoving(3, 1.e-7, 2, 1);
    }
        //