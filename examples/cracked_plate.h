
   // Defines the problem dimension
    const int dimension = 2;
{
auto yieldFunction = [](const REAL &plast, REAL &sigma_y, REAL &hardening){
    // sigma_y =.5/sqrt(3.)+500.*plast;
    hardening = 0.;
    sigma_y =.24e9+hardening*plast;
};


auto exactSolElasticity2D = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord;
    auto pi = M_PI;
    REAL gPoisson = 0.;

    TVar2 Kolosov = 3. - 4.*gPoisson;
    TVar2 mu = gE / (2.*(1.+gPoisson));
    TVar2 theta = atan2(x[1],x[0]);
    
    auto thetaval = shapeFAD::val(theta);
    if (thetaval < (0.)) theta += 2. * M_PI;
    TVar2 r = sqrt(x[0]*x[0]+x[1]*x[1]);
    TVar2 a = 1.;
    TVar2 sinfty = 1.;
    TVar2 tinfty = 1.;
    TVar2 K1 = sinfty * sqrt(M_PI*a);//Modo de abertura
    TVar2 K2 = tinfty * sqrt(M_PI*a);//Modo de deslizamento


    TVar2 cosT2 = cos(0.5*theta);
    TVar2 sinT2 = sin(0.5*theta);

    TVar2 mult1 = sqrt(r) / (2. * mu * sqrt(2.*M_PI));
    TVar2 mult2 = Kolosov - 1. + 2.*sinT2*sinT2;
    TVar2 mult3 = Kolosov + 1. + 2.*cosT2*cosT2;
    TVar2 mult4 = Kolosov + 1. - 2.*cosT2*cosT2;
    TVar2 mult5 = Kolosov - 1. - 2.*sinT2*sinT2;
    
    disp[0] = K1*mult1*cosT2*mult2 + K2*mult1*sinT2*mult3;
    disp[1] = K1*mult1*sinT2*mult3 + K2*mult1*cosT2*mult4;




    u[0] = cos(pi*x)*sin(2.*pi*y);
    u[1] = cos(pi*y)*sin(pi*x);
    gradU(0,0) = -pi*sin(pi*x)*sin(2.*pi*y);
    gradU(0,1) = pi*cos(pi*x)*cos(pi*y);
    gradU(1,0) = 2.*pi*cos(pi*x)*cos(2.*pi*y);
    gradU(1,1) = -pi*sin(pi*x)*sin(pi*y);
};

auto forcingFunctionElasticity2D = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    REAL E=1.;
    REAL poisson=0.;
    auto pi = M_PI;
    force[0] = -(E*pi*pi*cos(pi*x)*((1.+poisson)*sin(pi*y) + 2.*(3. - 2.*poisson)*sin(2.*pi*y)))/(2.*(-1. + poisson*poisson));
    force[1] = -(E*pi*pi*(-((-3. + poisson)*cos(pi*y)) + 2.*(1. + poisson)*cos(2*pi*y))*sin(pi*x))/(2.*(-1. + poisson*poisson));
};

        



    CompMesh* cmesh = new CompMesh(); 

    Elasticity2D * matelas = new Elasticity2D(11,210.e9,.3,true,0.05);
    // ElasticityPositional2D * matelas = new ElasticityPositional2D(10,1.e3,.0,false);

    //BC
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();
    //Encastre
    val2[0] = 1.;
    L2Projection * matbc1 = new L2Projection(10,2,3,val1,val2);

    val2.setZero();
    //Encastre point
    L2Projection * matbc2 = new L2Projection(8,2,0,val1,val2);

    //Load
    val2[1] = 000.;
    L2Projection * matbc3 = new L2Projection(9,2,1,val1,val2);
    
    cmesh->InsertMaterial(matbc1);
    cmesh->InsertMaterial(matbc2);
    cmesh->InsertMaterial(matbc3);
    
    VonMises *plastmodel = new VonMises(matelas);
    plastmodel->SetUniaxialYieldFunction(yieldFunction);
    // cmesh->InsertMaterial(matelas);
    cmesh->InsertMaterial(plastmodel); 
    std::vector<L2Projection *> bcIncrement = {matbc3};
   
    GmshTools::Read(*cmesh,"../cantilever.msh");

    IncrementalAnalysis an(cmesh,SolverType::ELU, 100, bcIncrement,1.e-7,300);
    VecDouble increment(2);
    // increment.setZero();    
    increment[1] = -1000.;
    // increment[1] = 0.02;
    an.SetIncrement(increment);
#ifdef RELEASE_BUILD
    an.SType() = SolverType::EUmfpack;
    // an.SType() = SolverType::ECholmod;
#endif

    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"PlasticStrain"};
    VectorNames = {"Displacement","Stress","Strain","RealStress"};
    // VectorNames = {"Displacement","Stress","Strain"};
    an.Run("plasticitytest",ScalarNames,VectorNames);
    
    // std::vector<std::string> integrate = {"Solution","DerivativeX","DerivativeY"};
    // std::set<int> matid = {6};
    // std::map<std::string,VecDouble> result;
    // cmesh->Integrate(matid,integrate,result);
    // REAL aux = 206.9/((1.29)*(1.-2.*.29));
    // REAL valX = aux * ((1.-.29)*result["DerivativeX"][0] + .29*result["DerivativeY"][0])/.45;
    // REAL valY = aux * (.29*result["DerivativeX"][0] + (1.-.29)*result["DerivativeY"][0])/.45;
    // int a = 0;
    
}           