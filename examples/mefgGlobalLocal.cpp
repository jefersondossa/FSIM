#include "GeoMesh.h"
#include "LinearAnalysis.h"
#include "MixedCompMesh.h"
#include "GmshTools.h"
#include "VTUGenerator.h"
#include "L2Projection.h"
#include "Elasticity2D.h"
#include "CompMeshTools.h"
#include "InterpolatedBC.h"
#include "EigenSpMatrix.h"
#include "ElementEnriched.h"
#include "Element.h"
#include "GlobalLocalEnrichment.h"
#include "MEFGGlobalLocalTools.h"

//Comentar erro ao debugar após incluir a weakform GlobalLocalEnrichment e o element ElementEnriched. 

// Defines the problem dimension
const int dimension = 2;
auto forcingFunction = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    
    force[0] = -1*y;
    //force[1] = -1;
};
auto forcingFunctionAB = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    
    REAL A = 1.0;
    REAL Q = 0.543075579;
    REAL lambda = 0.544483737;
    REAL r = sqrt(x*x + y*y);
    REAL theta = atan2(y, x);

    MatrixDouble stress(2,2);
    //stress x
    stress(0,0) = A*lambda*pow(r, lambda-1)*((2 - Q*(lambda+1))*cos((lambda-1)*theta) - (lambda-1)*cos((lambda-3)*theta));
    //stress y
    stress(1,1) = A*lambda*pow(r, lambda-1)*((2 + Q*(lambda+1))*cos((lambda-1)*theta) + (lambda-1)*cos((lambda-3)*theta));
    //stress xy
    stress(0,1) = A*lambda*pow(r, lambda-1)*((lambda-1)*sin((lambda-3)*theta) + Q*(lambda+1)*sin((lambda-1)*theta));
    stress(1,0) = stress(0,1);

    VecDouble n(2);
    n[0] = sqrt(2)/2;
    n[1] = sqrt(2)/2;

    force = stress * n;
};

auto forcingFunctionBC = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    
    REAL A = 1.0;
    REAL Q = 0.543075579;
    REAL lambda = 0.544483737;
    REAL r = sqrt(x*x + y*y);
    REAL theta = atan2(y, x);

    MatrixDouble stress(2,2);
    //stress x
    stress(0,0) = A*lambda*pow(r, lambda-1)*((2 - Q*(lambda+1))*cos((lambda-1)*theta) - (lambda-1)*cos((lambda-3)*theta));
    //stress y
    stress(1,1) = A*lambda*pow(r, lambda-1)*((2 + Q*(lambda+1))*cos((lambda-1)*theta) + (lambda-1)*cos((lambda-3)*theta));
    //stress xy
    stress(0,1) = A*lambda*pow(r, lambda-1)*((lambda-1)*sin((lambda-3)*theta) + Q*(lambda+1)*sin((lambda-1)*theta));
    stress(1,0) = stress(0,1);

    VecDouble n(2);
    n[0] = - sqrt(2)/2;
    n[1] = sqrt(2)/2;

    force = stress * n;
};

auto forcingFunctionEF = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    
    REAL A = 1.0;
    REAL Q = 0.543075579;
    REAL lambda = 0.544483737;
    REAL r = sqrt(x*x + y*y);
    REAL theta = atan2(y, x);

    MatrixDouble stress(2,2);
    //stress x
    stress(0,0) = A*lambda*pow(r, lambda-1)*((2 - Q*(lambda+1))*cos((lambda-1)*theta) - (lambda-1)*cos((lambda-3)*theta));
    //stress y
    stress(1,1) = A*lambda*pow(r, lambda-1)*((2 + Q*(lambda+1))*cos((lambda-1)*theta) + (lambda-1)*cos((lambda-3)*theta));
    //stress xy
    stress(0,1) = A*lambda*pow(r, lambda-1)*((lambda-1)*sin((lambda-3)*theta) + Q*(lambda+1)*sin((lambda-1)*theta));
    stress(1,0) = stress(0,1);

    VecDouble n(2);
    n[0] = -sqrt(2)/2;
    n[1] = -sqrt(2)/2;

    force = stress * n;
};

auto forcingFunctionFA = [](const VecDouble &coord, VecDouble &force){
    const auto &x=coord[0];
    const auto &y=coord[1];
    
    REAL A = 1.0;
    REAL Q = 0.543075579;
    REAL lambda = 0.544483737;
    REAL r = sqrt(x*x + y*y);
    REAL theta = atan2(y, x);

    MatrixDouble stress(2,2);
    //stress x
    stress(0,0) = A*lambda*pow(r, lambda-1)*((2 - Q*(lambda+1))*cos((lambda-1)*theta) - (lambda-1)*cos((lambda-3)*theta));
    //stress y
    stress(1,1) = A*lambda*pow(r, lambda-1)*((2 + Q*(lambda+1))*cos((lambda-1)*theta) + (lambda-1)*cos((lambda-3)*theta));
    //stress xy
    stress(0,1) = A*lambda*pow(r, lambda-1)*((lambda-1)*sin((lambda-3)*theta) + Q*(lambda+1)*sin((lambda-1)*theta));
    stress(1,0) = stress(0,1);

    VecDouble n(2);
    n[0] = sqrt(2)/2;
    n[1] = -sqrt(2)/2;

    force = stress * n;
};
auto exactSol = [](const VecDouble &coord, VecDouble &u, MatrixDouble &gradU){
    const auto &x=coord[0];
    const auto &y=coord[1];

    REAL A = 1.0;
    REAL Q = 0.543075579;
    REAL lambda = 0.544483737;
    REAL r = sqrt(x*x + y*y);
    REAL theta = atan2(y, x);

    REAL young = 1.0;
    REAL nu = 0.3;
    REAL kappa = 3 - 4*nu;
    REAL G = young / (2*(1+nu));

    u[0] = (A/(2*G)) * pow(r, lambda) * ((kappa - Q * (lambda+1))*cos(lambda*theta) - lambda*cos((lambda-2)*theta));
    u[1] = (A/(2*G)) * pow(r, lambda) * ((kappa + Q * (lambda+1))*sin(lambda*theta) + lambda*sin((lambda-2)*theta));
    
    REAL dUxdr = (A/(2*G)) * lambda * pow(r, (lambda-1)) * ((kappa - Q*(lambda+1))*cos(lambda*theta) - lambda*cos((lambda-2)*theta));
    REAL dUydr = (A/(2*G)) * lambda * pow(r, (lambda-1)) * ((kappa + Q*(lambda+1))*cos(lambda*theta) + lambda*cos((lambda-2)*theta));

    REAL dUxdtheta = (A/(2*G)) * pow(r, lambda) * (-(kappa - Q*(lambda+1)) * lambda * sin(lambda*theta) + lambda * (lambda-2) * sin((lambda-2)*theta));
    REAL dUydtheta = (A/(2*G)) * pow(r, lambda) * ((kappa + Q*(lambda+1)) * lambda * cos(lambda*theta) + lambda * (lambda-2) * cos((lambda-2)*theta));

    gradU(0,0) = dUxdr*(x/r) + dUxdtheta*(-y/pow(r, 2));
    gradU(0,1) = dUxdr*(y/r) + dUxdtheta*(x/pow(r, 2));
    gradU(1,0) = dUydr*(x/r) + dUydtheta*(-y/pow(r, 2));
    gradU(1,1) = dUydr*(y/r) + dUxdtheta*(x/pow(r, 2));
};

void CreateGlobalModel(CompMesh *cmeshG);
void SolveGlobalProblem(CompMesh *cmeshG);
void CreateLocalModel(CompMesh *cmeshG, CompMesh *cmeshL);
void SolveLocalProblem(CompMesh *cmeshL);
void LocalToGlobalCorrespondence(CompMesh *cmeshG, CompMesh *cmeshL);
void LocalToGlobalCorrespondenceBoundary(CompMesh *cmeshG, CompMesh *cmeshL);
void CreateEnrichedModel(CompMesh *cmeshG, CompMesh *cmeshL);
void SolveEnrichedProblem(CompMesh *cmeshG);
bool CheckConvergence(int it);

int overlappingRegion;
int overlappingNHDirichletBoundary;
int overlappingNHNeumannBoundary;
int globalLocalIterations;
REAL globalLocalTolerance;
REAL ModElasticity = 1.;
REAL PoissonRatio = 0.; //0.3;//0.49999;

VecDouble previousSol;
VecDouble currentSol;
//Local index to global index correspondence for elements
std::map<int,int> globalElementCorrespondence;
//Local node/integration point to global node/integration point correspondence 
std::map<int, MatrixDouble> globalNodeCorrespondence;
//Global connects to be enriched to the index of the new connect associated with the enriched degree of freedom
std::map<int,int> enrichedConnects;

int main(int argc, char **args) { 

    overlappingRegion = 2;
    globalLocalIterations = 2;
    globalLocalTolerance = 1e-4;
    overlappingNHDirichletBoundary = 2;
    overlappingNHNeumannBoundary = 4;

    //GLOBAL MODEL
    std::cout << "Solve Global Problem \n";

    //Create Global Model
    GeoMesh *gmeshG = new GeoMesh();
    GmshTools::Read(*gmeshG,"../global.msh");
    CompMesh *cmeshG = new CompMesh(gmeshG,ApproxType::EIsoparametric);
    CreateGlobalModel(cmeshG);
    gmeshG->Print("gmeshGlobal.txt");
    //cmeshG->Print("cmeshGlobal.txt");
    
    //Solve Global Problem
    SolveGlobalProblem(cmeshG);

    CompMesh *cmeshL = nullptr;
    GeoMesh * gmeshL = new GeoMesh();

    GlobalLocalEnrichment *globalLocal = new GlobalLocalEnrichment(overlappingRegion, dimension, ModElasticity, PoissonRatio, false, 1.0);
    globalLocal->SetForcingFunction(forcingFunction);

    int it = 0;
    //for(int it = 0; it < globalLocalIterations; it++){
    while(it < globalLocalIterations){
        std::cout << "\nGlobal-Local Iteration " << it + 1 << ": ";

        //LOCAL MODEL
        std::cout << "\nSolve Local Problem \n";
        if(it == 0){
            //Create Local Model
            GmshTools::Read(*gmeshL,"../local.msh");
            cmeshL = new CompMesh(gmeshL,ApproxType::EIsoparametric);
            CreateLocalModel(cmeshG, cmeshL);
            gmeshL->Print("gmeshLocal.txt");
            //cmeshL->Print("cmeshLocal.txt");

            MEFGGlobalLocalTools::LocalToGlobalCorrespondence(cmeshG,cmeshL, globalElementCorrespondence, globalNodeCorrespondence, overlappingRegion, overlappingNHDirichletBoundary);
            MEFGGlobalLocalTools::LocalToGlobalCorrespondenceBoundary(cmeshG,cmeshL, globalElementCorrespondence, globalNodeCorrespondence, overlappingNHNeumannBoundary);
        }
        //Solve Local Problem
        SolveLocalProblem(cmeshL);

        //ENRICHED GLOBAL MODEL
        std::cout << "\nSolve Enriched Global Problem \n";
        if(it == 0) MEFGGlobalLocalTools::CreateEnrichedModel(cmeshG,cmeshL, globalElementCorrespondence, globalNodeCorrespondence, enrichedConnects, overlappingNHNeumannBoundary, globalLocal);

        cmeshG->Print("cmeshEnrichedGlobal.txt");
        SolveEnrichedProblem(cmeshG);

        if(CheckConvergence(it)) break;
        it++;
   };
}   

void CreateGlobalModel(CompMesh *cmeshG){

    Elasticity2D *matelasticityG1 = new Elasticity2D(1, ModElasticity, PoissonRatio, false, 1.0); //domínio global
    Elasticity2D *matelasticityG2 = new Elasticity2D(2, ModElasticity, PoissonRatio, false, 1.0); //domínio local

    //matelasticityG1->SetExactSolution(exactSol);

    cmeshG->InsertMaterial(matelasticityG1);
    cmeshG->InsertMaterial(matelasticityG2);

    //Chapa retangular tracionada
    enrichedConnects[0]=-1;
    enrichedConnects[5]=-1;

    //Chapa L 
    //enrichedConnects[8]=-1;

    //Q4 element
    //enrichedConnects[28]=-1;
    //enrichedConnects[29]=-1;
    //enrichedConnects[51]=-1;
    //enrichedConnects[52]=-1;
    //enrichedConnects[62]=-1;
    //enrichedConnects[63]=-1;
    //enrichedConnects[64]=-1;

    //Q9 element
    /*enrichedConnects[42]=-1;
    enrichedConnects[114]=-1;
    enrichedConnects[216]=-1;
    enrichedConnects[44]=-1;
    enrichedConnects[116]=-1;
    enrichedConnects[222]=-1;
    enrichedConnects[217]=-1;
    enrichedConnects[223]=-1;
    enrichedConnects[111]=-1;
    enrichedConnects[207]=-1;
    enrichedConnects[209]=-1;
    enrichedConnects[213]=-1;
    enrichedConnects[113]=-1;
    enrichedConnects[214]=-1;
    enrichedConnects[45]=-1;
    enrichedConnects[198]=-1;
    enrichedConnects[200]=-1;
    enrichedConnects[202]=-1;
    enrichedConnects[46]=-1;
    enrichedConnects[203]=-1;*/

    //Q4 element
    /*enrichedConnects[44]=-1;
    enrichedConnects[45]=-1;
    enrichedConnects[113]=-1;
    enrichedConnects[116]=-1;
    enrichedConnects[206]=-1;
    enrichedConnects[215]=-1;
    enrichedConnects[221]=-1;*/

    /*enrichedConnects[43]=-1;
    enrichedConnects[44]=-1;
    enrichedConnects[45]=-1;
    enrichedConnects[46]=-1;
    enrichedConnects[112]=-1;
    enrichedConnects[113]=-1;
    enrichedConnects[115]=-1;
    enrichedConnects[116]=-1;
    enrichedConnects[198]=-1;
    enrichedConnects[200]=-1;
    enrichedConnects[204]=-1;
    enrichedConnects[206]=-1;
    enrichedConnects[207]=-1;
    enrichedConnects[209]=-1;
    enrichedConnects[213]=-1;
    enrichedConnects[215]=-1;
    enrichedConnects[216]=-1;
    enrichedConnects[217]=-1;
    enrichedConnects[220]=-1;
    enrichedConnects[221]=-1;*/

    //Q9 element
    /*enrichedConnects[691]=-1;
    enrichedConnects[708]=-1;
    enrichedConnects[77]=-1;
    enrichedConnects[714]=-1;
    enrichedConnects[716]=-1;
    enrichedConnects[80]=-1;
    enrichedConnects[233]=-1;
    enrichedConnects[237]=-1;
    enrichedConnects[76]=-1;
    enrichedConnects[72]=-1;
    enrichedConnects[767]=-1;
    enrichedConnects[775]=-1;
    enrichedConnects[244]=-1;
    enrichedConnects[819]=-1;
    enrichedConnects[815]=-1;
    enrichedConnects[742]=-1;
    enrichedConnects[773]=-1;
    enrichedConnects[240]=-1;
    enrichedConnects[818]=-1;
    enrichedConnects[791]=-1;*/

    /*enrichedConnects[686]=-1;
    enrichedConnects[702]=-1;   
    enrichedConnects[688]=-1;
    enrichedConnects[705]=-1;
    enrichedConnects[78]=-1;
    enrichedConnects[710]=-1;
    enrichedConnects[712]=-1;
    enrichedConnects[707]=-1;
    enrichedConnects[709]=-1;
    enrichedConnects[81]=-1;
    enrichedConnects[690]=-1;
    enrichedConnects[711]=-1;
    enrichedConnects[691]=-1;
    enrichedConnects[708]=-1;
    enrichedConnects[77]=-1;
    enrichedConnects[713]=-1;
    enrichedConnects[715]=-1;
    enrichedConnects[714]=-1;
    enrichedConnects[716]=-1;
    enrichedConnects[80]=-1;
    enrichedConnects[232]=-1;
    enrichedConnects[236]=-1;
    enrichedConnects[233]=-1;
    enrichedConnects[237]=-1;
    enrichedConnects[8]=-1;
    enrichedConnects[76]=-1;
    enrichedConnects[72]=-1;
    enrichedConnects[75]=-1;
    enrichedConnects[71]=-1;
    enrichedConnects[754]=-1;
    enrichedConnects[768]=-1;
    enrichedConnects[767]=-1;
    enrichedConnects[775]=-1;
    enrichedConnects[244]=-1;
    enrichedConnects[819]=-1;
    enrichedConnects[815]=-1;
    enrichedConnects[817]=-1;
    enrichedConnects[796]=-1;
    enrichedConnects[737]=-1;
    enrichedConnects[766]=-1;
    enrichedConnects[742]=-1;
    enrichedConnects[773]=-1;
    enrichedConnects[240]=-1;
    enrichedConnects[818]=-1;
    enrichedConnects[791]=-1;
    enrichedConnects[816]=-1;
    enrichedConnects[785]=-1;
    enrichedConnects[751]=-1;
    enrichedConnects[771]=-1;
    enrichedConnects[770]=-1;
    enrichedConnects[774]=-1;
    enrichedConnects[243]=-1;
    enrichedConnects[822]=-1;
    enrichedConnects[821]=-1;
    enrichedConnects[824]=-1;
    enrichedConnects[799]=-1;
    enrichedConnects[735]=-1;
    enrichedConnects[769]=-1;
    enrichedConnects[741]=-1;
    enrichedConnects[772]=-1;
    enrichedConnects[239]=-1;
    enrichedConnects[820]=-1;
    enrichedConnects[790]=-1;
    enrichedConnects[823]=-1;
    enrichedConnects[784]=-1;*/

    //Q4 element
    /*enrichedConnects[76]=-1;
    enrichedConnects[77]=-1;
    enrichedConnects[237]=-1;
    enrichedConnects[244]=-1;
    enrichedConnects[734]=-1;
    enrichedConnects[783]=-1;
    enrichedConnects[809]=-1;*/

    /*enrichedConnects[73]=-1;
    enrichedConnects[74]=-1;
    enrichedConnects[75]=-1;
    enrichedConnects[76]=-1;
    enrichedConnects[77]=-1;
    enrichedConnects[78]=-1;
    enrichedConnects[79]=-1;
    enrichedConnects[80]=-1;
    enrichedConnects[234]=-1;
    enrichedConnects[235]=-1;
    enrichedConnects[236]=-1;
    enrichedConnects[237]=-1;
    enrichedConnects[241]=-1;
    enrichedConnects[242]=-1;
    enrichedConnects[243]=-1;
    enrichedConnects[244]=-1;
    enrichedConnects[686]=-1;
    enrichedConnects[688]=-1;
    enrichedConnects[692]=-1;
    enrichedConnects[694]=-1;
    enrichedConnects[702]=-1;
    enrichedConnects[705]=-1;
    enrichedConnects[719]=-1;
    enrichedConnects[725]=-1;
    enrichedConnects[727]=-1;
    enrichedConnects[728]=-1;
    enrichedConnects[729]=-1;
    enrichedConnects[730]=-1;
    enrichedConnects[731]=-1;
    enrichedConnects[732]=-1;
    enrichedConnects[733]=-1;
    enrichedConnects[734]=-1;
    enrichedConnects[735]=-1;
    enrichedConnects[737]=-1;
    enrichedConnects[741]=-1;
    enrichedConnects[743]=-1;
    enrichedConnects[751]=-1;
    enrichedConnects[754]=-1;
    enrichedConnects[768]=-1;
    enrichedConnects[774]=-1;
    enrichedConnects[776]=-1;
    enrichedConnects[777]=-1;
    enrichedConnects[778]=-1;
    enrichedConnects[779]=-1;
    enrichedConnects[780]=-1;
    enrichedConnects[781]=-1;
    enrichedConnects[782]=-1;
    enrichedConnects[783]=-1;
    enrichedConnects[784]=-1;
    enrichedConnects[785]=-1;
    enrichedConnects[788]=-1;
    enrichedConnects[789]=-1;
    enrichedConnects[796]=-1;
    enrichedConnects[799]=-1;
    enrichedConnects[805]=-1;
    enrichedConnects[806]=-1;
    enrichedConnects[807]=-1;
    enrichedConnects[808]=-1;
    enrichedConnects[809]=-1;
    enrichedConnects[810]=-1;
    enrichedConnects[811]=-1;
    enrichedConnects[812]=-1;
    enrichedConnects[813]=-1;
    enrichedConnects[814]=-1;*/

    //Q9 element
    /*enrichedConnects[2556]=-1;
    enrichedConnects[2660]=-1;
    enrichedConnects[141]=-1;
    enrichedConnects[2664]=-1;
    enrichedConnects[2666]=-1;
    enrichedConnects[148]=-1;
    enrichedConnects[477]=-1;
    enrichedConnects[485]=-1;
    enrichedConnects[140]=-1;
    enrichedConnects[132]=-1;
    enrichedConnects[2929]=-1;
    enrichedConnects[2935]=-1;
    enrichedConnects[500]=-1;
    enrichedConnects[3133]=-1;
    enrichedConnects[3129]=-1;
    enrichedConnects[2791]=-1;
    enrichedConnects[2933]=-1;
    enrichedConnects[492]=-1;
    enrichedConnects[3132]=-1;
    enrichedConnects[3011]=-1;*/

    /*enrichedConnects[2526]=-1;
    enrichedConnects[2607]=-1;
    enrichedConnects[2542]=-1;
    enrichedConnects[2610]=-1;
    enrichedConnects[2528]=-1;
    enrichedConnects[2618]=-1;
    enrichedConnects[2545]=-1;
    enrichedConnects[2621]=-1;
    enrichedConnects[144]=-1;
    enrichedConnects[2638]=-1;
    enrichedConnects[2640]=-1;
    enrichedConnects[2635]=-1;
    enrichedConnects[2637]=-1;
    enrichedConnects[2626]=-1;
    enrichedConnects[2628]=-1;
    enrichedConnects[2623]=-1;
    enrichedConnects[2625]=-1;
    enrichedConnects[151]=-1;
    enrichedConnects[2550]=-1;
    enrichedConnects[2639]=-1;
    enrichedConnects[2552]=-1;
    enrichedConnects[2636]=-1;
    enrichedConnects[2547]=-1;
    enrichedConnects[2627]=-1;
    enrichedConnects[2549]=-1;
    enrichedConnects[2624]=-1;
    enrichedConnects[143]=-1;
    enrichedConnects[2641]=-1;
    enrichedConnects[2644]=-1;
    enrichedConnects[2643]=-1;
    enrichedConnects[2646]=-1;
    enrichedConnects[2629]=-1;
    enrichedConnects[2632]=-1;
    enrichedConnects[2631]=-1;
    enrichedConnects[2634]=-1;
    enrichedConnects[150]=-1;
    enrichedConnects[2530]=-1;
    enrichedConnects[2642]=-1;
    enrichedConnects[2551]=-1;
    enrichedConnects[2645]=-1;
    enrichedConnects[2531]=-1;
    enrichedConnects[2630]=-1;
    enrichedConnects[2548]=-1;
    enrichedConnects[2663]=-1;
    enrichedConnects[142]=-1;
    enrichedConnects[2651]=-1;
    enrichedConnects[2653]=-1;
    enrichedConnects[2647]=-1;
    enrichedConnects[2650]=-1;
    enrichedConnects[2649]=-1;
    enrichedConnects[2663]=-1;
    enrichedConnects[2659]=-1;
    enrichedConnects[2661]=-1;
    enrichedConnects[149]=-1;
    enrichedConnects[2553]=-1;
    enrichedConnects[2652]=-1;
    enrichedConnects[2555]=-1;
    enrichedConnects[2648]=-1;
    enrichedConnects[2554]=-1;
    enrichedConnects[2662]=-1;
    enrichedConnects[2556]=-1;
    enrichedConnects[2660]=-1;
    enrichedConnects[141]=-1;
    enrichedConnects[2654]=-1;
    enrichedConnects[2656]=-1;
    enrichedConnects[2655]=-1;
    enrichedConnects[2658]=-1;
    enrichedConnects[2657]=-1;
    enrichedConnects[2665]=-1;
    enrichedConnects[2664]=-1;
    enrichedConnects[2666]=-1;
    enrichedConnects[148]=-1;
    enrichedConnects[474]=-1;
    enrichedConnects[482]=-1;
    enrichedConnects[475]=-1;
    enrichedConnects[483]=-1;
    enrichedConnects[476]=-1;
    enrichedConnects[484]=-1;
    enrichedConnects[477]=-1;
    enrichedConnects[485]=-1;
    enrichedConnects[8]=-1;
    enrichedConnects[140]=-1;
    enrichedConnects[132]=-1;
    enrichedConnects[139]=-1;
    enrichedConnects[131]=-1;
    enrichedConnects[138]=-1;
    enrichedConnects[130]=-1;
    enrichedConnects[137]=-1;
    enrichedConnects[129]=-1;
    enrichedConnects[2846]=-1;
    enrichedConnects[2894]=-1;
    enrichedConnects[2893]=-1;
    enrichedConnects[2903]=-1;
    enrichedConnects[2902]=-1;
    enrichedConnects[2930]=-1;
    enrichedConnects[2929]=-1;
    enrichedConnects[2935]=-1;
    enrichedConnects[500]=-1;
    enrichedConnects[3133]=-1;
    enrichedConnects[3129]=-1;
    enrichedConnects[3131]=-1;
    enrichedConnects[3120]=-1;
    enrichedConnects[3122]=-1;
    enrichedConnects[3117]=-1;
    enrichedConnects[3119]=-1;
    enrichedConnects[3040]=-1;
    enrichedConnects[2770]=-1;
    enrichedConnects[2892]=-1;
    enrichedConnects[2784]=-1;
    enrichedConnects[2900]=-1;
    enrichedConnects[2783]=-1;
    enrichedConnects[2928]=-1;
    enrichedConnects[2791]=-1;
    enrichedConnects[2933]=-1;
    enrichedConnects[492]=-1;
    enrichedConnects[3132]=-1;
    enrichedConnects[3011]=-1;
    enrichedConnects[3130]=-1;
    enrichedConnects[3007]=-1;
    enrichedConnects[3121]=-1;
    enrichedConnects[3009]=-1;
    enrichedConnects[3118]=-1;
    enrichedConnects[2988]=-1;
    enrichedConnects[2843]=-1;
    enrichedConnects[2897]=-1;
    enrichedConnects[2896]=-1;
    enrichedConnects[2901]=-1;
    enrichedConnects[2899]=-1;
    enrichedConnects[2932]=-1;
    enrichedConnects[2931]=-1;
    enrichedConnects[2934]=-1;
    enrichedConnects[499]=-1;
    enrichedConnects[3136]=-1;
    enrichedConnects[3135]=-1;
    enrichedConnects[3138]=-1;
    enrichedConnects[3123]=-1;
    enrichedConnects[3126]=-1;
    enrichedConnects[3125]=-1;
    enrichedConnects[3128]=-1;
    enrichedConnects[3043]=-1;
    enrichedConnects[2753]=-1;
    enrichedConnects[2895]=-1;
    enrichedConnects[2782]=-1;
    enrichedConnects[2898]=-1;
    enrichedConnects[2758]=-1;
    enrichedConnects[2918]=-1;
    enrichedConnects[2789]=-1;
    enrichedConnects[2926]=-1;
    enrichedConnects[491]=-1;
    enrichedConnects[3134]=-1;
    enrichedConnects[3010]=-1;
    enrichedConnects[3137]=-1;
    enrichedConnects[2983]=-1;
    enrichedConnects[3124]=-1;
    enrichedConnects[3008]=-1;
    enrichedConnects[3127]=-1;
    enrichedConnects[2977]=-1;
    enrichedConnects[2835]=-1;
    enrichedConnects[2906]=-1;
    enrichedConnects[2905]=-1;
    enrichedConnects[2915]=-1;
    enrichedConnects[2914]=-1;
    enrichedConnects[2919]=-1;
    enrichedConnects[2917]=-1;
    enrichedConnects[2927]=-1;
    enrichedConnects[498]=-1;
    enrichedConnects[3142]=-1;
    enrichedConnects[3139]=-1;
    enrichedConnects[3142]=-1;
    enrichedConnects[3141]=-1;
    enrichedConnects[3155]=-1;
    enrichedConnects[3151]=-1;
    enrichedConnects[3153]=-1;
    enrichedConnects[3053]=-1;
    enrichedConnects[2767]=-1;
    enrichedConnects[2904]=-1;
    enrichedConnects[2787]=-1;
    enrichedConnects[2912]=-1;
    enrichedConnects[2786]=-1;
    enrichedConnects[2916]=-1;
    enrichedConnects[2790]=-1;
    enrichedConnects[2924]=-1;
    enrichedConnects[490]=-1;
    enrichedConnects[3143]=-1;
    enrichedConnects[3014]=-1;
    enrichedConnects[3140]=-1;
    enrichedConnects[3013]=-1;
    enrichedConnects[3154]=-1;
    enrichedConnects[3016]=-1;
    enrichedConnects[3152]=-1;
    enrichedConnects[2991]=-1;
    enrichedConnects[2832]=-1;
    enrichedConnects[2909]=-1;
    enrichedConnects[2908]=-1;
    enrichedConnects[2913]=-1;
    enrichedConnects[2911]=-1;
    enrichedConnects[2922]=-1;
    enrichedConnects[2921]=-1;
    enrichedConnects[2925]=-1;
    enrichedConnects[497]=-1;
    enrichedConnects[3147]=-1;
    enrichedConnects[3146]=-1;
    enrichedConnects[3150]=-1;
    enrichedConnects[3149]=-1;
    enrichedConnects[3158]=-1;
    enrichedConnects[3157]=-1;
    enrichedConnects[3160]=-1;
    enrichedConnects[3056]=-1;
    enrichedConnects[2751]=-1;
    enrichedConnects[2907]=-1;
    enrichedConnects[2785]=-1;
    enrichedConnects[2910]=-1;
    enrichedConnects[2757]=-1;
    enrichedConnects[2920]=-1;
    enrichedConnects[2788]=-1;
    enrichedConnects[2923]=-1;
    enrichedConnects[489]=-1;
    enrichedConnects[3145]=-1;
    enrichedConnects[3012]=-1;
    enrichedConnects[3148]=-1;
    enrichedConnects[2982]=-1;
    enrichedConnects[3156]=-1;
    enrichedConnects[3015]=-1;
    enrichedConnects[3159]=-1;
    enrichedConnects[2976]=-1;*/
    
    //BC 
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();

    //chapa retangular
    L2Projection * matbcG1 = new L2Projection(3,dimension-1,BoundaryConditionType::kDirichlet,val1,val2);
    val2[1] = 1.0;
    L2Projection * matbcG2 = new L2Projection(4,dimension-1,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    val2.setZero();
    L2Projection * matbcG3 = new L2Projection(5,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    matbcG3->SetForcingFunction(forcingFunction);

    //Chapa L teste: apoio fixo no lado BC e carregamento unitário uniforme ao longe de FA
    // L2Projection * matbcG1 = new L2Projection(3,dimension-1,BoundaryConditionType::kDirichlet,val1,val2);
    // L2Projection * matbcG2 = new L2Projection(4,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    // matbcG2->SetForcingFunction(forcingFunction);

    //Chapa L: apoio horizontal no pontos B e F, apoio vertical no ponto A
    // val2[0] = 1.0;
    // L2Projection * matbcG1 = new L2Projection(3,dimension-2,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    // val2.setZero();
    // val2[1] = 1.0;
    // L2Projection * matbcG2 = new L2Projection(4,dimension-2,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    
    //Chapa L: apoio fixo no ponto D, apoio vertical no ponto A
    //L2Projection * matbcG1 = new L2Projection(3,dimension-2,BoundaryConditionType::kDirichlet,val1,val2);
    //val2[0] = 1.0;
    //L2Projection * matbcG2 = new L2Projection(4,dimension-2,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);

    // val2.setZero();
    //L2Projection * matbcG3 = new L2Projection(5,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    // matbcG3->SetForcingFunction(forcingFunctionAB);
    // L2Projection * matbcG4 = new L2Projection(6,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    // matbcG4->SetForcingFunction(forcingFunctionBC);
    // L2Projection * matbcG5 = new L2Projection(7,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    // matbcG5->SetForcingFunction(forcingFunctionEF);
    // L2Projection * matbcG6 = new L2Projection(8,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    // matbcG6->SetForcingFunction(forcingFunctionFA);

    cmeshG->InsertMaterial(matbcG1);
    cmeshG->InsertMaterial(matbcG2);
    cmeshG->InsertMaterial(matbcG3);
    //cmeshG->InsertMaterial(matbcG4);
    //cmeshG->InsertMaterial(matbcG5);
    //cmeshG->InsertMaterial(matbcG6);
    
    cmeshG->AutoBuild();
}

void SolveGlobalProblem(CompMesh *cmeshG){

    LinearAnalysis anG(cmeshG,SolverType::ELDLt);
       
    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"SigmaX","SigmaY","TauXY"};
    VectorNames = {"Displacement"}; //, "ExactDisplacement"}

    anG.Run();
    //anG.PrintGlobalMatrix();
    anG.PrintSolution();
    //anG.PrintGlobalRhs();

    //VecDouble errors(4);
    //anG.PostProcessError(errors);

    EigenSpMatrix *spMat = dynamic_cast<EigenSpMatrix *>(anG.GlobalMatrix());
    if (!spMat) {
        std::cerr << "Error: GlobalMatrix is not of type EigenSpMatrix." << std::endl;
        return;
    }

    VecDouble sol = spMat->Solution();
    VecDouble rhs = spMat->Rhs();
    REAL strainEnergy = (sol.dot(rhs))/2;
    std::cout << std::fixed << std::setprecision(10) << "Strain Energy: "<< strainEnergy << std::endl;

    VTUGenerator::PrintResults(cmeshG,"globalResult",ScalarNames,VectorNames);
}

void CreateLocalModel(CompMesh *cmeshG, CompMesh * cmeshL){

    Elasticity2D* matelasticityL = new Elasticity2D(1, ModElasticity, PoissonRatio, false, 1.0);
    cmeshL->InsertMaterial(matelasticityL);

    //BC;
    MatrixDouble val1(2,2);
    val1.setZero();
    VecDouble val2(2);
    val2.setZero();
    
    //InterpolatedBC * matbcL1 = new InterpolatedBC(2,dimension-1,2,BoundaryConditionType::kDirichlet,&globalElementCorrespondence,&globalNodeCorrespondence,cmeshG);
    //Chapa L: apoio fixo no ponto D, apoio vertical no ponto A
    //L2Projection * matbcL2 = new L2Projection(3,dimension-2,BoundaryConditionType::kDirichlet,val1,val2);
    
    //chapa retangular
    val2[1] = 1.0;
    L2Projection * matbcL1 = new L2Projection(3,dimension-1,BoundaryConditionType::kDirectionalHomogeneousDirichlet,val1,val2);
    val2.setZero();
    InterpolatedBC * matbcL2 = new InterpolatedBC(2,dimension-1,2,BoundaryConditionType::kDirichlet,&globalElementCorrespondence,&globalNodeCorrespondence,cmeshG);
    //val2[0] = -2.;
    //L2Projection * matbcL2 = new L2Projection(2,dimension-1,BoundaryConditionType::kDirectionalNonHomogeneousDirichlet,val1,val2);
    L2Projection * matbcL3 = new L2Projection(4,dimension-1,BoundaryConditionType::kNeumann,val1,val2);
    matbcL3->SetForcingFunction(forcingFunction);

    cmeshL->InsertMaterial(matbcL1);
    cmeshL->InsertMaterial(matbcL2);
    cmeshL->InsertMaterial(matbcL3);
    
    cmeshL->AutoBuild();
}

void SolveLocalProblem(CompMesh *cmeshL){

    LinearAnalysis anL(cmeshL,SolverType::ELDLt);

    anL.Run();
    // anL.PrintGlobalMatrix();
    // anL.PrintGlobalRhs();
    anL.PrintSolution();

    std::vector<std::string> ScalarNames, VectorNames;
    ScalarNames = {"SigmaX","SigmaY","TauXY"};
    VectorNames = {"Displacement"};

    EigenSpMatrix *spMat = dynamic_cast<EigenSpMatrix *>(anL.GlobalMatrix());
    if (!spMat) {
        std::cerr << "Error: GlobalMatrix is not of type EigenSpMatrix." << std::endl;
        return;
    }

    VecDouble sol = spMat->Solution();
    VecDouble rhs = spMat->Rhs();

    REAL bignumber = 1.e20;
    for (int64_t i = 0; i < rhs.size(); i++){
        if (rhs[i]>1.e10) rhs[i] = 0.;
        if (rhs[i]<-1.e10) rhs[i] = 0.;
    }

    REAL strainEnergy = (sol.dot(rhs))/2;

    std::cout << "Strain Energy: "<< strainEnergy << std::endl;

    VTUGenerator::PrintResults(cmeshL,"localResult",ScalarNames,VectorNames); 
};

void SolveEnrichedProblem(CompMesh *cmeshG){

    LinearAnalysis anE(cmeshG,SolverType::ELU);

    anE.Run();
    // anE.PrintGlobalMatrix();
    // anE.PrintGlobalRhs();
    anE.PrintSolution();
    
    EigenSpMatrix *spMat = dynamic_cast<EigenSpMatrix *>(anE.GlobalMatrix());
    if (!spMat) {
        std::cerr << "Error: GlobalMatrix is not of type EigenSpMatrix." << std::endl;
        return;
    }

    currentSol = spMat->Solution();

    VecDouble rhs = spMat->Rhs();
    REAL strainEnergy = (currentSol.dot(rhs))/2;
    std::cout << "Strain Energy: "<< strainEnergy << std::endl;

    std::vector<std::string> ScalarNames, VectorNames;
    // ScalarNames = {"SigmaX","SigmaY","TauXY"};
    VectorNames = {"Displacement"};

    VTUGenerator::PrintResults(cmeshG,"Enriched",ScalarNames,VectorNames); 
};

bool CheckConvergence(int it){

    if(it == 0){
        previousSol = currentSol;
        return false;
    }else{
        MatrixDouble diff = currentSol - previousSol;

        if(diff.norm() <= globalLocalTolerance){
            return true;
        }else{
            previousSol = currentSol;
            return false;
        }
    }
};