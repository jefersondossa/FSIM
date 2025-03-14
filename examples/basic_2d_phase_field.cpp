#include <fstream>
#include <TransientAnalysis.h>
#include <GmshTools.h>
#include <VTUGenerator.h>
#include <L2Projection.h>
#include <PhaseField.h>
#include <Elasticity2D.h>
#include <ElementTransient.h>
#include <ShapePoint.h>
#include <Node.h>
#include <memory>
#include <LagrangeMultiplier.h>
#include <ElementLagrangeMultiplier.h>
#include <ShapeQuadrilateralLin.h>
#include <ShapeTriangleLin.h>
#include <NullWeakForm.h>

void SetupBoundaryConditionsElasticity2D(CompMesh& modelElasticity2D)
{
    MatrixDouble val1(2, 2);
    VecDouble val2(2);
    constexpr auto kVolumeMatId = 15;
    auto *govEquationElasticity2D = new Elasticity2D(kVolumeMatId, 1e3, 0.3);
    modelElasticity2D.InsertMaterial(govEquationElasticity2D);

    val1.setZero();
    val2.setZero();
    val2[0] = 0.;

    constexpr auto kEngasteMatId = 16;

    auto *engasteBC = new L2Projection(kEngasteMatId, 2, BoundaryConditionType::kDirichlet, val1, val2);
    modelElasticity2D.InsertMaterial(engasteBC);

    val1.setZero();
    val2.setZero();
    val2[0] = 0.;

    constexpr auto kFreeMatTopId = 18;

    auto *freeBC = new L2Projection(kFreeMatTopId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelElasticity2D.InsertMaterial(freeBC);

    val1.setZero();
    val2.setZero();
    val2[0] = 0.;

    constexpr auto kFreeMatBottomId = 19;

    auto *freeBCEl2D = new L2Projection(kFreeMatBottomId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelElasticity2D.InsertMaterial(freeBCEl2D);

    val2.setZero();
    val2[1] = -10;
    constexpr auto kLoadMatId = 17;
    auto *El2D = new L2Projection(kLoadMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelElasticity2D.InsertMaterial(El2D);
    
    GmshTools::Read(modelElasticity2D, "../../rectangle.msh");

    // TODO: Disabled because now we recalculate the stiffness matrix contribution per node (based on phi)
    // Disables memory on elements
    // (makes sure elemental stiffness is not recalculated)
    // for(auto& [_, pWeakForm] : modelElasticity2D.MaterialVector()) {
        // govEquationElasticity2D->SetHasMemory(false);
    // }
}

void SetupBoundaryConditionsPhaseField(CompMesh& modelPhaseField)
{
    MatrixDouble val1(1, 1);
    VecDouble val2(1);
    val1.setZero();
    val2.setZero();
    constexpr auto kPhaseFieldInternalMatId = 15;
    auto *govEquationPF = new PhaseField(kPhaseFieldInternalMatId, 2, 1.e-4, 20.00, 1e-1);
    modelPhaseField.InsertMaterial(govEquationPF);

    constexpr auto kPhaseFieldLeftMatId = 16;
    auto *leftBCPF = new L2Projection(kPhaseFieldLeftMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelPhaseField.InsertMaterial(leftBCPF);

    constexpr auto kPhaseFieldRightMatId = 17;
    auto *rightBCPF = new L2Projection(kPhaseFieldRightMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelPhaseField.InsertMaterial(rightBCPF);

    constexpr auto kPhaseFieldTopMatId = 18;
    auto *topBCPF = new L2Projection(kPhaseFieldTopMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelPhaseField.InsertMaterial(topBCPF);

    // val2[0] = 1.0;
    constexpr auto kPhaseFieldBotttomMatId = 19;
    auto *bottomBCPF = new L2Projection(kPhaseFieldBotttomMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelPhaseField.InsertMaterial(bottomBCPF);
    
    GmshTools::Read(modelPhaseField, "../../rectangle.msh");

    //For the phase field implementation, we also need to impose an additional restrain to the final volume.
    //To do so, we create a new node, which serves to store the additional Lagrange multiplier variable.
    //However, the LM has the contribution of all elements. Thus, the new node need to be set as a new connectivity for all elements.
    //Note that this completely changes the data structure, and additional care is needed to handle it.
    //For instance, changes in the shape functions and assemble algorithm may be needed.
    //First, create a new node and a new element to store the LM.
    VecDouble Coor(3);
    Coor[0] = 0.0;
    Coor[1] = 0.0;
    Coor[2] = 0.0;
    Node *newNode = new Node(Coor,modelPhaseField.NodeVec().size(),1);
    newNode->AllocateTimeDerivatives();
    modelPhaseField.NodeVec().push_back(newNode);
    int64_t index = newNode->Index();
    VecInt connect(1);
    connect[0] = index;
    // Point
    Element* gel = nullptr;
    int matnull = modelPhaseField.GetNewMaterialId();
    NullWeakForm *nullwf = new NullWeakForm(matnull, 1);
    gel = new ElementTransient<ShapePoint>(index,connect,&modelPhaseField,nullwf);
    modelPhaseField.InsertElement(gel);
    
    // Creates the Lagrange multiplier material
    int matlagmult = modelPhaseField.GetNewMaterialId();
    LagrangeMultiplier *lagmult = new LagrangeMultiplier(matlagmult, 1);
    //Sets the initial/final volume constraint as a forcing function
    // const auto diam = 2*100.0/700;
    // double finalVol = 2.-12.*M_PI*diam*diam/4.0;
    static constexpr auto diam = 2*100.0/700;
    static constexpr auto holesVol = 12.*M_PI*diam*diam/4.0;
    static constexpr auto finalVol = 2.-holesVol;
    static constexpr auto finalVolRel = finalVol / 2.0;
    auto forcingFunction = [](const VecDouble &coord, VecDouble &force){
        force[0] = finalVolRel;
        // force[0] = 0.0;
    };
    lagmult->SetForcingFunction(forcingFunction);
    modelPhaseField.InsertMaterial(lagmult);
    //Now, create new elements of the type LagrangeMultiplier, which will be used to impose the additional constrain.
    // Here, the point element will always be the second and the volumetric element, the first.
    std::cout << "NElements before LagrangeMultiplier: " << modelPhaseField.NElements() << std::endl;
    for (auto &element : modelPhaseField.ElementVec())
    {
        if (element->Dimension() != modelPhaseField.Dimension()) continue;

        Element* gelmult = nullptr;
        auto *wf = modelPhaseField.Material(matlagmult);
        auto newindex = modelPhaseField.NElements();
        auto type = element->Type();

        switch (type)
        {
        case ElementType::EQuadrilateral:
            if (modelPhaseField.GetDefaultOrder()==1){
                gelmult = new ElementLagrangeMultiplier<ShapeQuadrilateralLin>(newindex, element, gel, &modelPhaseField, wf);
            } else {
                PanicButton();
            }
            break;
        case ElementType::ETriangle:
            if (modelPhaseField.GetDefaultOrder()==1){
                gelmult = new ElementLagrangeMultiplier<ShapeTriangleLin>(newindex, element, gel, &modelPhaseField, wf);
            } else {
                PanicButton();
            }
            break;
        
        default:
            PanicButton();
            break;
        }
        modelPhaseField.InsertElement(gelmult);
    }
    std::cout << "NElements after LagrangeMultiplier: " << modelPhaseField.NElements() << std::endl;
    delete [] modelPhaseField.part_elem;
    modelPhaseField.part_elem = new int[modelPhaseField.NElements()]();

}

int main()
{
    std::unique_ptr<CompMesh> modelElasticity2D = std::make_unique<CompMesh>();
    SetupBoundaryConditionsElasticity2D(*modelElasticity2D);
    LinearAnalysis anElasticity2D(modelElasticity2D.get(), SolverType::ELU);

    std::unique_ptr<CompMesh> modelPhaseField = std::make_unique<CompMesh>();
    SetupBoundaryConditionsPhaseField(*modelPhaseField);
    TransientAnalysis anPhaseField(modelPhaseField.get(), SolverType::ELU,false);
    // NonLinearAnalysis anPhaseField(modelPhaseField.get(), SolverType::ELU);
    anPhaseField.SetMaxIter(10);
    anPhaseField.SetTolerance(1e-6);

    const auto diam = 2*100.0/700;
    const auto espacamento_horizontal = 10.0/700;
    const auto espacamento_vertical = 70.0/700;
    
    const std::vector<VecDouble> hole_positions = {
        Eigen::Vector3d({diam + espacamento_horizontal, 0.0, 0.0}),
        Eigen::Vector3d({3.3*(diam + espacamento_horizontal), 0.0, 0.0}),
        Eigen::Vector3d({5.5*(diam + espacamento_horizontal), 0.0, 0.0}),

        Eigen::Vector3d({0, 1.0*diam, 0.0}),
        Eigen::Vector3d({0.0, 2.5*diam, 0.0}),

        Eigen::Vector3d({diam + espacamento_horizontal, espacamento_vertical + 1.5*diam, 0.0}),
        Eigen::Vector3d({3.3*(diam + espacamento_horizontal), espacamento_vertical + 1.5*diam, 0.0}),
        Eigen::Vector3d({5.5*(diam + espacamento_horizontal), espacamento_vertical + 1.5*diam, 0.0}),

        Eigen::Vector3d({2.0*(diam + espacamento_horizontal), 1.0*diam, 0.0}),
        Eigen::Vector3d({2.0*(diam + espacamento_horizontal), 2.5*diam, 0.0}),

        Eigen::Vector3d({diam + espacamento_horizontal, espacamento_vertical + 3*diam, 0.0}),
        Eigen::Vector3d({3.3*(diam + espacamento_horizontal), espacamento_vertical + 3*diam, 0.0}),
        Eigen::Vector3d({5.5*(diam + espacamento_horizontal), espacamento_vertical + 3*diam, 0.0}),

        Eigen::Vector3d({4.3*(diam + espacamento_horizontal), 1.0*diam, 0.0}),
        Eigen::Vector3d({4.3*(diam + espacamento_horizontal), 2.5*diam, 0.0}),

        Eigen::Vector3d({2.0, 1.0*diam, 0.0}),
        Eigen::Vector3d({2.0, 2.5*diam, 0.0}),
    };

    constexpr auto min_val = 1e-3;
    constexpr auto max_val = 1.0;

    for(size_t i = 0; i < anPhaseField.MeshVector()[0]->NNodes(); i++)
    {
        auto& node = anPhaseField.MeshVector()[0]->NodeVec()[i];

        const auto& node_pos = node->getCoordinates();

        bool has_void = false;
        for(const auto& hole : hole_positions)
        {
            if((node_pos - hole).norm() <= 0.8*diam/2.0)
            {
                has_void = true;
                break;
            }
        }

        auto& fSol = node->Solution();
        auto& fSolPrev = node->PrevSolution();
        fSol.resize(1);
        fSolPrev.resize(1);

        fSol[0] = min_val;//(has_void) ? 0.0 : max_val;
        fSolPrev[0] = fSol[0];
    }

    std::vector<std::string> ScalarNamesPhaseField, VectorNamesPhaseField;
    ScalarNamesPhaseField = {"Solution","TimeDerivative"};
    VectorNamesPhaseField = {"Derivative"};

    anPhaseField.PrintVariables("basic_2d_phase_field_dt", ScalarNamesPhaseField, VectorNamesPhaseField);

    if(modelElasticity2D->NElements() != modelPhaseField->NElements())
    {   
        // We assume that we are using the same mesh for phase field and elasticity!
        //Yes, but the check need to be different since now we have elements corresponding to the Lagrange multiplier
        // PanicButton();
    }

    anElasticity2D.Run();

    int i = 0;
    std::vector<std::string> ScalarNamesElasticity2D, VectorNamesElasticity2D;
    ScalarNamesElasticity2D = {"Compliance", "ComplianceSensibility","WeightFunction"};
    VectorNamesElasticity2D = {"Displacement", "Stress","Strain"};

    for (int64_t i_el = 0; i_el < modelElasticity2D->NElements(); i_el++)
    {
        auto elemPhaseField = modelPhaseField->ElementVec()[i_el];

        elemPhaseField->IntegrationData().fJ = 0;
    }

    VTUGenerator::PrintResults(modelElasticity2D.get(), "cantilever_2d_beam", ScalarNamesElasticity2D, VectorNamesElasticity2D, {}, 1);
    
    // Solving
    while(true)
    {
        // For the elasticity problem, transfer phase field solution to the weightFunction data structure.
        for (size_t inode = 0; inode < modelElasticity2D->NNodes(); inode++){
            modelElasticity2D->NodeVec()[inode]->setWeightFunction(std::max(min_val, pow(modelPhaseField->NodeVec()[inode]->Solution()[0], 3)));
        }
        // For each elasticity element, compute the compliance and then transfer the information
        // to the phase field corresponding element. The value of J need to be computed for each integration point.
        // Thus the //J data structure need to be a vector with lenght equal to the number of integration points.
        // Otherwhise, the value of J can also be stored in the weight function data structure.
        for (int64_t i_el = 0; i_el < modelElasticity2D->NElements() 
            && i_el < modelPhaseField->NElements(); i_el++)
        {
            auto elemElas2D = modelElasticity2D->ElementVec()[i_el];
            auto elemPhaseField = modelPhaseField->ElementVec()[i_el];
            
            if (elemElas2D->Dimension() != modelElasticity2D->Dimension())
            {
                continue;
            }
            elemElas2D->setIntegPointWeightFunction();
            
            //Trasfer complience to the phase field element
            const auto compliance_var_idx = elemElas2D->GetWeakForm()->VariableIndex("ComplianceSensibility");
            const auto nvar = elemElas2D->GetWeakForm()->NSolutionVariables(compliance_var_idx);
            VecDouble compl_vec(nvar);
            elemElas2D->Solution(compliance_var_idx, compl_vec);

            elemPhaseField->IntegrationData().fJ = compl_vec[0];//(compl_vec[0] > 0 ? 1.0 : -1.0);// ;
        }

        anElasticity2D.Run();
        VTUGenerator::PrintResults(modelElasticity2D.get(), "cantilever_2d_beam", ScalarNamesElasticity2D, VectorNamesElasticity2D, {}, i);

        anPhaseField.Run(1);
        VTUGenerator::PrintResults(modelPhaseField.get(), "phase_field_2d_", ScalarNamesPhaseField, VectorNamesPhaseField, {}, i);
        std::cout << "LAGRANGE MULTIPLIER = " << modelPhaseField->NodeVec().back()->Solution()[0] << std::endl;
        i++;
    }


    return 0;
}