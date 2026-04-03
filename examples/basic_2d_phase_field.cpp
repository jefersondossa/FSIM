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
#include <eso.hpp>

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <limits>

using namespace std;

const int Nx = (111) + 1, Ny = (37) + 1;
const double Lx = 3.0, Ly = 1.0;
const double dx = Lx / (Nx - 1);
const double dy = Ly / (Ny - 1);
const double dt = 4e-4;//5e-3;

constexpr double KK_ = 1.8;
const int PHASEFIELD_STEPS_PER_ELASTICITY = 1;
const int SAVE_EVERY = 8*PHASEFIELD_STEPS_PER_ELASTICITY * 50;
const int MAX_STEPS = 12000*SAVE_EVERY;

const double M = 1.0;
const double gamma_ = 0.01*8 * KK_;
const double ksi = 4e-4*0.6; //4e-4*1.4;

const double beta_ = 1.5 * KK_;

const double V = 0.4;

using Grid = vector<vector<double>>;
using GridStress = vector<vector<double>>;

void write_vtu(const Grid &phi, const std::string &filename) {
    const int Nx = phi.size();
    const int Ny = phi[0].size();
    const int numPoints = Nx * Ny;
    const int numCells = (Nx - 1) * (Ny - 1);

    std::ofstream out(filename);
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << numPoints << "\" NumberOfCells=\"" << numCells << "\">\n";

    // Point coordinates
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            out << i * dx << " " << j * dy << " 0.0\n";
        }
    }
    out << "        </DataArray>\n";
    out << "      </Points>\n";

    // Cell connectivity
    out << "      <Cells>\n";
    // Connectivity
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int j = 0; j < Ny - 1; ++j) {
        for (int i = 0; i < Nx - 1; ++i) {
            int p0 =  i      +  j      * Nx;
            int p1 = (i + 1) +  j      * Nx;
            int p2 = (i + 1) + (j + 1) * Nx;
            int p3 =  i      + (j + 1) * Nx;
            out << p0 << " " << p1 << " " << p2 << " " << p3 << "\n";
        }
    }
    out << "        </DataArray>\n";
    // Offsets
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int c = 1; c <= numCells; ++c) {
        out << 4 * c << "\n";
    }
    out << "        </DataArray>\n";
    // Cell types: 9 = VTK_QUAD
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int c = 0; c < numCells; ++c) {
        out << "9\n";
    }
    out << "        </DataArray>\n";
    out << "      </Cells>\n";

    // Point data: phase field
    out << "      <PointData Scalars=\"phi\">\n";
    out << "        <DataArray type=\"Float32\" Name=\"phi\" format=\"ascii\">\n";
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            out << phi[i][j] << "\n";
        }
    }
    out << "        </DataArray>\n";
    out << "      </PointData>\n";

    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";
}

double w_prime(double phi) {
    return 2.0 * phi * (1.0 - phi) * (1.0 - 2.0 * phi);
}

void initialize(Grid &phi) {
    constexpr auto A = 0.2;
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            double x = i * dx - 0.5;
            double y = j * dy - 0.5;
            phi[i][j] = V + A * sin(8 * M_PI * x) * sin(8 * M_PI * y);
        }
    }
}

double dC_dphi(double phi) {
    const double C0 = 1.0;   // stiffness of phase 0
    const double C1 = 10.0;  // stiffness of phase 1
    return C1 - C0;
}

double compute_lagrange_multiplier(const Grid &phi, const Grid &strain_energy, const Grid &mech_coupling_normalized) {
    double sum = 0.0;
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            double lap_phi = (phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1] - 4.0 * phi[i][j]) / (dx * dx);

            double wphi = w_prime(phi[i][j]);

            double mech = beta_ * mech_coupling_normalized[i][j];

            sum += ksi * lap_phi - wphi + (1.0 / gamma_) * mech;
        }
    }

    // Normalize by domain area
    double domain_area = (Nx - 2) * (Ny - 2);  // inner domain
    return sum / domain_area;
}

void calculate_mech_coupling(const Grid &phi, const Grid& strain_energy, Grid& mech_coupling_normalized) {
    double max_mech_coupling = std::numeric_limits<double>::min();
    for (int i = 1; i < Nx - 1; ++i) {
        for(int j = 1; j < Ny - 1; ++j) {
            const double mech_coupling = dC_dphi(phi[i][j]) * strain_energy[i][j];
            mech_coupling_normalized[i][j] = mech_coupling;
            max_mech_coupling = std::max(mech_coupling, max_mech_coupling);
        }
    }

    for (int i = 1; i < Nx - 1; ++i) {
        for(int j = 1; j < Ny - 1; ++j) {
            mech_coupling_normalized[i][j] /= max_mech_coupling;
            mech_coupling_normalized[i][j] = std::min(mech_coupling_normalized[i][j], 1.0);
        }
    }
}

void update(Grid &phi, Grid &phi_new, Grid& strain_energy, Grid& mech_coupling_normalized, int step) {
    calculate_mech_coupling(phi, strain_energy, mech_coupling_normalized);

    const double lambda = compute_lagrange_multiplier(phi, strain_energy, mech_coupling_normalized);

    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {

            const double x = i * dx;
            const double y = j * dy;

            const double lap_phi = (phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1] - 4.0 * phi[i][j]) / (dx * dx);
            const double wphi = w_prime(phi[i][j]);

            const double mech_coupling = beta_ * mech_coupling_normalized[i][j];

            const auto rhs = ksi * lap_phi - wphi + (1.0 / gamma_) * mech_coupling - lambda;
            phi_new[i][j] = phi[i][j] + dt * M * gamma_ * rhs;

            // phi_new[i][j] = std::min(1.0, std::max(0.0, phi_new[i][j]));
        }
    }

    // simple Neumann BC (zero gradient)
    for (int i = 0; i < Nx; ++i) {
        phi_new[i][0] = phi_new[i][1];
        phi_new[i][Ny-1] = phi_new[i][Ny-2];
    }
    for (int j = 0; j < Ny; ++j) {
        phi_new[0][j] = phi_new[1][j];
        phi_new[Nx-1][j] = phi_new[Nx-2][j];
    }
}

void SetupBoundaryConditionsElasticity2DCantileverRightBottom(CompMesh& modelElasticity2D)
{
    MatrixDouble val1(2, 2);
    VecDouble val2(2);
    constexpr auto kVolumeMatId = 15;
    auto *govEquationElasticity2D = new Elasticity2D(kVolumeMatId, 1e3, 0.3);
    modelElasticity2D.InsertMaterial(govEquationElasticity2D);
    val1.setZero();
    val2.setZero();
    constexpr auto kEngasteMatId = 18;
    auto *engasteBC = new L2Projection(kEngasteMatId, 2, BoundaryConditionType::kDirichlet, val1, val2);
    modelElasticity2D.InsertMaterial(engasteBC);
    val1.setZero();
    val2.setZero();
    val2[1] = -10;
    constexpr auto kLoadMatId = 17;
    auto *El2D2 = new L2Projection(kLoadMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelElasticity2D.InsertMaterial(El2D2);
    GmshTools::Read(modelElasticity2D, "../../rectangle.msh");

    // TODO: Disabled because now we recalculate the stiffness matrix contribution per node (based on phi)
    // Disables memory on elements
    // (makes sure elemental stiffness is not recalculated)
    // for(auto& [_, pWeakForm] : modelElasticity2D.MaterialVector()) {
        // govEquationElasticity2D->SetHasMemory(false);
    // }
}

void SetupBoundaryConditionsElasticity2DMBB(CompMesh &modelElasticity2D)
{
    MatrixDouble val1(2, 2);
    VecDouble val2(2);
    constexpr auto kVolumeMatId = 15;
    auto *govEquationElasticity2D = new Elasticity2D(kVolumeMatId, 1e3, 0.3);
    modelElasticity2D.InsertMaterial(govEquationElasticity2D);

    val1.setZero();
    val2.setZero();
    val2[0] = 1.0;
    constexpr auto kEngasteMatId = 16;
    auto *engasteBC = new L2Projection(kEngasteMatId, 2, BoundaryConditionType::kDirectionalHomogeneousDirichlet, val1, val2);
    modelElasticity2D.InsertMaterial(engasteBC);

    val1.setZero();
    val2.setZero();
    val2[1] = 1.0;
    constexpr auto kBottomRight = 17;
    auto *bottomRightBC = new L2Projection(kBottomRight, 2, BoundaryConditionType::kDirectionalHomogeneousDirichlet, val1, val2);
    modelElasticity2D.InsertMaterial(bottomRightBC);

    val2.setZero();
    val1.setZero();
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

    val1.setZero();
    val2.setZero();
    val2[1] = -1;
    constexpr auto kLoadMatId = 23;
    auto *El2D2 = new L2Projection(kLoadMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelElasticity2D.InsertMaterial(El2D2);
    GmshTools::Read(modelElasticity2D, "../../rectangle_mbb.msh");
}

void SetupBoundaryConditionsElasticity2DBridge(CompMesh &modelElasticity2D)
{
    MatrixDouble val1(2, 2);
    VecDouble val2(2);
    constexpr auto kVolumeMatId = 15;
    auto *govEquationElasticity2D = new Elasticity2D(kVolumeMatId, 1e3, 0.3);
    modelElasticity2D.InsertMaterial(govEquationElasticity2D);

    // val1.setZero();
    // val2.setZero();
    // constexpr auto kLeftMatId = 16;
    // auto *leftBC = new L2Projection(kLeftMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    // modelElasticity2D.InsertMaterial(leftBC);

    val1.setZero();
    val2.setZero();
    val2[1] = 1.0;
    constexpr auto kBottomRight = 17;
    auto *bottomRightBC = new L2Projection(kBottomRight, 2, BoundaryConditionType::kDirectionalHomogeneousDirichlet, val1, val2);
    modelElasticity2D.InsertMaterial(bottomRightBC);

    // val2.setZero();
    // val1.setZero();
    // val2[0] = 0.;
    // constexpr auto kFreeMatTopId = 18;
    // auto *freeBC = new L2Projection(kFreeMatTopId, 2, BoundaryConditionType::kNeumann, val1, val2);
    // modelElasticity2D.InsertMaterial(freeBC);

    // val1.setZero();
    // val2.setZero();
    // constexpr auto kFreeMatBottomId = 19;
    // auto *freeBCEl2D = new L2Projection(kFreeMatBottomId, 2, BoundaryConditionType::kNeumann, val1, val2);
    // modelElasticity2D.InsertMaterial(freeBCEl2D);

    val1.setZero();
    val2.setZero();
    val2[1] = -40;
    constexpr auto kMiddleBottom = 23;
    auto *El2D2 = new L2Projection(kMiddleBottom, 2, BoundaryConditionType::kNeumann, val1, val2);
    modelElasticity2D.InsertMaterial(El2D2);

    val1.setZero();
    val2.setZero();
    constexpr auto kLeftBottom = 24;
    auto *bottomLeftBC = new L2Projection(kLeftBottom, 2, BoundaryConditionType::kDirichlet, val1, val2);
    modelElasticity2D.InsertMaterial(bottomLeftBC);

    GmshTools::Read(modelElasticity2D, "../../rectangle_bridge.msh");
}


int main()
{
    Grid phi(Nx, vector<double>(Ny)), phi_new(Nx, vector<double>(Ny));
    initialize(phi);

    std::unique_ptr<CompMesh> modelElasticity2D = std::make_unique<CompMesh>();
    // SetupBoundaryConditionsElasticity2DBridge(*modelElasticity2D);
    SetupBoundaryConditionsElasticity2DMBB(*modelElasticity2D);
    LinearAnalysis anElasticity2D(modelElasticity2D.get(), SolverType::ELU);

    // RunEso(*modelElasticity2D, anElasticity2D, 0.55);
    // return 0;

    int i = 0;
    std::vector<std::string> ScalarNamesElasticity2D, VectorNamesElasticity2D;
    ScalarNamesElasticity2D = {"Compliance", "ComplianceSensibility", "WeightFunction"};
    VectorNamesElasticity2D = {"Displacement", "Stress","Strain"};

    int step = 0;
    Grid strain_energy(Nx, vector<double>(Ny));
    Grid mech_coupling_normalized(Nx, vector<double>(Ny));

    write_vtu(phi, "phi_" + std::to_string(0) + ".vtu");

    while(true)
    {
        // std::cout << "********** OPTIMIZATION STEP " << step << " **********\n";
        for (size_t inode = 0; inode < modelElasticity2D->NNodes(); inode++){
            const auto nodeVector = modelElasticity2D->NodeVec()[inode]->getCoordinates();
            const int i = nodeVector.x()/dx;
            const int j = nodeVector.y()/dy;
            if(i < 0 || i > Nx || j < 0 || j> Ny) {
                std::cout << "i: " << i << " j: " << j << " inode: " << inode << std::endl;
                PanicButton();
            }
            modelElasticity2D->NodeVec()[inode]->setWeightFunction(std::min(std::max(std::pow(phi[i][j], 3), 1e-2), 1.0));
        }
        for (int64_t i_el = 0; i_el < modelElasticity2D->NElements(); i_el++){
            auto &elemElas2D = modelElasticity2D->ElementVec()[i_el];            
            if (elemElas2D->Dimension() != modelElasticity2D->Dimension()) continue;
            elemElas2D->setIntegPointWeightFunction();
        }
        anElasticity2D.Run();
        if(i % SAVE_EVERY == 0) {
            VTUGenerator::PrintResults(modelElasticity2D.get(), "cantilever_2d_beam", ScalarNamesElasticity2D, VectorNamesElasticity2D, {}, i);
        }


        for (int64_t i_el = 0; i_el < modelElasticity2D->NElements(); i_el++){
            auto elemElas2D = modelElasticity2D->ElementVec()[i_el];

            const auto geoCenter = elemElas2D->GetGeometricCenter(*modelElasticity2D);
            
            if (elemElas2D->Dimension() != modelElasticity2D->Dimension())
            {
                continue;
            }
                    
            // for (int inode = 0; inode < elemElas2D->NElNodes(); inode++){
            //     auto xparametric = elemElas2D->NodeCoord(inode);
            //     auto connect = modelElasticity2D->ElementVec()[i_el]->getConnectivity();
            //     elemElas2D->IntegrationData().fAdimCoord = xparametric;
            //     if (!elemElas2D->IntegrationData().fNeedsSol || !elemElas2D->IntegrationData().fNeedsDSol){
            //         elemElas2D->IntegrationData().fNeedsSol = true;
            //         elemElas2D->IntegrationData().fNeedsDSol = true;
            //         elemElas2D->IntegrationData().fSol.resize(elemElas2D->GetWeakForm()->NState());
            //         elemElas2D->IntegrationData().fDSolDx.resize(elemElas2D->GetWeakForm()->NState(),elemElas2D->Dimension());
            //     }
            //     elemElas2D->ComputeJacobian();
            //     elemElas2D->ComputeSpatialDerivatives();
            //     elemElas2D->interpolateSolution();
            //     // const auto compliance_var_idx = elemElas2D->GetWeakForm()->VariableIndex("Compliance");
            //     // const auto nvar = elemElas2D->GetWeakForm()->NSolutionVariables(compliance_var_idx);
            //     // VecDouble compl_vec(nvar);
            //     // elemElas2D->Solution(compliance_var_idx, compl_vec);
            // }
            // const auto compliance_var_idx = elemElas2D->GetWeakForm()->VariableIndex("ComplianceSensibility");
            const auto compliance_var_idx = elemElas2D->GetWeakForm()->VariableIndex("ComplianceSensibility");
            const auto nvar = elemElas2D->GetWeakForm()->NSolutionVariables(compliance_var_idx);
            VecDouble compl_vec(nvar);

            const int i = (geoCenter.x()/dx) + 1;
            const int j = (geoCenter.y()/dy) + 1;
            elemElas2D->Solution(compliance_var_idx, compl_vec);
            if(i < 0 || i > Nx || j < 0 || j> Ny) {
                std::cout << "i: " << i << " j: " << j << " inode: " << i_el << std::endl;
                PanicButton();
            }
            strain_energy[i][j] = compl_vec[0];
        }

        for(int i = 0; i < PHASEFIELD_STEPS_PER_ELASTICITY; i++) {
            update(phi, phi_new, strain_energy, mech_coupling_normalized, step);
            phi.swap(phi_new);
            step++;
        }

        if(step % SAVE_EVERY == 0) {
            write_vtu(phi, "phi_" + std::to_string(step) + ".vtu");
            write_vtu(mech_coupling_normalized, "mech_coupling_normalized_" + std::to_string(step) + ".vtu");
            // VTUGenerator::PrintResults(modelElasticity2D.get(), "cantilever_2d_beam", ScalarNamesElasticity2D, VectorNamesElasticity2D, {}, i);
        }

        if(step > MAX_STEPS) {
            break;
        }

        i++;
    }

    return 0;
}
