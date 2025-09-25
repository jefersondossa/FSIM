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
#include <algorithm>
#include <bits/stdc++.h>

static double neumann_g(double x, double y, double t) { return 0.0; }
static double u0(double x, double y)
{
    // double cx = 0.5, cy = 0.5, sigma = 0.08;
    // double r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
    // return exp(-r2 / (2 * sigma * sigma));
    constexpr auto A = 0.2;
    const double V = 0.55;
    return V + A * sin(8 * M_PI * x) * sin(8 * M_PI * y);
}

// Utilities -----------------------------------------------------------------
struct Mesh
{
    int Nx, Ny;
    int numNodesX, numNodesY;
    int numNodes;
    int numElems;
    double hx, hy;
    std::vector<std::pair<double, double>> nodes;
    std::vector<std::array<int, 4>> elems;
    double minH;
    
    Mesh(const std::string& filename) {
        minH = 1e30;
        std::unique_ptr<CompMesh> meshData = std::make_unique<CompMesh>();
        MatrixDouble val1(1, 1);
        VecDouble val2(1);
        val1.setZero();
        val2.setZero();
        meshData->InsertMaterial(new PhaseField(15, 2));
        meshData->InsertMaterial(new L2Projection(16, 2, BoundaryConditionType::kDirichlet, val1, val2));
        meshData->InsertMaterial(new L2Projection(18, 2, BoundaryConditionType::kNeumann, val1, val2));
        meshData->InsertMaterial(new L2Projection(19, 2, BoundaryConditionType::kNeumann, val1, val2));
        meshData->InsertMaterial(new L2Projection(45, 2, BoundaryConditionType::kNeumann, val1, val2));
        GmshTools::Read(*meshData, filename);

        nodes.resize(meshData->NNodes());
        elems.resize(meshData->NElements());
        numNodes = meshData->NNodes();
        numElems = meshData->NElements();

        std::array<int, 4> elemNodesIds;

        for (int64_t iEle = 0; iEle < meshData->NElements(); iEle++)
        {
            auto elem = meshData->ElementVec()[iEle];

            assert(elem->NElNodes() == elemNodesIds.size());
            const auto connect = meshData->ElementVec()[iEle]->getConnectivity();

            for (int iNode = 0; iNode < elem->NElNodes(); iNode++)
            {
                const auto x = meshData->NodeVec()[connect[iNode]]->getCoordinates()[0];
                const auto y = meshData->NodeVec()[connect[iNode]]->getCoordinates()[1];

                nodes[connect[iNode]] = std::make_pair(x, y);

                elemNodesIds[iNode] = connect[iNode];
            }

            minH = std::min(std::fabs( - nodes[connect[1]].first), minH);
            const double x0 = nodes[connect[0]].first;
            const double y0 = nodes[connect[0]].second;

            const double x1 = nodes[connect[1]].first;
            const double y1 = nodes[connect[1]].second;

            const double x2 = nodes[connect[2]].first;
            const double y2 = nodes[connect[2]].second;

            const double x3 = nodes[connect[3]].first;
            const double y3 = nodes[connect[3]].second;

            const double area = 0.5 * fabs(
                x0*y1 + x1*y2 + x2*y3 + x3*y0
            - y0*x1 - y1*x2 - y2*x3 - y3*x0
            );

            minH = std::min(minH, std::sqrt(area));

            elems[iEle] = {elemNodesIds[0], elemNodesIds[1], elemNodesIds[2], elemNodesIds[3]};
        }
    }
    Mesh(int Nx_, int Ny_) : Nx(Nx_), Ny(Ny_)
    {
        numNodesX = Nx + 1;
        numNodesY = Ny + 1;
        numNodes = numNodesX * numNodesY;
        numElems = Nx * Ny;
        hx = 1.0 / Nx;
        hy = 1.0 / Ny;
        minH = std::min(hx, hy);
        nodes.reserve(numNodes);
        for (int j = 0; j < numNodesY; ++j)
        {
            for (int i = 0; i < numNodesX; ++i)
                nodes.emplace_back(i * hx, j * hy);
        }
        elems.reserve(numElems);
        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                int n0 = j * numNodesX + i;
                int n1 = n0 + 1;
                int n3 = n0 + numNodesX;
                int n2 = n3 + 1;
                elems.push_back({n0, n1, n2, n3});
            }
        }
    }
};

inline void shape_functions(double xi, double eta, double N[4], double dN_dxi[4], double dN_deta[4])
{
    N[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
    N[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
    N[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
    N[3] = 0.25 * (1.0 - xi) * (1.0 + eta);
    dN_dxi[0] = -0.25 * (1.0 - eta);
    dN_dxi[1] = 0.25 * (1.0 - eta);
    dN_dxi[2] = 0.25 * (1.0 + eta);
    dN_dxi[3] = -0.25 * (1.0 + eta);
    dN_deta[0] = -0.25 * (1.0 - xi);
    dN_deta[1] = -0.25 * (1.0 + xi);
    dN_deta[2] = 0.25 * (1.0 + xi);
    dN_deta[3] = 0.25 * (1.0 - xi);
}

void write_vtu(const Mesh &mesh, const std::vector<double> &u, int step, double time)
{
    std::ostringstream fname;
    fname << "solution_" << std::setw(5) << std::setfill('0') << step << ".vtu";
    std::ofstream out(fname.str());
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "<UnstructuredGrid>\n";
    out << "<Piece NumberOfPoints=\"" << mesh.numNodes << "\" NumberOfCells=\"" << mesh.numElems << "\">\n";

    // Points
    out << "<Points>\n";
    out << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (auto [x, y] : mesh.nodes)
    {
        out << x << " " << y << " 0\n";
    }
    out << "</DataArray>\n";
    out << "</Points>\n";

    // Cells
    out << "<Cells>\n";
    out << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (auto en : mesh.elems)
    {
        out << en[0] << " " << en[1] << " " << en[2] << " " << en[3] << "\n";
    }
    out << "</DataArray>\n";
    out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int off = 0;
    for (int e = 0; e < mesh.numElems; ++e)
    {
        off += 4;
        out << off << "\n";
    }
    out << "</DataArray>\n";
    out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int e = 0; e < mesh.numElems; ++e)
        out << 9 << "\n"; // VTK_QUAD=9
    out << "</DataArray>\n";
    out << "</Cells>\n";

    // Point data (solution)
    out << "<PointData Scalars=\"u\">\n";
    out << "<DataArray type=\"Float64\" Name=\"u\" format=\"ascii\">\n";
    for (double val : u)
        out << val << "\n";
    out << "</DataArray>\n";
    out << "</PointData>\n";

    out << "</Piece>\n";
    out << "</UnstructuredGrid>\n";
    out << "</VTKFile>\n";
    std::cerr << "Wrote " << fname.str() << " at time=" << time << "\n";
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
    // val2[0] = 1.;
    // constexpr auto kEngasteMatId = 22;
    constexpr auto kEngasteMatId = 16;
    // auto *engasteBC = new L2Projection(kEngasteMatId, 2, BoundaryConditionType::kDirectionalHomogeneousDirichlet, val1, val2);
    auto *engasteBC = new L2Projection(kEngasteMatId, 2, BoundaryConditionType::kDirichlet, val1, val2);
    modelElasticity2D.InsertMaterial(engasteBC);
    // {
    //     val1.setZero();
    //     val2.setZero();
    //     constexpr auto kEngasteMatId = 22;
    //     auto *engasteBC = new L2Projection(kEngasteMatId, 2, BoundaryConditionType::kDirichlet, val1, val2);
    //     modelElasticity2D.InsertMaterial(engasteBC);
    // }
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
    {
        val2[0] = 10.;
        constexpr auto kRightMatBottomId = 45;
        auto *freeBCEl2D = new L2Projection(kRightMatBottomId, 2, BoundaryConditionType::kNeumann, val1, val2);
        modelElasticity2D.InsertMaterial(freeBCEl2D);
    }

    // {
    //     constexpr auto kFreeMatBottomId = 23;
    //     auto *freeBCEl2D = new L2Projection(kFreeMatBottomId, 2, BoundaryConditionType::kNeumann, val1, val2);
    //     modelElasticity2D.InsertMaterial(freeBCEl2D);
    // }
    // val1.setZero();
    // val2.setZero();
    // val2[0] = 0.;
    // constexpr auto kExLoadMatId = 21;
    // auto *El2D = new L2Projection(kExLoadMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    // modelElasticity2D.InsertMaterial(El2D);
    // val1.setZero();
    // val2.setZero();
    // // val2[0] = -10;
    // val2[1] = -10;
    // constexpr auto kLoadMatId = 17;
    // auto *El2D2 = new L2Projection(kLoadMatId, 2, BoundaryConditionType::kNeumann, val1, val2);
    // modelElasticity2D.InsertMaterial(El2D2);
    GmshTools::Read(modelElasticity2D, "../../rectangle.msh");

    // TODO: Disabled because now we recalculate the stiffness matrix contribution per node (based on phi)
    // Disables memory on elements
    // (makes sure elemental stiffness is not recalculated)
    // for(auto& [_, pWeakForm] : modelElasticity2D.Materialstd::Vector()) {
        // govEquationElasticity2D->SetHasMemory(false);
    // }
}

int main()
{
    std::unique_ptr<CompMesh> modelElasticity2D = std::make_unique<CompMesh>();
    SetupBoundaryConditionsElasticity2DCantileverRightBottom(*modelElasticity2D);
    LinearAnalysis anElasticity2D(modelElasticity2D.get(), SolverType::ELU);

    Mesh mesh("../../rectangle.msh");
    double h = mesh.minH;//min(mesh.hx, mesh.hy);
    double dt = 0.25 * h * h; // conservative for stability
    double T = 0.5;
    double kappa = 0.05; // example diffusivity
    int numSteps = std::max(1, (int)ceil(T / dt));
    dt = T / numSteps;

    std::vector<double> u(mesh.numNodes, 0.0), u_new(mesh.numNodes, 0.0), Mdiag(mesh.numNodes, 0.0), rhs(mesh.numNodes, 0.0);

    // initial condition
    for (int a = 0; a < mesh.numNodes; ++a)
    {
        auto [x, y] = mesh.nodes[a];
        u[a] = u0(x, y);
    }

    // compute lumped mass diag and stiffness per element
    const int nq = 2;
    double gp[nq] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)}, gw[nq] = {1.0, 1.0};
    for (int e = 0; e < mesh.numElems; ++e)
    {
        auto enodes = mesh.elems[e];
        double x[4], y[4];
        for (int i = 0; i < 4; ++i)
        {
            auto [xx, yy] = mesh.nodes[enodes[i]];
            x[i] = xx;
            y[i] = yy;
        }
        double Me[4][4] = {0};
        for (int ig = 0; ig < nq; ++ig)
            for (int jg = 0; jg < nq; ++jg)
            {
                double xi = gp[ig], eta = gp[jg], w = gw[ig] * gw[jg];
                double N[4], dNdxi[4], dNdeta[4];
                shape_functions(xi, eta, N, dNdxi, dNdeta);
                double J11 = 0, J12 = 0, J21 = 0, J22 = 0;
                for (int a = 0; a < 4; ++a)
                {
                    J11 += dNdxi[a] * x[a];
                    J12 += dNdxi[a] * y[a];
                    J21 += dNdeta[a] * x[a];
                    J22 += dNdeta[a] * y[a];
                }
                double detJ = J11 * J22 - J12 * J21;
                for (int a = 0; a < 4; ++a)
                    for (int b = 0; b < 4; ++b)
                        Me[a][b] += N[a] * N[b] * detJ * w;
            }
        for (int a = 0; a < 4; ++a)
        {
            double rowSum = 0;
            for (int b = 0; b < 4; ++b)
                rowSum += Me[a][b];
            Mdiag[enodes[a]] += rowSum;
        }
    }

    // store stiffness matrices per element
    std::vector<std::array<double, 16>> Ke_list(mesh.numElems);
    for (int e = 0; e < mesh.numElems; ++e)
    {
        auto enodes = mesh.elems[e];
        double x[4], y[4];
        for (int i = 0; i < 4; ++i)
        {
            auto [xx, yy] = mesh.nodes[enodes[i]];
            x[i] = xx;
            y[i] = yy;
        }
        double Ke[4][4];
        memset(Ke, 0, sizeof(Ke));
        for (int ig = 0; ig < nq; ++ig)
            for (int jg = 0; jg < nq; ++jg)
            {
                double xi = gp[ig], eta = gp[jg], w = gw[ig] * gw[jg];
                double N[4], dNdxi[4], dNdeta[4];
                shape_functions(xi, eta, N, dNdxi, dNdeta);
                double J11 = 0, J12 = 0, J21 = 0, J22 = 0;
                for (int a = 0; a < 4; ++a)
                {
                    J11 += dNdxi[a] * x[a];
                    J12 += dNdxi[a] * y[a];
                    J21 += dNdeta[a] * x[a];
                    J22 += dNdeta[a] * y[a];
                }
                double detJ = J11 * J22 - J12 * J21;
                double invJ11 = J22 / detJ, invJ12 = -J12 / detJ, invJ21 = -J21 / detJ, invJ22 = J11 / detJ;
                double dNdx[4], dNdy[4];
                for (int a = 0; a < 4; ++a)
                {
                    dNdx[a] = invJ11 * dNdxi[a] + invJ12 * dNdeta[a];
                    dNdy[a] = invJ21 * dNdxi[a] + invJ22 * dNdeta[a];
                }
                for (int a = 0; a < 4; ++a)
                    for (int b = 0; b < 4; ++b)
                        Ke[a][b] += kappa * (dNdx[a] * dNdx[b] + dNdy[a] * dNdy[b]) * detJ * w;
            }
        std::array<double, 16> Kflat;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                Kflat[a * 4 + b] = Ke[a][b];
        Ke_list[e] = Kflat;
    }

    // compute total mass to enforce (V0). We'll choose V0 = initial mass (conserve mass)
    double sumM = 0.0;
    for (double m : Mdiag)
        sumM += m;
    auto compute_total_mass = [&](const std::vector<double> &vec)
    { double s=0; for(int i=0;i<(int)vec.size(); ++i) s += Mdiag[i]*vec[i]; return s; };
    double V0 = compute_total_mass(u);
    std::cerr << "Initial total mass V0 = " << V0 << "   sumM=" << sumM << "";

    // write initial
    double t = 0.0;
    write_vtu(mesh, u, 0, t);

    // time loop with explicit Euler + global Lagrange multiplier (enforced as constant correction)
    for (int step = 0; step < numSteps; ++step)
    {
        for (size_t inode = 0; inode < modelElasticity2D->NNodes(); inode++)
        {
            const auto solutionPhaseField = u_new[inode];
            modelElasticity2D->NodeVec()[inode]->setWeightFunction(std::min(std::max(std::pow(solutionPhaseField, 3), 1e-2), 1.0));
        }

        for (int64_t i_el = 0; i_el < modelElasticity2D->NElements(); i_el++){
            auto &elemElas2D = modelElasticity2D->ElementVec()[i_el];            
            if (elemElas2D->Dimension() != modelElasticity2D->Dimension()) continue;
            elemElas2D->setIntegPointWeightFunction();
        }
        std::vector<double> compliances;
        for (int64_t i_el = 0; i_el < modelElasticity2D->NElements(); i_el++)
        {
            auto elemElas2D = modelElasticity2D->ElementVec()[i_el];
            auto elemPhaseField = mesh.elems[i_el];

            if (elemElas2D->Dimension() != modelElasticity2D->Dimension())
            {
                continue;
            }

            for (int inode = 0; inode < elemElas2D->NElNodes(); inode++)
            {
                auto xparametric = elemElas2D->NodeCoord(inode);
                auto connect = modelElasticity2D->ElementVec()[i_el]->getConnectivity();
                elemElas2D->IntegrationData().fAdimCoord = xparametric;
                if (!elemElas2D->IntegrationData().fNeedsSol || !elemElas2D->IntegrationData().fNeedsDSol)
                {
                    elemElas2D->IntegrationData().fNeedsSol = true;
                    elemElas2D->IntegrationData().fNeedsDSol = true;
                    elemElas2D->IntegrationData().fSol.resize(elemElas2D->GetWeakForm()->NState());
                    elemElas2D->IntegrationData().fDSolDx.resize(elemElas2D->GetWeakForm()->NState(), elemElas2D->Dimension());
                }
                elemElas2D->ComputeJacobian();
                elemElas2D->ComputeSpatialDerivatives();
                elemElas2D->interpolateSolution();
                const auto compliance_var_idx = elemElas2D->GetWeakForm()->VariableIndex("ComplianceSensibility");
                const auto nvar = elemElas2D->GetWeakForm()->NSolutionVariables(compliance_var_idx);
                VecDouble compl_vec(nvar);
                elemElas2D->Solution(compliance_var_idx, compl_vec);
                compliances.push_back(compl_vec[0]);
            }
        }

        const auto maxCompliance = *std::max_element(compliances.begin(), compliances.end());

        for(auto& compliance : compliances)
        {
            compliance /= maxCompliance;
        }

        assert(compliances.size() == u.size());

        const auto f_source = [&compliances](int64_t iNode, double u) -> double
        {
            double res = 15 * 2 * u * (1.0 * u) * (1.0 - u);
            res += compliances[iNode];
            return res;
        };

        fill(rhs.begin(), rhs.end(), 0.0);
        // K*u contribution
        for (int e = 0; e < mesh.numElems; ++e)
        {
            auto enodes = mesh.elems[e];
            auto Kflat = Ke_list[e];
            double u_local[4];
            for (int a = 0; a < 4; ++a)
                u_local[a] = u[enodes[a]];
            double Ku[4] = {0, 0, 0, 0};
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    Ku[a] += Kflat[a * 4 + b] * u_local[b];
            for (int a = 0; a < 4; ++a)
                rhs[enodes[a]] -= Ku[a];
        }
        // source term (note: evaluate u at quadrature point by interpolation)
        for (int e = 0; e < mesh.numElems; ++e)
        {
            auto enodes = mesh.elems[e];
            double x[4], y[4];
            for (int i = 0; i < 4; ++i)
            {
                auto [xx, yy] = mesh.nodes[enodes[i]];
                x[i] = xx;
                y[i] = yy;
            }
            for (int ig = 0; ig < nq; ++ig)
                for (int jg = 0; jg < nq; ++jg)
                {
                    double xi = gp[ig], eta = gp[jg], w = gw[ig] * gw[jg];
                    double N[4], dNdxi[4], dNdeta[4];
                    shape_functions(xi, eta, N, dNdxi, dNdeta);
                    double J11 = 0, J12 = 0, J21 = 0, J22 = 0;
                    for (int a = 0; a < 4; ++a)
                    {
                        J11 += dNdxi[a] * x[a];
                        J12 += dNdxi[a] * y[a];
                        J21 += dNdeta[a] * x[a];
                        J22 += dNdeta[a] * y[a];
                    }
                    double detJ = J11 * J22 - J12 * J21;
                    double X = 0, Y = 0, uq = 0;
                    for (int a = 0; a < 4; ++a)
                    {
                        X += N[a] * x[a];
                        Y += N[a] * y[a];
                        uq += N[a] * u[enodes[a]];
                    }
                    for (int a = 0; a < 4; ++a) {
                        const double fval = f_source(enodes[a], uq);
                        rhs[enodes[a]] += N[a] * fval * detJ * w;
                    }
                }
            // Neumann flux (as before) -- omitted here for brevity (zero flux in current setup)
        }
        // explicit Euler provisional update
        for (int a = 0; a < mesh.numNodes; ++a)
        {
            double m = Mdiag[a];
            if (m <= 0)
                m = 1e-16;
            u_new[a] = u[a] + dt * (rhs[a] / m);
        }

        // Enforce global mass conservation using a single Lagrange multiplier that acts as a constant correction
        // We look for constant c such that sum_i M_i*(u_new_i + c) = V0  --> c = (V0 - mass_new)/sumM
        double mass_new = compute_total_mass(u_new);
        double c = (V0 - mass_new) / sumM; // this is the constant correction added to all nodal values
        for (int a = 0; a < mesh.numNodes; ++a)
            u_new[a] += c;

        // swap and advance
        u.swap(u_new);
        t += dt;

        // diagnostics
        double mass_after = compute_total_mass(u);
        std::cerr << "step " << step + 1 << " t=" << t << " mass(after)=" << mass_after << " c=" << c << " max(u)=" << (*max_element(u.begin(), u.end())) << "\n";

        // write vtu
        write_vtu(mesh, u, step + 1, t);
    }

    return 0;
}
