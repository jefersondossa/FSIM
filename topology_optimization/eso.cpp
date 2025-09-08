#include "eso.hpp"
#include <LinearAnalysis.h>
#include <GmshTools.h>
#include <VTUGenerator.h>
#include <L2Projection.h>
#include <Elasticity3D.h>
#include <memory>
#include <unordered_set>

void RunEso(CompMesh& model, LinearAnalysis& an, double target_final_vol)
{
    // Disables memory on elements
    // (makes sure elemental stiffness is not recalculated)
    for(auto& [_, pWeakForm] : model.MaterialVector()) {
            pWeakForm->SetHasMemory(false);
    }

    an.Run();

    std::vector<std::string> ScalarNames, VectorNames;
    VectorNames = {"Displacement", "Stress"};
    ScalarNames = {"Compliance"};
    VTUGenerator::PrintResults(&model, "cantilever_2d_beam", ScalarNames, VectorNames);

    int64_t already_removed = 0;
    const double min_percentage_vol = target_final_vol;

    const auto n_elements = model.NElements();

    std::vector<std::pair<int64_t, double>> compliances{n_elements, std::make_pair<int64_t, double>(0, 0.0)};
    std::vector<double> eles_rho(model.NElements(), 1.0);
    const auto n_mat_eles = std::count_if(model.ElementVec().begin(), model.ElementVec().end(), 
        [model_dim = model.Dimension()](Element* pElement) {return pElement && pElement->Dimension() == model_dim; });
    double vol_percentage = (n_mat_eles - already_removed) / ((double)n_mat_eles);
    int index_it = 0;
    std::vector<double> elemental_compliances;
    elemental_compliances.resize(model.NElements());

    const double filter_radius = 0.05;

    while (vol_percentage > min_percentage_vol)
    {
        for (int64_t i_el = 0; i_el < model.NElements(); i_el++)
        {
            auto elem = model.ElementVec()[i_el];

            if (elem->Dimension() != model.Dimension())
            {
                elemental_compliances[i_el] = std::numeric_limits<double>::max();
                continue;
            }
            const auto compliance_var_idx = elem->GetWeakForm()->VariableIndex("Compliance");
            const auto nvar = elem->GetWeakForm()->NSolutionVariables(compliance_var_idx);
            VecDouble sol(nvar);
            elem->Solution(compliance_var_idx, sol);
            elemental_compliances[i_el] = sol[0];
        }

        for (int64_t i_el = 0; i_el < model.NElements(); i_el++)
        {
            if(elemental_compliances[i_el] == std::numeric_limits<double>::max())
            {
                compliances[i_el].first = i_el;
                compliances[i_el].second = elemental_compliances[i_el];
                continue;
            }

            const auto elem = model.ElementVec()[i_el];
            const VecDouble center_of_mass = elem->GetGeometricCenter(model);

            double sum = 0.0;

            std::vector<std::size_t> neighbors_to_visit{};
            std::unordered_set<std::size_t> neighbors_visited{};

            neighbors_to_visit.push_back(i_el);

            double compliance_val = 0.0;
            while (!neighbors_to_visit.empty())
            {
                const auto iel_neighbor = neighbors_to_visit.front();
                if(neighbors_visited.find(iel_neighbor) != neighbors_visited.end())
                {
                    neighbors_to_visit.front() = neighbors_to_visit.back();
                    neighbors_to_visit.pop_back();
                    continue;
                }

                if(elemental_compliances[iel_neighbor] != std::numeric_limits<double>::max())
                {
                    const auto elem_visiting = model.ElementVec()[iel_neighbor];
                    const VecDouble center_of_mass_neighbor = elem_visiting->GetGeometricCenter(model);

                    const auto distance = (center_of_mass_neighbor - center_of_mass).norm();
                    if (distance < filter_radius)
                    {
                        const auto fac = std::max(0.0, filter_radius - distance);
                        sum += fac;
                        compliance_val += fac * elemental_compliances[iel_neighbor] * eles_rho[iel_neighbor];

                        for (std::size_t i_neighbor_of_neighbor = 0; i_neighbor_of_neighbor < elem_visiting->getNumberOfNeighborElements(); i_neighbor_of_neighbor++)
                        {
                            neighbors_to_visit.push_back(elem_visiting->getNeighborElement(i_neighbor_of_neighbor));
                        }
                    }
                }
                neighbors_visited.insert(iel_neighbor);
                neighbors_to_visit.front() = neighbors_to_visit.back();
                neighbors_to_visit.pop_back();
            }

            const auto compliance = compliance_val / sum;

            compliances[i_el].first = i_el;
            compliances[i_el].second = compliance;
        }

        std::sort(compliances.begin(), compliances.end(), [](const auto &a, const auto &b)
                  { return a.second < b.second; });

        size_t n_eles_to_be_removed = static_cast<size_t>((n_mat_eles - already_removed) * 0.01);

        if (n_eles_to_be_removed == 1)
            break;

        std::cout << "Will remove " << n_eles_to_be_removed << " elements!" << '\n';

        for (size_t i = 0; i < compliances.size() && n_eles_to_be_removed > 0; i++)
        {
            const auto &compliance = compliances[i];
            const auto i_el = compliance.first;
            const auto elem = model.ElementVec()[i_el];
            const auto rho = eles_rho[i_el];

            if (rho != 0.01)
            {
                eles_rho[i_el] = 0.01;
                if (!elem->IntegrationData().fStiffnessMatrix.has_value())
                {
                    PanicButton();
                }

                elem->IntegrationData().fStiffnessMatrix.value() *= 0.01;
                already_removed++;
                n_eles_to_be_removed--;
            }
        }

        an.Run();

        vol_percentage = (n_mat_eles - already_removed) / ((double)n_mat_eles);
        VTUGenerator::PrintResults(&model, std::string{"cantilever_2d_beam"} + std::to_string(index_it), ScalarNames, VectorNames, { {"rho", eles_rho} });
        index_it++;
    }
}