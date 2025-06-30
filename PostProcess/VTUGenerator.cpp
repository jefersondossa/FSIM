#include "VTUGenerator.h"
#include "LagrangeMultiplier.h"
#include "NullWeakForm.h"

void VTUGenerator::PrintResultsGraph(CompMesh *cmesh, std::string filename, std::vector<std::string> &scalnames, std::vector<std::string> &vecnames, int step){

    auto graphmesh = cmesh->GetGraphMesh();

    //    std::cout << "Printing Velocity Results" << '\n';
    std::string s = filename + std::to_string(step) + ".vtu";
    
    std::fstream output_v(s.c_str(), std::ios_base::out);

    output_v << "<?xml version=\"1.0\"?>" << '\n'
             << "<VTKFile type=\"UnstructuredGrid\">" << '\n'
             << "  <UnstructuredGrid>" << '\n'
             << "  <Piece NumberOfPoints=\"" << graphmesh->NNodes()
             << "\"  NumberOfCells=\"" << graphmesh->NElements()
             << "\">" << '\n';

    //WRITE NODAL COORDINATES
    output_v << "    <Points>" << '\n'
             << "      <DataArray type=\"Float64\" "
             << "NumberOfComponents=\"3\" format=\"ascii\">" << '\n';

    for (int i=0; i<graphmesh->NNodes(); i++){
        auto x = graphmesh->Node(i);
        output_v << x[0] << " " << x[1] << " " << x[2] << '\n';        
    };

    output_v << "      </DataArray>" << '\n'
             << "    </Points>" << '\n';
    
    //WRITE ELEMENT CONNECTIVITY
    output_v << "    <Cells>" << '\n'
             << "      <DataArray type=\"Int32\" "
             << "Name=\"connectivity\" format=\"ascii\">" << '\n';
    
    for (int i=0; i<graphmesh->NElements(); i++){
        auto connec=graphmesh->Connect(i);
        for (int k = 0; k < connec.size(); k++)
        {
            output_v << connec[k] << " ";
        }
        output_v << '\n';
    };
    output_v << "      </DataArray>" << '\n';
  
    //WRITE OFFSETS IN DATA ARRAY
    output_v << "      <DataArray type=\"Int32\""
             << " Name=\"offsets\" format=\"ascii\">" << '\n';
    
    int aux = 0;
    for (int i=0; i<graphmesh->NElements(); i++){
        aux += graphmesh->Connect(i).size();
        output_v << aux << '\n';
    };
    output_v << "      </DataArray>" << '\n';
  
    //WRITE ELEMENT TYPES
    output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
             << "format=\"ascii\">" << '\n';

    for (int i=0; i<graphmesh->NElements(); i++){
        output_v << graphmesh->ElType(i) << '\n';
    };

    output_v << "      </DataArray>" << '\n'
             << "    </Cells>" << '\n';

    //WRITE NODAL RESULTS
    output_v << "    <PointData>" << '\n';

    //Scalar values
    std::map<int64_t,std::vector<VecDouble>> scalSol;
    std::map<int64_t,std::vector<VecDouble>> vectSol;
    // std::vector<std::vector<VecDouble,scalnames.size()>,graphmesh->NNodes()> scalSol;
    for (int64_t iel = 0; iel < graphmesh->NElements(); iel++){
        auto compel = cmesh->ElementVec()[iel];
        auto graphconnect = graphmesh->Connect(iel);
        if (compel->Dimension() != cmesh->Dimension()) continue;
        auto wf = compel->GetWeakForm();
        auto lagm = dynamic_cast<LagrangeMultiplier*>(wf);
        auto null = dynamic_cast<NullWeakForm*>(wf);
        if (lagm || null) continue;

        for (int inode = 0; inode < compel->NElNodes(); inode++){
            auto xparametric = compel->NodeCoord(inode);
            compel->IntegrationData().fAdimCoord = xparametric;
            if (!compel->IntegrationData().fNeedsSol || !compel->IntegrationData().fNeedsDSol){
                compel->IntegrationData().fNeedsSol = true;
                compel->IntegrationData().fNeedsDSol = true;
                compel->IntegrationData().fSol.resize(compel->GetWeakForm()->NState());
                compel->IntegrationData().fDSolDx.resize(compel->GetWeakForm()->NState(),compel->Dimension());
            }
            compel->ComputeJacobian();
            compel->ComputeSpatialDerivatives();
            compel->interpolateSolution();
            compel->interpolateSolDerivatives();
            
            for (int iscal = 0; iscal < scalnames.size(); iscal++){
                int varindex = compel->GetWeakForm()->VariableIndex(scalnames[iscal]);
                int nvar = compel->GetWeakForm()->NSolutionVariables(varindex);
                VecDouble Sol(nvar);
                compel->Solution(varindex,Sol);
                scalSol[graphconnect[inode]].push_back(Sol);
            }
            for (int ivect = 0; ivect < vecnames.size(); ivect++){
                int varindex = compel->GetWeakForm()->VariableIndex(vecnames[ivect]);
                int nvar = compel->GetWeakForm()->NSolutionVariables(varindex);
                VecDouble Sol(nvar);
                compel->Solution(varindex,Sol);
                vectSol[graphconnect[inode]].push_back(Sol);
            }
        }
    }
    for (int64_t iel = 0; iel < graphmesh->NElements(); iel++){
        auto compel = cmesh->ElementVec()[iel];
        if (compel->Dimension() == cmesh->Dimension()) continue;
        auto graphconnect = graphmesh->Connect(iel);
        for (int i = 0; i < graphconnect.size(); i++){
            auto meshnode = graphmesh->GraphNodeToMeshNode(graphconnect[i]);
            auto graphnode = graphmesh->MeshNodeToGraphNode(meshnode);
            scalSol[graphconnect[i]] = scalSol[graphnode];
            vectSol[graphconnect[i]] = vectSol[graphnode];
        }
    }
    
    for (int iscal = 0; iscal < scalnames.size(); iscal++){
        output_v<< "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
            << "Name=\"" << scalnames[iscal] << "\" format=\"ascii\">" << '\n';
        for (int i=0; i<graphmesh->NNodes(); i++){
            if (scalSol[i].size()==0){
                output_v <<std::scientific<< "0" << '\n';
            }else{
                output_v <<std::scientific<< scalSol[i][iscal][0] << '\n';
            }
        }
        output_v << "      </DataArray> " << '\n';
    }
    
    for (int iscal = 0; iscal < vecnames.size(); iscal++){
        output_v<< "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
            << "Name=\"" << vecnames[iscal] << "\" format=\"ascii\">" << '\n';
        for (int i=0; i<graphmesh->NNodes(); i++){
            if (vectSol[i].size()==0){
                output_v <<std::scientific<< "0 0 0" << '\n';
            }else{
                output_v<<std::scientific << vectSol[i][iscal][0] << " " << vectSol[i][iscal][1] << " " << vectSol[i][iscal][2] << '\n';    
            }
        }
        output_v << "      </DataArray> " << '\n';
    }
    
    




    output_v << "    </PointData>" << '\n'; 

    //WRITE ELEMENT RESULTS
    output_v << "    <CellData>" << '\n';
    
    //Some element wise result
    output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                << "Name=\"Material\" format=\"ascii\">" << '\n';
    for (int i=0; i<graphmesh->NElements(); i++){
        auto compel = cmesh->ElementVec()[i];
        int matid = compel->GetWeakForm()->Id();
        output_v << matid << '\n';
    };
    output_v << "      </DataArray> " << '\n';
    
    output_v << "    </CellData>" << '\n'; 

    //FINALIZE OUTPUT FILE
    output_v << "  </Piece>" << '\n';
    output_v << "  </UnstructuredGrid>" << '\n'
             << "</VTKFile>" << '\n';

}

void VTUGenerator::PrintResults(CompMesh *cmesh, std::string filename, std::vector<std::string> &scalnames, std::vector<std::string> &vecnames, const std::vector<CustomCellField>& custom_cel_fields, int step){

    //    std::cout << "Printing Velocity Results" << '\n';
    std::string s = filename + std::to_string(step) + ".vtu";
    
    std::fstream output_v(s.c_str(), std::ios_base::out);

    output_v << "<?xml version=\"1.0\"?>" << '\n'
             << "<VTKFile type=\"UnstructuredGrid\">" << '\n'
             << "  <UnstructuredGrid>" << '\n'
             << "  <Piece NumberOfPoints=\"" << cmesh->NNodes()
             << "\"  NumberOfCells=\"" << cmesh->NElements()
             << "\">" << '\n';

    //WRITE NODAL COORDINATES
    output_v << "    <Points>" << '\n'
             << "      <DataArray type=\"Float64\" "
             << "NumberOfComponents=\"3\" format=\"ascii\">" << '\n';

    for (int i=0; i<cmesh->NNodes(); i++){
        auto x = cmesh->NodeVec()[i]->getCoordinates();
        for(std::size_t i = 0; i < 3; i++) {
            if(i < cmesh->Dimension())
            {
                output_v << x[i] << " ";    
            }
            else {
                output_v << 0.0 << " ";
            }
        }
        output_v << '\n';
    };

    output_v << "      </DataArray>" << '\n'
             << "    </Points>" << '\n';
    
    //WRITE ELEMENT CONNECTIVITY
    output_v << "    <Cells>" << '\n'
             << "      <DataArray type=\"Int32\" "
             << "Name=\"connectivity\" format=\"ascii\">" << '\n';
    
    for (int i=0; i<cmesh->NElements(); i++){
        auto connec=cmesh->ElementVec()[i]->getConnectivity();
        for (int k = 0; k < connec.size(); k++)
        {
            output_v << connec[k] << " ";
        }
        output_v << '\n';
    };
    output_v << "      </DataArray>" << '\n';
  
    //WRITE OFFSETS IN DATA ARRAY
    output_v << "      <DataArray type=\"Int32\""
             << " Name=\"offsets\" format=\"ascii\">" << '\n';
    
    int aux = 0;
    for (int i=0; i<cmesh->NElements(); i++){
        aux += cmesh->ElementVec()[i]->getConnectivity().size();
        output_v << aux << '\n';
    };
    output_v << "      </DataArray>" << '\n';
  
    //WRITE ELEMENT TYPES
    output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
             << "format=\"ascii\">" << '\n';

    for (int i=0; i<cmesh->NElements(); i++){
        output_v << cmesh->ElementVec()[i]->PrintType() << '\n';
    };

    output_v << "      </DataArray>" << '\n'
             << "    </Cells>" << '\n';

    //WRITE NODAL RESULTS
    output_v << "    <PointData>" << '\n';

    //Scalar values
    std::map<int64_t,std::vector<VecDouble>> scalSol;
    std::map<int64_t,std::vector<VecDouble>> vectSol;
    // std::vector<std::vector<VecDouble,scalnames.size()>,graphmesh->NNodes()> scalSol;
    for (int64_t iel = 0; iel < cmesh->NElements(); iel++){
        auto compel = cmesh->ElementVec()[iel];
        auto graphconnect = cmesh->ElementVec()[iel]->getConnectivity();
        if (compel->Dimension() != cmesh->Dimension()) continue;
        auto wf = compel->GetWeakForm();
        auto lagm = dynamic_cast<LagrangeMultiplier*>(wf);
        auto null = dynamic_cast<NullWeakForm*>(wf);
        if (lagm || null) continue;

        for (int inode = 0; inode < compel->NElNodes(); inode++){
            auto xparametric = compel->NodeCoord(inode);
            compel->IntegrationData().fAdimCoord = xparametric;
            if (!compel->IntegrationData().fNeedsSol || !compel->IntegrationData().fNeedsDSol){
                compel->IntegrationData().fNeedsSol = true;
                compel->IntegrationData().fNeedsDSol = true;
                compel->IntegrationData().fSol.resize(compel->GetWeakForm()->NState());
                compel->IntegrationData().fDSolDx.resize(compel->GetWeakForm()->NState(),compel->Dimension());
            }
            compel->ComputeJacobian();
            compel->ComputeSpatialDerivatives();
            compel->interpolateSolution();
            compel->interpolateSolDerivatives();
            if (compel->IntegrationData().fNeedsTimeDerivatives){
                compel->interpolateSolDTimeDerivatives();
            }
            
            for (int iscal = 0; iscal < scalnames.size(); iscal++){
                int varindex = compel->GetWeakForm()->VariableIndex(scalnames[iscal]);
                int nvar = compel->GetWeakForm()->NSolutionVariables(varindex);
                VecDouble Sol(nvar);
                compel->Solution(varindex,Sol);
                if (scalSol[graphconnect[inode]].size()==0){
                    scalSol[graphconnect[inode]].resize(scalnames.size());
                }
                scalSol[graphconnect[inode]][iscal]=Sol;
            }
            for (int ivect = 0; ivect < vecnames.size(); ivect++){
                int varindex = compel->GetWeakForm()->VariableIndex(vecnames[ivect]);
                int nvar = compel->GetWeakForm()->NSolutionVariables(varindex);
                VecDouble Sol(nvar);
                compel->Solution(varindex,Sol);
                if (vectSol[graphconnect[inode]].size()==0){
                    vectSol[graphconnect[inode]].resize(vecnames.size());
                }
                vectSol[graphconnect[inode]][ivect]=Sol;
            }
        }
    }
    // for (int64_t iel = 0; iel < cmesh->NElements(); iel++){
    //     auto compel = cmesh->ElementVec()[iel];
    //     if (compel->Dimension() == cmesh->Dimension()) continue;
    //     auto graphconnect = cmesh->ElementVec()[iel]->getConnectivity();
    //     for (int i = 0; i < graphconnect.size(); i++){
    //         auto meshnode = graphconnect[i];//graphmesh->GraphNodeToMeshNode(graphconnect[i]);
    //         auto graphnode = graphmesh->MeshNodeToGraphNode(meshnode);
    //         scalSol[graphconnect[i]] = scalSol[graphnode];
    //         vectSol[graphconnect[i]] = vectSol[graphnode];
    //     }
    // }
    
    for (int iscal = 0; iscal < scalnames.size(); iscal++){
        output_v<< "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
            << "Name=\"" << scalnames[iscal] << "\" format=\"ascii\">" << '\n';
        for (int i=0; i<cmesh->NNodes(); i++){
            if (scalSol[i].size()==0){
                output_v <<std::scientific<< "0" << '\n';
            }else{
                output_v <<std::scientific<< scalSol[i][iscal][0] << '\n';
            }
        }
        output_v << "      </DataArray> " << '\n';
    }
    
    for (int iscal = 0; iscal < vecnames.size(); iscal++){
        output_v<< "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
            << "Name=\"" << vecnames[iscal] << "\" format=\"ascii\">" << '\n';
        for (int i=0; i<cmesh->NNodes(); i++){
            if (vectSol[i].size()==0){
                output_v <<std::scientific<< "0 0 0" << '\n';
            }else{
                output_v<<std::scientific << vectSol[i][iscal][0] << " " << vectSol[i][iscal][1] << " " << vectSol[i][iscal][2] << '\n';    
            }
        }
        output_v << "      </DataArray> " << '\n';
    }
    
    




    output_v << "    </PointData>" << '\n'; 

    //WRITE ELEMENT RESULTS
    output_v << "    <CellData>" << '\n';
    
    //Some element wise result
    output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                << "Name=\"Material\" format=\"ascii\">" << '\n';
    for (int i=0; i<cmesh->NElements(); i++){
        auto compel = cmesh->ElementVec()[i];
        int matid = compel->GetWeakForm()->Id();
        output_v << matid << '\n';
    };
    output_v << "      </DataArray> " << '\n';

    for(const auto& custom_cell_field : custom_cel_fields) {
        if(custom_cell_field.data.size() != cmesh->NElements()) {
            // Custom field does not apply to all elements!
            PanicButton();
        }
        output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                << "Name=\""<< custom_cell_field.name <<"\" format=\"ascii\">" << '\n';
        for (int i=0; i<cmesh->NElements(); i++){
            output_v << custom_cell_field.data[i] << '\n';
        };
        output_v << "      </DataArray> " << '\n';
    }
    
    output_v << "    </CellData>" << '\n'; 

    //FINALIZE OUTPUT FILE
    output_v << "  </Piece>" << '\n';
    output_v << "  </UnstructuredGrid>" << '\n'
             << "</VTKFile>" << '\n';

}

void VTUGenerator::PrintResults(Arlequin *arl, std::string filename){

    //    std::cout << "Printing Velocity Results" << '\n';
    std::string s = filename + "Coarse.vtu";
    
    std::fstream output_v(s.c_str(), std::ios_base::out);

    output_v << "<?xml version=\"1.0\"?>" << '\n'
             << "<VTKFile type=\"UnstructuredGrid\">" << '\n'
             << "  <UnstructuredGrid>" << '\n'
             << "  <Piece NumberOfPoints=\"" << arl->MeshVec()[0]->NNodes()
             << "\"  NumberOfCells=\"" << arl->MeshVec()[0]->NElements()
             << "\">" << '\n';

    //WRITE NODAL COORDINATES
    output_v << "    <Points>" << '\n'
             << "      <DataArray type=\"Float64\" "
             << "NumberOfComponents=\"3\" format=\"ascii\">" << '\n';

    for (int i=0; i<arl->MeshVec()[0]->NNodes(); i++){
        auto x = arl->MeshVec()[0]->NodeVec()[i]->getCoordinates();
        output_v << x[0] << " " << x[1] << " " << x[2] << '\n';        
    };

    output_v << "      </DataArray>" << '\n'
             << "    </Points>" << '\n';
    
    //WRITE ELEMENT CONNECTIVITY
    output_v << "    <Cells>" << '\n'
             << "      <DataArray type=\"Int32\" "
             << "Name=\"connectivity\" format=\"ascii\">" << '\n';
    
    for (int i=0; i<arl->MeshVec()[0]->NElements(); i++){
        auto connec=arl->MeshVec()[0]->ElementVec()[i]->getConnectivity();
        for (int k = 0; k < connec.size(); k++)
        {
            output_v << connec[k] << " ";
        }
        output_v << '\n';
    };
    output_v << "      </DataArray>" << '\n';
  
    //WRITE OFFSETS IN DATA ARRAY
    output_v << "      <DataArray type=\"Int32\""
             << " Name=\"offsets\" format=\"ascii\">" << '\n';
    
    int aux = 0;
    for (int i=0; i<arl->MeshVec()[0]->NElements(); i++){
        aux += arl->MeshVec()[0]->ElementVec()[i]->getConnectivity().size();
        output_v << aux << '\n';
    };
    output_v << "      </DataArray>" << '\n';
  
    //WRITE ELEMENT TYPES
    output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
             << "format=\"ascii\">" << '\n';

    for (int i=0; i<arl->MeshVec()[0]->NElements(); i++){
        output_v << arl->MeshVec()[0]->ElementVec()[i]->PrintType() << '\n';
    };

    output_v << "      </DataArray>" << '\n'
             << "    </Cells>" << '\n';

    //WRITE NODAL RESULTS
    output_v << "    <PointData>" << '\n';

    output_v<< "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
        << "Name=\"" << "Signaled Distance" << "\" format=\"ascii\">" << '\n';
    for (int i=0; i<arl->MeshVec()[0]->NNodes(); i++){
        output_v << arl->GlobalNodeSignaledDistance(i) << '\n';
    }
    output_v << "      </DataArray> " << '\n';

    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
        << "Name=\"" << "Weight Function" << "\" format=\"ascii\">" << '\n';
    for (int i=0; i<arl->MeshVec()[0]->NNodes(); i++){
        output_v << arl->MeshVec()[0]->NodeVec()[i]->getWeightFunction() << '\n';
    }
    output_v << "      </DataArray> " << '\n';


    output_v << "    </PointData>" << '\n'; 

    //WRITE ELEMENT RESULTS
    output_v << "    <CellData>" << '\n';
    
    //Some element wise result
    output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                << "Name=\"Material\" format=\"ascii\">" << '\n';
    for (int i=0; i<arl->MeshVec()[0]->NElements(); i++){
        auto compel = arl->MeshVec()[0]->ElementVec()[i];
        int matid = compel->GetWeakForm()->Id();
        output_v << matid << '\n';
    };
    output_v << "      </DataArray> " << '\n';
    
    output_v << "    </CellData>" << '\n'; 

    //FINALIZE OUTPUT FILE
    output_v << "  </Piece>" << '\n';
    output_v << "  </UnstructuredGrid>" << '\n'
             << "</VTKFile>" << '\n';




    //Fine 
    //    std::cout << "Printing Velocity Results" << '\n';
    std::string sfine = filename + "Fine.vtu";
    
    std::fstream output_fine(sfine.c_str(), std::ios_base::out);

    output_fine << "<?xml version=\"1.0\"?>" << '\n'
                << "<VTKFile type=\"UnstructuredGrid\">" << '\n'
                << "  <UnstructuredGrid>" << '\n'
                << "  <Piece NumberOfPoints=\"" << arl->MeshVec()[1]->NNodes()
                << "\"  NumberOfCells=\"" << arl->MeshVec()[1]->NElements()
                << "\">" << '\n';

    //WRITE NODAL COORDINATES
    output_fine << "    <Points>" << '\n'
                << "      <DataArray type=\"Float64\" "
                << "NumberOfComponents=\"3\" format=\"ascii\">" << '\n';

    for (int i=0; i<arl->MeshVec()[1]->NNodes(); i++){
        auto x = arl->MeshVec()[1]->NodeVec()[i]->getCoordinates();
        output_fine << x[0] << " " << x[1] << " " << x[2] << '\n';        
    };

    output_fine << "      </DataArray>" << '\n'
             << "    </Points>" << '\n';
    
    //WRITE ELEMENT CONNECTIVITY
    output_fine << "    <Cells>" << '\n'
             << "      <DataArray type=\"Int32\" "
             << "Name=\"connectivity\" format=\"ascii\">" << '\n';
    
    for (int i=0; i<arl->MeshVec()[1]->NElements(); i++){
        auto connec=arl->MeshVec()[1]->ElementVec()[i]->getConnectivity();
        for (int k = 0; k < connec.size(); k++)
        {
            output_fine << connec[k] << " ";
        }
        output_fine << '\n';
    };
    output_fine << "      </DataArray>" << '\n';
  
    //WRITE OFFSETS IN DATA ARRAY
    output_fine << "      <DataArray type=\"Int32\""
             << " Name=\"offsets\" format=\"ascii\">" << '\n';
    
    aux = 0;
    for (int i=0; i<arl->MeshVec()[1]->NElements(); i++){
        aux += arl->MeshVec()[1]->ElementVec()[i]->getConnectivity().size();
        output_fine << aux << '\n';
    };
    output_fine << "      </DataArray>" << '\n';
  
    //WRITE ELEMENT TYPES
    output_fine << "      <DataArray type=\"UInt8\" Name=\"types\" "
             << "format=\"ascii\">" << '\n';

    for (int i=0; i<arl->MeshVec()[1]->NElements(); i++){
        output_fine << arl->MeshVec()[1]->ElementVec()[i]->PrintType() << '\n';
    };

    output_fine << "      </DataArray>" << '\n'
             << "    </Cells>" << '\n';

    //WRITE NODAL RESULTS
    output_fine << "    <PointData>" << '\n';

    output_fine << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
        << "Name=\"" << "Signaled Distance" << "\" format=\"ascii\">" << '\n';
    for (int i=0; i<arl->MeshVec()[1]->NNodes(); i++){
        output_fine << arl->LocalNodeSignaledDistance(i) << '\n';
    }
    output_fine << "      </DataArray> " << '\n';

    output_fine << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
        << "Name=\"" << "Weight Function" << "\" format=\"ascii\">" << '\n';
    for (int i=0; i<arl->MeshVec()[1]->NNodes(); i++){
        output_fine << arl->MeshVec()[1]->NodeVec()[i]->getWeightFunction() << '\n';
    }
    output_fine << "      </DataArray> " << '\n';


    output_fine << "    </PointData>" << '\n'; 

    //WRITE ELEMENT RESULTS
    output_fine << "    <CellData>" << '\n';
    
    //Some element wise result
    output_fine <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                << "Name=\"Material\" format=\"ascii\">" << '\n';
    for (int i=0; i<arl->MeshVec()[1]->NElements(); i++){
        auto compel = arl->MeshVec()[1]->ElementVec()[i];
        int matid = compel->GetWeakForm()->Id();
        output_fine << matid << '\n';
    };
    output_fine << "      </DataArray> " << '\n';

    output_fine <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                << "Name=\"Gluing Zone\" format=\"ascii\">" << '\n';
    for (int i=0; i<arl->MeshVec()[1]->NElements(); i++){
        int isgluing = 0;
        if (arl->fGluingElementIndex.find(i) != arl->fGluingElementIndex.end())isgluing = 1;
        output_fine << isgluing << '\n';
    };
    output_fine << "      </DataArray> " << '\n';
    
    output_fine << "    </CellData>" << '\n'; 

    //FINALIZE OUTPUT FILE
    output_fine << "  </Piece>" << '\n';
    output_fine << "  </UnstructuredGrid>" << '\n'
             << "</VTKFile>" << '\n';


}