#include "VTUGenerator.h"


void VTUGenerator::PrintResultsGraph(CompMesh *cmesh, std::string filename, std::vector<std::string> &scalnames, std::vector<std::string> &vecnames, int step){

    auto graphmesh = cmesh->GetGraphMesh();

    //    std::cout << "Printing Velocity Results" << std::endl;
    std::string s = filename + std::to_string(step) + ".vtu";
    
    std::fstream output_v(s.c_str(), std::ios_base::out);

    output_v << "<?xml version=\"1.0\"?>" << std::endl
             << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
             << "  <UnstructuredGrid>" << std::endl
             << "  <Piece NumberOfPoints=\"" << graphmesh->NNodes()
             << "\"  NumberOfCells=\"" << graphmesh->NElements()
             << "\">" << std::endl;

    //WRITE NODAL COORDINATES
    output_v << "    <Points>" << std::endl
             << "      <DataArray type=\"Float64\" "
             << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

    for (int i=0; i<graphmesh->NNodes(); i++){
        auto x = graphmesh->Node(i);
        output_v << x[0] << " " << x[1] << " " << x[2] << std::endl;        
    };

    output_v << "      </DataArray>" << std::endl
             << "    </Points>" << std::endl;
    
    //WRITE ELEMENT CONNECTIVITY
    output_v << "    <Cells>" << std::endl
             << "      <DataArray type=\"Int32\" "
             << "Name=\"connectivity\" format=\"ascii\">" << std::endl;
    
    for (int i=0; i<graphmesh->NElements(); i++){
        auto connec=graphmesh->Connect(i);
        for (int k = 0; k < connec.size(); k++)
        {
            output_v << connec[k] << " ";
        }
        output_v << std::endl;
    };
    output_v << "      </DataArray>" << std::endl;
  
    //WRITE OFFSETS IN DATA ARRAY
    output_v << "      <DataArray type=\"Int32\""
             << " Name=\"offsets\" format=\"ascii\">" << std::endl;
    
    int aux = 0;
    for (int i=0; i<graphmesh->NElements(); i++){
        aux += graphmesh->Connect(i).size();
        output_v << aux << std::endl;
    };
    output_v << "      </DataArray>" << std::endl;
  
    //WRITE ELEMENT TYPES
    output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
             << "format=\"ascii\">" << std::endl;

    for (int i=0; i<graphmesh->NElements(); i++){
        output_v << graphmesh->ElType(i) << std::endl;
    };

    output_v << "      </DataArray>" << std::endl
             << "    </Cells>" << std::endl;

    //WRITE NODAL RESULTS
    output_v << "    <PointData>" << std::endl;

    //Scalar values
    std::map<int64_t,std::vector<VecDouble>> scalSol;
    std::map<int64_t,std::vector<VecDouble>> vectSol;
    // std::vector<std::vector<VecDouble,scalnames.size()>,graphmesh->NNodes()> scalSol;
    for (int64_t iel = 0; iel < graphmesh->NElements(); iel++){
        auto compel = cmesh->ElementVec()[iel];
        auto graphconnect = graphmesh->Connect(iel);
        if (compel->Dimension() != cmesh->Dimension()) continue;

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
            << "Name=\"" << scalnames[iscal] << "\" format=\"ascii\">" << std::endl;
        for (int i=0; i<graphmesh->NNodes(); i++){
            // if (scalSol[i].size()==0){
                // output_v <<std::scientific<< "0" << std::endl;
            // }else{
                output_v <<std::scientific<< scalSol[i][iscal][0] << std::endl;
            // }
        }
        output_v << "      </DataArray> " << std::endl;
    }
    
    for (int iscal = 0; iscal < vecnames.size(); iscal++){
        output_v<< "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
            << "Name=\"" << vecnames[iscal] << "\" format=\"ascii\">" << std::endl;
        for (int i=0; i<graphmesh->NNodes(); i++){
            // if (vectSol[i].size()==0){
                // output_v <<std::scientific<< "0,0,0" << std::endl;
            // }else{
                output_v<<std::scientific << vectSol[i][iscal][0] << " " << vectSol[i][iscal][1] << " " << vectSol[i][iscal][2] << std::endl;    
            // }
        }
        output_v << "      </DataArray> " << std::endl;
    }
    
    




    output_v << "    </PointData>" << std::endl; 

    //WRITE ELEMENT RESULTS
    output_v << "    <CellData>" << std::endl;
    
    //Some element wise result
    output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                << "Name=\"Material\" format=\"ascii\">" << std::endl;
    for (int i=0; i<graphmesh->NElements(); i++){
        auto compel = cmesh->ElementVec()[i];
        int matid = compel->GetWeakForm()->Id();
        output_v << matid << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;
    
    output_v << "    </CellData>" << std::endl; 

    //FINALIZE OUTPUT FILE
    output_v << "  </Piece>" << std::endl;
    output_v << "  </UnstructuredGrid>" << std::endl
             << "</VTKFile>" << std::endl;

}

void VTUGenerator::PrintResults(CompMesh *cmesh, std::string filename, std::vector<std::string> &scalnames, std::vector<std::string> &vecnames, const std::vector<CustomCellField>& custom_cel_fields, int step){

    //    std::cout << "Printing Velocity Results" << std::endl;
    std::string s = filename + std::to_string(step) + ".vtu";
    
    std::fstream output_v(s.c_str(), std::ios_base::out);

    output_v << "<?xml version=\"1.0\"?>" << std::endl
             << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
             << "  <UnstructuredGrid>" << std::endl
             << "  <Piece NumberOfPoints=\"" << cmesh->NNodes()
             << "\"  NumberOfCells=\"" << cmesh->NElements()
             << "\">" << std::endl;

    //WRITE NODAL COORDINATES
    output_v << "    <Points>" << std::endl
             << "      <DataArray type=\"Float64\" "
             << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

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
        output_v << std::endl;
    };

    output_v << "      </DataArray>" << std::endl
             << "    </Points>" << std::endl;
    
    //WRITE ELEMENT CONNECTIVITY
    output_v << "    <Cells>" << std::endl
             << "      <DataArray type=\"Int32\" "
             << "Name=\"connectivity\" format=\"ascii\">" << std::endl;
    
    for (int i=0; i<cmesh->NElements(); i++){
        auto connec=cmesh->ElementVec()[i]->getConnectivity();
        for (int k = 0; k < connec.size(); k++)
        {
            output_v << connec[k] << " ";
        }
        output_v << std::endl;
    };
    output_v << "      </DataArray>" << std::endl;
  
    //WRITE OFFSETS IN DATA ARRAY
    output_v << "      <DataArray type=\"Int32\""
             << " Name=\"offsets\" format=\"ascii\">" << std::endl;
    
    int aux = 0;
    for (int i=0; i<cmesh->NElements(); i++){
        aux += cmesh->ElementVec()[i]->getConnectivity().size();
        output_v << aux << std::endl;
    };
    output_v << "      </DataArray>" << std::endl;
  
    //WRITE ELEMENT TYPES
    output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
             << "format=\"ascii\">" << std::endl;

    for (int i=0; i<cmesh->NElements(); i++){
        output_v << cmesh->ElementVec()[i]->PrintType() << std::endl;
    };

    output_v << "      </DataArray>" << std::endl
             << "    </Cells>" << std::endl;

    //WRITE NODAL RESULTS
    output_v << "    <PointData>" << std::endl;

    //Scalar values
    std::map<int64_t,std::vector<VecDouble>> scalSol;
    std::map<int64_t,std::vector<VecDouble>> vectSol;
    // std::vector<std::vector<VecDouble,scalnames.size()>,graphmesh->NNodes()> scalSol;
    for (int64_t iel = 0; iel < cmesh->NElements(); iel++){
        auto compel = cmesh->ElementVec()[iel];
        auto graphconnect = cmesh->ElementVec()[iel]->getConnectivity();
        if (compel->Dimension() != cmesh->Dimension()) continue;

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
            << "Name=\"" << scalnames[iscal] << "\" format=\"ascii\">" << std::endl;
        for (int i=0; i<cmesh->NNodes(); i++){
            // if (scalSol[i].size()==0){
                // output_v <<std::scientific<< "0" << std::endl;
            // }else{
                output_v <<std::scientific<< scalSol[i][iscal][0] << std::endl;
            // }
        }
        output_v << "      </DataArray> " << std::endl;
    }
    
    for (int iscal = 0; iscal < vecnames.size(); iscal++){
        output_v<< "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
            << "Name=\"" << vecnames[iscal] << "\" format=\"ascii\">" << std::endl;
        for (int i=0; i<cmesh->NNodes(); i++){
            // if (vectSol[i].size()==0){
                // output_v <<std::scientific<< "0,0,0" << std::endl;
            // }else{
                output_v<<std::scientific << vectSol[i][iscal][0] << " " << vectSol[i][iscal][1] << " " << vectSol[i][iscal][2] << std::endl;    
            // }
        }
        output_v << "      </DataArray> " << std::endl;
    }
    
    




    output_v << "    </PointData>" << std::endl; 

    //WRITE ELEMENT RESULTS
    output_v << "    <CellData>" << std::endl;
    
    //Some element wise result
    output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                << "Name=\"Material\" format=\"ascii\">" << std::endl;
    for (int i=0; i<cmesh->NElements(); i++){
        auto compel = cmesh->ElementVec()[i];
        int matid = compel->GetWeakForm()->Id();
        output_v << matid << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;

    for(const auto& custom_cell_field : custom_cel_fields) {
        if(custom_cell_field.data.size() != cmesh->NElements()) {
            // Custom field does not apply to all elements!
            PanicButton();
        }
        output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                << "Name=\""<< custom_cell_field.name <<"\" format=\"ascii\">" << std::endl;
        for (int i=0; i<cmesh->NElements(); i++){
            output_v << custom_cell_field.data[i] << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    }
    
    output_v << "    </CellData>" << std::endl; 

    //FINALIZE OUTPUT FILE
    output_v << "  </Piece>" << std::endl;
    output_v << "  </UnstructuredGrid>" << std::endl
             << "</VTKFile>" << std::endl;

}

void VTUGenerator::PrintResults(Arlequin *arl, std::string filename){

    //    std::cout << "Printing Velocity Results" << std::endl;
    std::string s = filename + "Coarse.vtu";
    
    std::fstream output_v(s.c_str(), std::ios_base::out);

    output_v << "<?xml version=\"1.0\"?>" << std::endl
             << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
             << "  <UnstructuredGrid>" << std::endl
             << "  <Piece NumberOfPoints=\"" << arl->MeshVec()[0]->NNodes()
             << "\"  NumberOfCells=\"" << arl->MeshVec()[0]->NElements()
             << "\">" << std::endl;

    //WRITE NODAL COORDINATES
    output_v << "    <Points>" << std::endl
             << "      <DataArray type=\"Float64\" "
             << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

    for (int i=0; i<arl->MeshVec()[0]->NNodes(); i++){
        auto x = arl->MeshVec()[0]->NodeVec()[i]->getCoordinates();
        output_v << x[0] << " " << x[1] << " " << x[2] << std::endl;        
    };

    output_v << "      </DataArray>" << std::endl
             << "    </Points>" << std::endl;
    
    //WRITE ELEMENT CONNECTIVITY
    output_v << "    <Cells>" << std::endl
             << "      <DataArray type=\"Int32\" "
             << "Name=\"connectivity\" format=\"ascii\">" << std::endl;
    
    for (int i=0; i<arl->MeshVec()[0]->NElements(); i++){
        auto connec=arl->MeshVec()[0]->ElementVec()[i]->getConnectivity();
        for (int k = 0; k < connec.size(); k++)
        {
            output_v << connec[k] << " ";
        }
        output_v << std::endl;
    };
    output_v << "      </DataArray>" << std::endl;
  
    //WRITE OFFSETS IN DATA ARRAY
    output_v << "      <DataArray type=\"Int32\""
             << " Name=\"offsets\" format=\"ascii\">" << std::endl;
    
    int aux = 0;
    for (int i=0; i<arl->MeshVec()[0]->NElements(); i++){
        aux += arl->MeshVec()[0]->ElementVec()[i]->getConnectivity().size();
        output_v << aux << std::endl;
    };
    output_v << "      </DataArray>" << std::endl;
  
    //WRITE ELEMENT TYPES
    output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
             << "format=\"ascii\">" << std::endl;

    for (int i=0; i<arl->MeshVec()[0]->NElements(); i++){
        output_v << arl->MeshVec()[0]->ElementVec()[i]->PrintType() << std::endl;
    };

    output_v << "      </DataArray>" << std::endl
             << "    </Cells>" << std::endl;

    //WRITE NODAL RESULTS
    output_v << "    <PointData>" << std::endl;

    output_v<< "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
        << "Name=\"" << "Signaled Distance" << "\" format=\"ascii\">" << std::endl;
    for (int i=0; i<arl->MeshVec()[0]->NNodes(); i++){
        output_v << arl->GlobalNodeSignaledDistance(i) << std::endl;
    }
    output_v << "      </DataArray> " << std::endl;

    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
        << "Name=\"" << "Weight Function" << "\" format=\"ascii\">" << std::endl;
    for (int i=0; i<arl->MeshVec()[0]->NNodes(); i++){
        output_v << arl->MeshVec()[0]->NodeVec()[i]->getWeightFunction() << std::endl;
    }
    output_v << "      </DataArray> " << std::endl;


    output_v << "    </PointData>" << std::endl; 

    //WRITE ELEMENT RESULTS
    output_v << "    <CellData>" << std::endl;
    
    //Some element wise result
    output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                << "Name=\"Material\" format=\"ascii\">" << std::endl;
    for (int i=0; i<arl->MeshVec()[0]->NElements(); i++){
        auto compel = arl->MeshVec()[0]->ElementVec()[i];
        int matid = compel->GetWeakForm()->Id();
        output_v << matid << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;
    
    output_v << "    </CellData>" << std::endl; 

    //FINALIZE OUTPUT FILE
    output_v << "  </Piece>" << std::endl;
    output_v << "  </UnstructuredGrid>" << std::endl
             << "</VTKFile>" << std::endl;




    //Fine 
    //    std::cout << "Printing Velocity Results" << std::endl;
    std::string sfine = filename + "Fine.vtu";
    
    std::fstream output_fine(sfine.c_str(), std::ios_base::out);

    output_fine << "<?xml version=\"1.0\"?>" << std::endl
                << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                << "  <UnstructuredGrid>" << std::endl
                << "  <Piece NumberOfPoints=\"" << arl->MeshVec()[1]->NNodes()
                << "\"  NumberOfCells=\"" << arl->MeshVec()[1]->NElements()
                << "\">" << std::endl;

    //WRITE NODAL COORDINATES
    output_fine << "    <Points>" << std::endl
                << "      <DataArray type=\"Float64\" "
                << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

    for (int i=0; i<arl->MeshVec()[1]->NNodes(); i++){
        auto x = arl->MeshVec()[1]->NodeVec()[i]->getCoordinates();
        output_fine << x[0] << " " << x[1] << " " << x[2] << std::endl;        
    };

    output_fine << "      </DataArray>" << std::endl
             << "    </Points>" << std::endl;
    
    //WRITE ELEMENT CONNECTIVITY
    output_fine << "    <Cells>" << std::endl
             << "      <DataArray type=\"Int32\" "
             << "Name=\"connectivity\" format=\"ascii\">" << std::endl;
    
    for (int i=0; i<arl->MeshVec()[1]->NElements(); i++){
        auto connec=arl->MeshVec()[1]->ElementVec()[i]->getConnectivity();
        for (int k = 0; k < connec.size(); k++)
        {
            output_fine << connec[k] << " ";
        }
        output_fine << std::endl;
    };
    output_fine << "      </DataArray>" << std::endl;
  
    //WRITE OFFSETS IN DATA ARRAY
    output_fine << "      <DataArray type=\"Int32\""
             << " Name=\"offsets\" format=\"ascii\">" << std::endl;
    
    aux = 0;
    for (int i=0; i<arl->MeshVec()[1]->NElements(); i++){
        aux += arl->MeshVec()[1]->ElementVec()[i]->getConnectivity().size();
        output_fine << aux << std::endl;
    };
    output_fine << "      </DataArray>" << std::endl;
  
    //WRITE ELEMENT TYPES
    output_fine << "      <DataArray type=\"UInt8\" Name=\"types\" "
             << "format=\"ascii\">" << std::endl;

    for (int i=0; i<arl->MeshVec()[1]->NElements(); i++){
        output_fine << arl->MeshVec()[1]->ElementVec()[i]->PrintType() << std::endl;
    };

    output_fine << "      </DataArray>" << std::endl
             << "    </Cells>" << std::endl;

    //WRITE NODAL RESULTS
    output_fine << "    <PointData>" << std::endl;

    output_fine << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
        << "Name=\"" << "Signaled Distance" << "\" format=\"ascii\">" << std::endl;
    for (int i=0; i<arl->MeshVec()[1]->NNodes(); i++){
        output_fine << arl->LocalNodeSignaledDistance(i) << std::endl;
    }
    output_fine << "      </DataArray> " << std::endl;

    output_fine << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
        << "Name=\"" << "Weight Function" << "\" format=\"ascii\">" << std::endl;
    for (int i=0; i<arl->MeshVec()[1]->NNodes(); i++){
        output_fine << arl->MeshVec()[1]->NodeVec()[i]->getWeightFunction() << std::endl;
    }
    output_fine << "      </DataArray> " << std::endl;


    output_fine << "    </PointData>" << std::endl; 

    //WRITE ELEMENT RESULTS
    output_fine << "    <CellData>" << std::endl;
    
    //Some element wise result
    output_fine <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                << "Name=\"Material\" format=\"ascii\">" << std::endl;
    for (int i=0; i<arl->MeshVec()[1]->NElements(); i++){
        auto compel = arl->MeshVec()[1]->ElementVec()[i];
        int matid = compel->GetWeakForm()->Id();
        output_fine << matid << std::endl;
    };
    output_fine << "      </DataArray> " << std::endl;

    output_fine <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                << "Name=\"Gluing Zone\" format=\"ascii\">" << std::endl;
    for (int i=0; i<arl->MeshVec()[1]->NElements(); i++){
        int isgluing = 0;
        if (arl->fGluingElementIndex.find(i) != arl->fGluingElementIndex.end())isgluing = 1;
        output_fine << isgluing << std::endl;
    };
    output_fine << "      </DataArray> " << std::endl;
    
    output_fine << "    </CellData>" << std::endl; 

    //FINALIZE OUTPUT FILE
    output_fine << "  </Piece>" << std::endl;
    output_fine << "  </UnstructuredGrid>" << std::endl
             << "</VTKFile>" << std::endl;


}