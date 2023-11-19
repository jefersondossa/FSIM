#include "VTUGenerator.h"


void VTUGenerator::PrintResults(CompMesh *cmesh, std::string filename, std::vector<std::string> &scalnames, std::vector<std::string> &vecnames, int step){

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

        for (int inode = 0; inode < compel->NElNodes(); inode++){
            auto xparametric = compel->NodeCoord(inode);
            compel->IntegrationData().fAdimCoord = xparametric;
            if (!compel->IntegrationData().fNeedsSol || compel->IntegrationData().fNeedsDSol){
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
            output_v << scalSol[i][iscal][0] << std::endl;
        }
        output_v << "      </DataArray> " << std::endl;
    }
    
    for (int iscal = 0; iscal < vecnames.size(); iscal++){
        output_v<< "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
            << "Name=\"" << vecnames[iscal] << "\" format=\"ascii\">" << std::endl;
        for (int i=0; i<graphmesh->NNodes(); i++){
            output_v << vectSol[i][iscal][0] << " " << vectSol[i][iscal][1] << " " << vectSol[i][iscal][2] << std::endl;
        }
        output_v << "      </DataArray> " << std::endl;
    }
    
    




    output_v << "    </PointData>" << std::endl; 

    //WRITE ELEMENT RESULTS
    output_v << "    <CellData>" << std::endl;
    
    //Some element wise result
    // output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
    //             << "Name=\"Process\" format=\"ascii\">" << std::endl;
    // for (int i=0; i<NElements(); i++){
    //     output_v << part_elem[i] << std::endl;
    // };
    // output_v << "      </DataArray> " << std::endl;
    
    output_v << "    </CellData>" << std::endl; 

    //FINALIZE OUTPUT FILE
    output_v << "  </Piece>" << std::endl;
    output_v << "  </UnstructuredGrid>" << std::endl
             << "</VTKFile>" << std::endl;

}