#include "VTUGenerator.h"


void VTUGenerator::PrintResults(CompMesh *cmesh, std::string filename, int step){

    //    std::cout << "Printing Velocity Results" << std::endl;
    std::string s = filename + std::to_string(step) + ".vtu";
    
    std::fstream output_v(s.c_str(), std::ios_base::out);
    int DIM = cmesh->Dimension();
    int DEG = cmesh->GetDefaultOrder();

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
        output_v << x[0] << " " << x[1] << " " << x[2] << std::endl;        
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
        output_v << aux +cmesh-> NElNodes() << std::endl;
        aux += cmesh->NElNodes();
    };
    output_v << "      </DataArray>" << std::endl;
  
    //WRITE ELEMENT TYPES
    output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
             << "format=\"ascii\">" << std::endl;
    
    int val = 0;
    if (DEG == 1)val = 5;
    if (DEG == 2)val = 22;
    if (DEG == 3)val = 69;

    for (int i=0; i<cmesh->NElements(); i++){
        output_v << val << std::endl;
    };

    output_v << "      </DataArray>" << std::endl
             << "    </Cells>" << std::endl;

    //WRITE NODAL RESULTS
    output_v << "    <PointData>" << std::endl;
    if (cmesh->getProblemParameters().ProbType()==EPoisson){
        output_v<< "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
            << "Name=\"Solution\" format=\"ascii\">" << std::endl;
        for (int i=0; i<cmesh->NNodes(); i++){
            output_v << cmesh->NodeVec()[i] -> GetSolution(0) << std::endl;
        }
        output_v << "      </DataArray> " << std::endl;
    } else if (cmesh->getProblemParameters().ProbType()==EElastic || cmesh->getProblemParameters().ProbType() == ESolidPositional){
        output_v<< "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
            << "Name=\"Displacement\" format=\"ascii\">" << std::endl;
        for (int i=0; i<cmesh->NNodes(); i++){           
            output_v << cmesh->NodeVec()[i] -> GetSolution(0) << " "             
                    << cmesh->NodeVec()[i] -> GetSolution(1) << " ";
            if (DIM == 2) {
                output_v <<  0. << std::endl;
            } else {
                output_v <<  cmesh->NodeVec()[i] -> GetSolution(2) << std::endl;
            }
        }
        output_v << "      </DataArray> " << std::endl;
    } else if (cmesh->getProblemParameters().ProbType()==EStokes || cmesh->getProblemParameters().ProbType()==ENavierStokes){
        output_v<< "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
            << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
        for (int i=0; i<cmesh->NNodes(); i++){
            // std::cout << "Solution - " << NodeVec()[i] -> GetSolution(0) << " " << NodeVec()[i] -> GetSolution(0) << std::endl;
            output_v << cmesh->NodeVec()[i] -> GetSolution(0) << " "             
                    << cmesh->NodeVec()[i] -> GetSolution(1) << " ";
            if (DIM == 2) {
                output_v <<  0. << std::endl;
            } else {
                output_v <<  cmesh->NodeVec()[i] -> GetSolution(2) << std::endl;
            }
        }
        output_v << "      </DataArray> " << std::endl;
        output_v<< "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
            << "Name=\"Pressure\" format=\"ascii\">" << std::endl;
        for (int i=0; i<cmesh->NNodes(); i++){
            output_v << cmesh->NodeVec()[i] -> GetSolution(DIM) << std::endl;
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