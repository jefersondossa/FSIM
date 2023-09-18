#include "Fluid.h"
#include "petscpartitioner.h"
#include "metis.h"

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------------READS THE .TXT INPUT FILE--------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Fluid<DIM,DEG>::readInputFile(const std::string& inputFile, std::ofstream& mirrorData){

    std::ifstream inputData(inputFile.c_str());
    std::string line;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);      
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    velocityInf.resize(3);
    fieldForces.resize(3);

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++READING PROBLEM VARIABLES+++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
    //Read number of nodes, elements, time steps and printing frequence
    inputData >> numTimeSteps >> printFreq;
    mirrorData << "Number of Time Steps   = " << numTimeSteps << std::endl;
    mirrorData << "Printing Frequence     = " << printFreq << std::endl;
    
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
    //Read undisturbed velocity and pressure components
    inputData >> velocityInf[0] >> velocityInf[1] >> velocityInf[2] >> pressInf;

    mirrorData << "Undisturbed Velocity x = " << velocityInf[0] << std::endl;
    mirrorData << "Undisturbed Velocity y = " << velocityInf[1] << std::endl;
    mirrorData << "Undisturbed Velocity z = " << velocityInf[2] << std::endl;
    mirrorData << "Undisturbed Pressure   = " << pressInf << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read undisturbed density, temperature, viscosity and thermal condutivity
    inputData >> rhoInf >> viscInf;

    mirrorData << "Undisturbed Density    = " << rhoInf << std::endl;
    mirrorData << "Undisturbed Viscosity  = " << viscInf << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read time step lenght
    inputData >> dTime >> integScheme;

    mirrorData << "Time Step              = " << dTime << std::endl;
    mirrorData << "Time Integration Scheme= " << integScheme << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);getline(inputData,line);

    //Read field forces
    inputData >> fieldForces[0] >> fieldForces[1] >> fieldForces[2];

    mirrorData << "Field Forces x         = " << fieldForces[0] << std::endl;
    mirrorData << "Field Forces y         = " << fieldForces[1] << std::endl;
    mirrorData << "Field Forces z         = " << fieldForces[2] << std::endl \
               << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);getline(inputData,line);




    //Read Arlequin variables
    double k1,k2;
    inputData >> glueZoneThickness >> arlequinEpsilon >> weightFunctionBehavior
              >> k1 >> k2;

    mirrorData << "Glue Zone Thickness    = " << glueZoneThickness << std::endl;
    mirrorData << "Epsilon                = " << arlequinEpsilon << std::endl;
    mirrorData << "Energy Weight Function = " << weightFunctionBehavior
               << std::endl;
    mirrorData << "K1                     = " << k1 << std::endl;
    mirrorData << "K2                     = " << k2 << std::endl << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    fluidParameters.setViscosity(viscInf);
    fluidParameters.setDensity(rhoInf);
    fluidParameters.setTimeStep(dTime);
    fluidParameters.setSpectralRadius(integScheme);
    fluidParameters.setFieldForce(fieldForces);
    fluidParameters.setArlequinOperatorConstants(k1,k2);
    fluidParameters.setVelocityInf(velocityInf);

    //Drag and lift
    inputData >> computeDragAndLift >> numberOfLines; 
    dragAndLiftBoundary.reserve(numberOfLines);
    for (int i = 0; i < numberOfLines; ++i)
    {
        int aux;
        inputData >> aux; 
        dragAndLiftBoundary.push_back(aux);
    }
    

    mirrorData << "Compute Drag and Lift  = " << computeDragAndLift<< std::endl;
    mirrorData << "Number of Lines  = " << numberOfLines << std::endl;
    for (int i = 0; i < numberOfLines; ++i)
    {
        mirrorData << "Lines  = " << dragAndLiftBoundary[i] << std::endl;
    }

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Printing results
    inputData >> printVelocity;              getline(inputData,line);
    inputData >> printAcceleration;              getline(inputData,line);
    inputData >> printRealVelocity;          getline(inputData,line);
    inputData >> printLagrangeMultipliers;   getline(inputData,line);
    inputData >> printElementCorrespondence; getline(inputData,line);
    inputData >> printDistFunction;          getline(inputData,line);
    inputData >> printEnergyWeightFunction;  getline(inputData,line);
    inputData >> printPressure;              getline(inputData,line);
    inputData >> printRealPressure;          getline(inputData,line);
    inputData >> printVorticity;             getline(inputData,line);
    inputData >> printInnerNormal;           getline(inputData,line);
    inputData >> printMeshVelocity;          getline(inputData,line);
    inputData >> printMeshDisplacement;      getline(inputData,line);
    inputData >> printGlueZone;              getline(inputData,line);
    inputData >> printJacobian;              getline(inputData,line);
    inputData >> printProcess;               getline(inputData,line);
    inputData >> printLines;              

    mirrorData << "PrintVelocity              = " << printVelocity << std::endl;
    mirrorData << "printAcceleration          = " << printAcceleration << std::endl;
    mirrorData << "PrintRealVelocity          = " << printRealVelocity << std::endl;
    mirrorData << "PrintLagrangeMultipliers   = " << printLagrangeMultipliers << std::endl;
    mirrorData << "PrintElementCorrespondence = " << printElementCorrespondence << std::endl;
    mirrorData << "PrintDistFunction          = " << printDistFunction << std::endl;
    mirrorData << "PrintEnergyWeightFunction  = " << printEnergyWeightFunction << std::endl;
    mirrorData << "PrintPressure              = " << printPressure << std::endl;
    mirrorData << "PrintRealPressure          = " << printRealPressure << std::endl;
    mirrorData << "PrintVorticity             = " << printVorticity << std::endl;
    mirrorData << "PrintInnerNormal           = " << printInnerNormal << std::endl;
    mirrorData << "PrintMeshVelocity          = " << printMeshVelocity << std::endl;
    mirrorData << "PrintMeshDisplacement      = " << printMeshDisplacement << std::endl;
    mirrorData << "PrintGlueZone              = " << printGlueZone << std::endl;
    mirrorData << "PrintJacobian              = " << printJacobian << std::endl;
    mirrorData << "PrintProcess               = " << printProcess << std::endl;
    mirrorData << "PrintLines                 = " << printLines << std::endl << std::endl;

    return;
}

//------------------------------------------------------------------------------
//------------------------------READS THE MESH NODES----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Fluid<DIM,DEG>::readNodes(std::ifstream &file, std::ofstream& mirrorData){

    if (rank == 0) std::cout << "2/9 Reading nodes..." << std::endl;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++NODES++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    std::string line;
    file >> numNodes;
    nodes_.reserve(numNodes);
    std::getline(file, line);
    int index = 0;
    if (rank == 0) std::cout << "Number of Nodes " << " " << numNodes << std::endl;
    for (int i = 0; i < numNodes; i++){
        VecDouble x(DIM);
        std::getline(file, line);
        std::vector<std::string> tokens = split2(line, " ");
        
        for (int j = 0; j < DIM; j++) std::istringstream(tokens[j+1]) >> x[j];
        
        Node<DIM,DEG> *node = new Node<DIM,DEG>(x,index);
        nodes_.push_back(node);
        index++;
    }
    std::getline(file, line); std::getline(file, line);

    mirrorData << "Nodal Coordinates " << numNodes << std::endl;
    for (int i = 0 ; i<numNodes; i++){
        VecDouble x = nodes_[i]->getCoordinates();       
        for (int j=0; j<DIM; j++){
            mirrorData << x[j] << " ";
        };
        mirrorData << std::endl;
        nodes_[i] -> setVelocity(fluidParameters.getVelocityInf());
        nodes_[i] -> setPreviousVelocity(fluidParameters.getVelocityInf());
        double zero = 0.;
        for (int k=0; k<DIM; k++){
            nodes_[i] -> setMeshVelocityComponent(k,zero);
            nodes_[i] -> setPreviousMeshVelocityComponent(k,zero);
            nodes_[i] -> setPreviousCoordinates(k,x[k]);
        };
    };
    return;
}

//------------------------------------------------------------------------------
//---------------------------READS THE MESH ELEMENTS----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Fluid<DIM,DEG>::readElements(Geometry* &geometry_, std::ifstream &file, std::ofstream& mirrorData, std::vector<Element<DIM,DEG>*> &elementsAux_, std::unordered_map<int, std::string> &physicalEntities){

    if (rank == 0) std::cout << "3/9 Reading elements..." << std::endl;

    //defyning the maps that are used to store the elements information
    std::unordered_map<int, std::string> gmshElement = { {1, "line"}, {2, "triangle"}, {3, "quadrilateral"}, {8, "line3"}, {9, "triangle6"}, {10, "quadrilateral9"}, {15, "vertex"}, {16, "quadrilateral8"}, {20, "triangle9"}, {21, "triangle10"}, {26, "line4"}, {36, "quadrilateral16"}, {39, "quadrilateral12"} };
    std::unordered_map<std::string, int> numNodes2 = { {"vertex", 1}, {"line", 2}, {"triangle", 3}, {"quadrilateral", 4}, {"line3", 3}, {"triangle6", 6}, {"quadrilateral8", 8}, {"quadrilateral9", 9}, {"line4", 4}, {"triangle", 9}, {"triangle10", 10}, {"quadrilateral12", 12}, {"quadrilateral16", 16}};
    std::unordered_map<std::string, std::string> supportedElements = { {"triangle", "T3"}, {"triangle6", "T6"}, {"triangle10", "T10"}, {"quadrilateral", "Q4"}, {"quadrilateral8", "Q8"}, {"quadrilateral9", "Q9"}, {"quadrilateral12", "Q12"}, {"quadrilateral16", "Q16"}, {"tetrahedron4", "TET4"}, {"tetrahedron10", "TET10"}, {"tetrahedron20", "TET20"} };
    std::unordered_map<Line*, std::vector< std::vector<int> >> lineElements;

    std::string line;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++ELEMENTS++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int number_elements;
    file >> number_elements;
    //elements_.reserve(number_elements);
    elementsAux_.reserve(number_elements);

    boundary_.reserve(number_elements/10);
    int index = 0;
    std::getline(file, line);
    int cont = 0;

    numBoundElems = 0;
    numElem = 0;
    numFSIInterfaces = 0;

    std::vector<BoundaryCondition*> dirichlet, neumann, glue, FSinterface;
    dirichlet = geometry_->getBoundaryCondition("DIRICHLET"); 
    neumann = geometry_->getBoundaryCondition("NEUMANN"); 
    glue = geometry_->getBoundaryCondition("GLUE");
    FSinterface = geometry_->getBoundaryCondition("FSINTERFACE");

    numFSIInterfaces = FSinterface.size();

    for (int i = 0; i < number_elements; i++)
    {
        std::getline(file, line);
        std::vector<std::string> tokens = split2(line, " ");
        std::vector<int> values(tokens.size(), 0);
        for (size_t j = 0; j < tokens.size(); j++)
            std::istringstream(tokens[j]) >> values[j];
        std::string elementType = gmshElement[values[1]];
        int number_nodes_per_element = numNodes2[elementType];
        std::vector<int> elementNodes;
        elementNodes.reserve(number_nodes_per_element);

        for (size_t j = 5 ; j < values.size(); j++)
            elementNodes.push_back(values[j]-1);
 
        std::string name = physicalEntities[values[3]];
        //Adding domain elements
        if (name[0] == 'v'){
            // if(rank == 0){
                Volume* object = geometry_ -> getVolume(name);
                numElem++;

                VecInt connect(nElNodes);

                if (DEG == 2){
                    connect[2] = elementNodes[2];
                    connect[1] = elementNodes[1];
                    connect[3] = elementNodes[3];
                    connect[0] = elementNodes[0];
                    connect[5] = elementNodes[5];
                    connect[9] = elementNodes[8];
                    connect[6] = elementNodes[6];
                    connect[7] = elementNodes[7];
                    connect[4] = elementNodes[4];
                    connect[8] = elementNodes[9];
                } else {
                    for (int k = 0; k < nElNodes; k++) connect[k] = elementNodes[k];
                }

                Element<DIM,DEG> *el = new Element<DIM,DEG>(index++,connect,this);
                elementsAux_.push_back(el);

                for (int k = 0; k < nElNodes; k++){
                    nodes_[connect[k]] -> pushInverseIncidence(index);
                };
            // }
        }
        else if (name[0] == 's') {
            if (DIM == 3){
                VecInt connectB(nBdNodes);

                for (int i = 0; i < nBdNodes; i++) connectB[i] = elementNodes[i];

                int ibound;

                std::string::size_type sz;   // alias of size_t
                ibound = std::stoi (&name[1],nullptr,10);

                VecInt constrain(3);
                VecDouble value(3);

                for (int i = 0; i < dirichlet.size(); i++){
                    if (name == dirichlet[i] -> getLineName()){
                        if ((dirichlet[i] -> getComponentX()).size() == 0){
                            constrain[0] = 0; value[0] = 0;
                        }else{
                            std::vector<double> c = dirichlet[i] -> getComponentX();
                            constrain[0] = 1;
                            value[0] = c[0];
                        }
                        if ((dirichlet[i] -> getComponentY()).size() == 0){
                            constrain[1] = 0; value[1] = 0;
                        }else{
                            std::vector<double> c = dirichlet[i] -> getComponentY();
                            constrain[1] = 1;
                            value[1] = c[0];
                        }
                        if ((dirichlet[i] -> getComponentZ()).size() == 0){
                            constrain[2] = 0; value[2] = 0;
                        }else{
                            std::vector<double> c = dirichlet[i] -> getComponentZ();
                            constrain[2] = 1;
                            value[2] = c[0];
                        }
                    }
                }
                for (int i = 0; i < neumann.size(); i++){
                    if (name == neumann[i] -> getLineName()){
                        if ((neumann[i] -> getComponentX()).size() == 0){
                            constrain[0] = 0; value[0] = 0;
                        }else{
                            std::vector<double> c = neumann[i] -> getComponentX();
                            constrain[0] = 0;
                            value[0] = c[0];
                        }
                        if ((neumann[i] -> getComponentY()).size() == 0){
                            constrain[1] = 0; value[1] = 0;
                        }else{
                            std::vector<double> c = neumann[i] -> getComponentY();
                            constrain[1] = 0;
                            value[1] = c[0];
                        }
                        if ((neumann[i] -> getComponentZ()).size() == 0){
                            constrain[2] = 0; value[2] = 0;
                        }else{
                            std::vector<double> c = neumann[i] -> getComponentZ();
                            constrain[2] = 0;
                            value[2] = c[0];
                        }
                    }
                }  
                for (int i = 0; i < glue.size(); i++){
                    if (name == glue[i] -> getLineName()){
                        if ((glue[i] -> getComponentX()).size() == 0){
                            constrain[0] = 2; value[0] = 0;
                        }else{
                            std::vector<double> c = glue[i] -> getComponentX();
                            constrain[0] = 2;
                            value[0] = c[0];
                        }
                        if ((glue[i] -> getComponentY()).size() == 0){
                            constrain[1] = 2; value[1] = 0;
                        }else{
                            std::vector<double> c = glue[i] -> getComponentY();
                            constrain[1] = 2;
                            value[1] = c[0];
                        }
                        if ((glue[i] -> getComponentZ()).size() == 0){
                            constrain[2] = 2; value[2] = 0;
                        }else{
                            std::vector<double> c = glue[i] -> getComponentZ();
                            constrain[2] = 2;
                            value[2] = c[0];
                        }
                    }
                }              
                for (int i = 0; i < FSinterface.size(); i++){
                    if (name == FSinterface[i] -> getLineName()){
                        if ((FSinterface[i] -> getComponentX()).size() == 0){
                            constrain[0] = 3; value[0] = 0;
                        }else{
                            std::vector<double> c = FSinterface[i] -> getComponentX();
                            constrain[0] = 3;
                            value[0] = c[0];
                        }
                        if ((FSinterface[i] -> getComponentY()).size() == 0){
                            constrain[1] = 3; value[1] = 0;
                        }else{
                            std::vector<double> c = FSinterface[i] -> getComponentY();
                            constrain[1] = 3;
                            value[1] = c[0];
                        }
                        if ((FSinterface[i] -> getComponentZ()).size() == 0){
                            constrain[2] = 3; value[2] = 0;
                        }else{
                            std::vector<double> c = FSinterface[i] -> getComponentZ();
                            constrain[2] = 3;
                            value[2] = c[0];
                        }
                    }
                }        
                Boundaries * bound = new Boundaries(connectB, numBoundElems++, constrain, value, ibound);
                // std::cout << "asdasd " << rank << " " << ibound << std::endl;
                boundary_.push_back(bound);
            } else {
                // if(rank == 0){
                    Surface* object = geometry_ -> getSurface(name);
                    numElem++;

                    VecInt connect(nElNodes);
                    for (int j = 0 ; j < nElNodes; j++) connect[j] = elementNodes[j];

                    Element<DIM,DEG> *el = new Element<DIM,DEG>(index++,connect,this);
                    elementsAux_.push_back(el);

                    for (int k = 0; k < nElNodes; k++){
                        nodes_[connect[k]] -> pushInverseIncidence(index);
                    };
                // }
            }
        } else if ((name[0] == 'l') && (DIM == 2)) {
            VecInt connectB(nBdNodes);

            for (int i = 0; i < nBdNodes; i++) connectB[i] = elementNodes[i];
            
            int ibound;

            std::string::size_type sz;   // alias of size_t
            ibound = std::stoi (&name[1],nullptr,10);

            VecInt constrain(3);
            VecDouble value(3);

            for (int i = 0; i < dirichlet.size(); i++){
                if (name == dirichlet[i] -> getLineName()){
                    if ((dirichlet[i] -> getComponentX()).size() == 0){
                        constrain[0] = 0; value[0] = 0;
                    }else{
                        std::vector<double> c = dirichlet[i] -> getComponentX();
                        constrain[0] = 1;
                        value[0] = c[0];
                    }
                    if ((dirichlet[i] -> getComponentY()).size() == 0){
                        constrain[1] = 0; value[1] = 0;
                    }else{
                        std::vector<double> c = dirichlet[i] -> getComponentY();
                        constrain[1] = 1;
                        value[1] = c[0];
                    }
                }
            }
            for (int i = 0; i < neumann.size(); i++){
                if (name == neumann[i] -> getLineName()){
                    if ((neumann[i] -> getComponentX()).size() == 0){
                        constrain[0] = 0; value[0] = 0;
                    }else{
                        std::vector<double> c = neumann[i] -> getComponentX();
                        constrain[0] = 0;
                        value[0] = c[0];
                    }
                    if ((neumann[i] -> getComponentY()).size() == 0){
                        constrain[1] = 0; value[1] = 0;
                    }else{
                        std::vector<double> c = neumann[i] -> getComponentY();
                        constrain[1] = 0;
                        value[1] = c[0];
                    }
                }
            }  
            for (int i = 0; i < glue.size(); i++){
                if (name == glue[i] -> getLineName()){
                    if ((glue[i] -> getComponentX()).size() == 0){
                        constrain[0] = 2; value[0] = 0;
                    }else{
                        std::vector<double> c = glue[i] -> getComponentX();
                        constrain[0] = 2;
                        value[0] = c[0];
                    }
                    if ((glue[i] -> getComponentY()).size() == 0){
                        constrain[1] = 2; value[1] = 0;
                    }else{
                        std::vector<double> c = glue[i] -> getComponentY();
                        constrain[1] = 2;
                        value[1] = c[0];
                    }//std::cout <<"aqui " << std::endl;
                }
            }              
            for (int i = 0; i < FSinterface.size(); i++){
                if (name == FSinterface[i] -> getLineName()){
                    if ((FSinterface[i] -> getComponentX()).size() == 0){
                        constrain[0] = 3; value[0] = 0;
                    }else{
                        std::vector<double> c = FSinterface[i] -> getComponentX();
                        constrain[0] = 3;
                        value[0] = c[0];
                    }
                    if ((FSinterface[i] -> getComponentY()).size() == 0){
                        constrain[1] = 3; value[1] = 0;
                    }else{
                        std::vector<double> c = FSinterface[i] -> getComponentY();
                        constrain[1] = 3;
                        value[1] = c[0];
                    }
                }
            }        
            Boundaries * bound = new Boundaries(connectB, numBoundElems++, constrain, value, ibound);
            // std::cout << "asdasd " << rank << " " << ibound << std::endl;
            boundary_.push_back(bound);           
        }   
    }

    return;
}

//------------------------------------------------------------------------------
//---------------------------READS THE MESH ELEMENTS----------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Fluid<DIM,DEG>::renumberConnectivity(){
    // Renumber nodes - start
    std::vector<int > neighborNodes;

    int* xadj;
    int numNd = numNodes;
    std::vector<int> adjncy;
    xadj = new int[numNd+1]();
    adjncy.reserve(10*numNodes);

    for (int iNode = 0; iNode < numNodes; iNode++){
        neighborNodes.reserve(nodes_[iNode] -> getNumberOfElements()*3);
        neighborNodes.push_back(iNode);

        for (int j = 0; j < nodes_[iNode] -> getNumberOfElements(); j++){
            int elem = nodes_[iNode] -> getInverseIncidenceElement(j);
            VecInt connec = elements_[elem-1] -> getConnectivity();

            // std::cout << "COMM " << connec[0] << " " << connec[4] << std::endl;
            bool flag = false;
            for (int i = 0; i < nElNodes; i++){
                for (int iNeig = 0; iNeig < neighborNodes.size(); iNeig++){
                    if (connec[i] == neighborNodes[iNeig]){
                        flag = true;
                        break;
                    }
                }
                if (flag == false) neighborNodes.push_back(connec[i]);
                flag = false;
            }
        }
        //Save nodal adjacency for domain partitioning
        xadj[iNode+1] = xadj[iNode] + neighborNodes.size() - 1;
        for (int i = 1; i < neighborNodes.size(); i++){
            adjncy.push_back(neighborNodes[i]);
        }
        neighborNodes.clear();
        neighborNodes.shrink_to_fit();
    }

    //Save a second adjacency vector in idx_t format
    int *adjncy2;
    int adj_size = adjncy.size();
    adjncy2 = new int[adj_size];
    for (int i = 0; i < adjncy.size(); i++) adjncy2[i] = adjncy[i];
    
    adjncy.clear(); adjncy.shrink_to_fit();

    int* perm;
    int* iperm;
    perm = new int[numNd];
    iperm = new int[numNd];

    // Call METIS for node renumbering

    // MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
    //                         numNodes, numNodes,
    //                         1,NULL,1,NULL,&A); CHKERRQ(ierr);
    
    // // MatGetOrdering(A, MATORDERINGMETISND, IS *rperm, IS *cperm)

    // MatDestroy(&A);

    METIS_NodeND(&numNd, xadj, adjncy2, NULL, NULL, perm, iperm);

    //Reorder nodes
    for (int i = 0, j; i < numNodes; ++i) {
        for (j = iperm[i]; j < i; j = iperm[j]);
        if (j == i) while (j = iperm[j],j != i) std::swap(nodes_[i],nodes_[j]);
    }

    // Update connectivity
    for (int i = 0; i < elements_.size(); i++){
        VecInt connect = elements_[i] -> getConnectivity();

        //Reorder connectivity
        for (int k = 0; k < nElNodes; k++) connect[k] = iperm[connect[k]];
        elements_[i] -> setConnectivity(connect);
    }
    // Update boundary connectivity
    for (int ibound = 0; ibound < numBoundElems; ibound++){
        VecInt connectB = boundary_[ibound] -> getBoundaryConnectivity();

        for (int k = 0; k < nBdNodes; k++) connectB[k] = iperm[connectB[k]];
        boundary_[ibound] -> setBoundaryConnectivity(connectB);
    }
    
    for (int i = 0; i < numNodes; i++) nodes_[i] -> clearInverseIncidence();

    for (int i = 0; i < elements_.size(); i++){
        VecInt connect = elements_[i] -> getConnectivity();

        for (int k = 0; k < nElNodes; k++) nodes_[connect[k]] -> pushInverseIncidence(i);
    }

    PetscLogDouble bytes = 0;
    PetscMemoryGetCurrentUsage(&bytes);
    PetscPrintf(PETSC_COMM_WORLD,"Memory used-1 %g M\n",bytes/(1024*1024));
    
    for (int i = 0; i < numNodes; i++){
        for (int j = 0; j < nodes_[i] -> getNumberOfElements(); j++){
            int elJ = nodes_[i] -> getInverseIncidenceElement(j);
            for (int k = 0; k < nodes_[i] -> getNumberOfElements(); k++)
                elements_[elJ] -> pushNeighborElement(nodes_[i] -> getInverseIncidenceElement(k)); 
        }
    }

    PetscMemoryGetCurrentUsage(&bytes);
    PetscPrintf(PETSC_COMM_WORLD,"Memory used00 %g M\n",bytes/(1024*1024));

    for (int i = 0; i < elements_.size(); i++) elements_[i] -> sortEraseNeighborElements();

    PetscMemoryGetCurrentUsage(&bytes);
    PetscPrintf(PETSC_COMM_WORLD,"Memory used11 %g M\n",bytes/(1024*1024));        

    // for (int i = 0; i < elements_.size(); i++){
    //     if (rank == 0) std::cout << "Neighbor " << i << " ";
    //     for (int j = 0; j < elements_[i] -> getNumberOfNeighborElements(); j++){
    //         std::cout << elements_[i] -> getNeighborElement(j) << " ";
    //     }
    //     std::cout << std::endl;
    // }

    delete [] perm;
    delete [] iperm;
    delete [] adjncy2;
    delete [] xadj;
    // Renumber nodes - end

    


    return;
}

//------------------------------------------------------------------------------
//-------------------------SETS THE BOUNDARY CONDITIONS-------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Fluid<DIM,DEG>::setBoundaryConstrains(){

    if (rank == 0) std::cout << "8/9 Setting boundary conditions..." << std::endl;

    for (int ibound = 0; ibound < numBoundElems; ibound++){
        
        VecInt connectB = boundary_[ibound] -> getBoundaryConnectivity();

        for (int k = 0; k < DIM; k++){
            if ((boundary_[ibound] -> getConstrain(k) == 3)){
                for (int j = 0; j < nBdNodes; j++) nodes_[connectB[j]] -> setConstrainsLaplace(k,1,0);
            };
        };

        for (int k = 0; k < DIM; k++){
            if ((boundary_[ibound] -> getConstrain(k) == 1) || (boundary_[ibound] -> getConstrain(k) == 3)){
                for (int j = 0; j < nBdNodes; j++) 
                    nodes_[connectB[j]] -> setConstrains(k,boundary_[ibound] -> getConstrain(k),
                                                         boundary_[ibound] -> getConstrainValue(k));
            };
        }
    };


    return;
}

//------------------------------------------------------------------------------
//DEFINES IF THE ELEMENT BELONGS TO THE BOUNDARY AND THE RESPECTIVE ELEMENT SIDE
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Fluid<DIM,DEG>::setBoundarySides(){

    if (rank == 0) std::cout << "9/9 Setting boundary sides..." << std::endl;

    //Sets fluid elements and sides on interface boundaries
    for (int i=0; i<numBoundElems; i++){

        int group = boundary_[i] -> getBoundaryGroup();

       if ((boundary_[i] -> getConstrain(0) > 0) || (boundary_[i] -> getConstrain(1) > 0)) {
            
            VecInt connectB = boundary_[i] -> getBoundaryConnectivity();

            for (int j=0; j<numElem; j++){
                VecInt connect = elements_[j] -> getConnectivity();

                int flag = 0;
            
                int side[nBdNodes];
                for (int k=0; k<nElNodes; k++){
                    for (int l = 0; l<nBdNodes; l++){
                        if (connectB[l] == connect[k]){
                            side[flag] = k;
                            flag++;
                        }
                    };
                };
                if (flag == nBdNodes){
                    boundary_[i] -> setElement(j);
                    //Sets element index and side
                    // for (int k=0; k<nBdNodes; k++) std::cout << "BD NODES " << i << " " << j << " " << k << " " << side[k] << std::endl;
                    
                    for (int k=0; k<DIM+1; k++){
                        // std::cout << "DDDDD " << k << std::endl;

                        int* end = side + nBdNodes;
                        int* foo = std::find(side, end, k);

                        if ((foo == end) && (elements_[j] -> getElemSideInBoundary() < 0)){
                            int aux = boundary_[i] -> getBoundaryGroup();
                            boundary_[i] -> setBoundaryGroup(aux);
                            boundary_[i] -> setElementSide(k);
                            boundary_[i] -> setElement(j);
                            elements_[j] -> setElemSideInBoundary(k);
                            if (boundary_[i] -> getConstrain(0) == 3) elements_[j] -> setFSIInterface();
                        }
                    }
                };
            };
        }
    };
    
    return;
}

//------------------------------------------------------------------------------
//---------------------SUBDIVIDES THE FINITE ELEMENT DOMAIN---------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Fluid<DIM,DEG>::domainDecompositionMETIS() {
    
    std::string mirror2;
    mirror2 = "domain_decomposition.txt";
    std::ofstream mirrorData(mirror2.c_str());
    
    int size;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    idx_t objval;
    idx_t numEl = numElem;
    idx_t numNd = numNodes;
    idx_t ssize = size;
    idx_t one = 1;
    int elem_start[numElem+1], elem_connec[nElNodes*numElem];
    part_elem = new int[numElem];
    part_nodes = new int[numNodes];


    for (int i = 0; i < numElem+1; i++){
        elem_start[i]=nElNodes*i;
    };
    for (int jel = 0; jel < numElem; jel++){
        auto connec=elements_[jel]->getConnectivity();        
        
        for (int i=0; i<nElNodes; i++){
        elem_connec[nElNodes*jel+i] = connec[i];
        };
    };

    //Performs the domain decomposition
    if (size == 1){
        for (int i = 0; i < numNodes; i++) part_nodes[i] = 0;
        for (int i = 0; i < numElem; i++) part_elem[i] = 0;
    } else {
        METIS_PartMeshDual(&numEl, &numNd, elem_start, elem_connec, \
                                NULL, NULL, &one, &ssize, NULL, NULL,    \
                                &objval, part_elem, part_nodes);
    }
    
    mirrorData << std::endl \
               << "FLUID MESH DOMAIN DECOMPOSITION - ELEMENTS" << std::endl;
    for(int i = 0; i < numElem; i++){
        mirrorData << "process = " << part_elem[i] \
                   << ", element = " << i << std::endl;
    };

    mirrorData << std::endl \
               << "FLUID MESH DOMAIN DECOMPOSITION - NODES" << std::endl;
    for(int i = 0; i < numNodes; i++){
        mirrorData << "process = " << part_nodes[i] \
                   << ", node = " << i << std::endl;
    };


};

//------------------------------------------------------------------------------
//----------------------COMPUTES DRAG AND LIFT COEFFICIENTS---------------------
//------------------------------------------------------------------------------
template<>
void Fluid<2,2>::dragAndLiftCoefficients(std::ofstream& dragLift){

    double dragCoefficient = 0.;
    double liftCoefficient = 0.;
    double pressureDragCoefficient = 0.;
    double pressureLiftCoefficient = 0.;
    double frictionDragCoefficient = 0.;
    double frictionLiftCoefficient = 0.;
    
    for (int jel = 0; jel < numBoundElems; jel++){   
        
        double rhoInf = 1.0;
        double velocityInf[2];
        velocityInf[0] = -1.;
        velocityInf[1] = 0.;
        
        double dForce = 0.;
        double lForce = 0.;
        double pDForce = 0.;
        double pLForce = 0.;
        double fDForce = 0.;
        double fLForce = 0.;
        double aux_Mom = 0.;
        double aux_Per = 0.;
        
       for (int i=0; i<numberOfLines; i++){
            //std::cout << "Bound group " << boundary_[jel] -> getBoundaryGroup() << std::endl;
            if (boundary_[jel] -> getBoundaryGroup() == dragAndLiftBoundary[i]){
                //std::cout << "AQUI " << numberOfLines<< " " << i << " " << dragAndLiftBoundary[i] << std::endl;
                int iel = boundary_[jel] -> getElement();
                // elements_[iel] -> computeDragAndLiftForces(pDForce, pLForce, fDForce, fLForce, dForce, lForce, aux_Mom, aux_Per);
                // elements_[iel] -> computeSeparationAngle();
            };
        };
        
        pressureDragCoefficient += pDForce ;/// 
            //(0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        pressureLiftCoefficient += pLForce;// / 
            //(0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        
        frictionDragCoefficient += fDForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        frictionLiftCoefficient += fLForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        
        dragCoefficient += dForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        liftCoefficient += lForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        
    };
    // std::cout << "vazao " << pressureDragCoefficient << " " << pressureLiftCoefficient << std::endl;
    if (rank == 0) {
        const int timeWidth = 11;
        const int numWidth = 11;
        dragLift << std::setprecision(3) << std::scientific;
        dragLift << std::left << std::setw(timeWidth) << iTimeStep * dTime;
        dragLift << std::setw(numWidth) << pressureDragCoefficient;
        dragLift << std::setw(numWidth) << pressureLiftCoefficient;
        dragLift << std::setw(numWidth) << frictionDragCoefficient;
        dragLift << std::setw(numWidth) << frictionLiftCoefficient;
        dragLift << std::setw(numWidth) << dragCoefficient;
        dragLift << std::setw(numWidth) << liftCoefficient;
        dragLift << std::endl;
    }
}


template<int DIM, int DEG>
void Fluid<DIM,DEG>::meshReading(Geometry* &geometry_, const std::string& inputFile, const std::string& inputMesh, const std::string& mirror, const bool& deleteFiles) {

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++OPPENING FILES+++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //opening the .msh file
    
    std::ofstream mirrorData(mirror.c_str());
    std::ifstream file(inputMesh);
    std::string line;
    std::getline(file, line); std::getline(file, line); std::getline(file, line); std::getline(file, line);
  

    readInputFile(inputFile,mirrorData);
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++READIN MESH+++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++PHYSICAL ENTITIES++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int number_physical_entities;
    file >> number_physical_entities;
    std::getline(file, line);
    std::unordered_map<int, std::string> physicalEntities;
    physicalEntities.reserve(number_physical_entities);

    for (int i = 0; i < number_physical_entities; i++)
    {
        std::getline(file, line);
        std::vector<std::string> tokens = split2(line, " ");
        int index;
        std::istringstream(tokens[1]) >> index;
        physicalEntities[index] = tokens[2].substr(1, tokens[2].size() - 2);
    }
    std::getline(file, line); std::getline(file, line);

    numIntegration = new DIntegration();

    readNodes(file,mirrorData);
    readElements(geometry_,file,mirrorData,elements_, physicalEntities);

    if (fProbType == ProblemType::ENavierStokes || fProbType == ProblemType::EStokes){
        numDOF = (DIM+1) * numNodes;
    } else if (fProbType == ProblemType::EPoisson) {
        numDOF = numNodes;
    } else if (fProbType == ProblemType::EElastic){
        numDOF = numNodes * DIM;
    } else {
        PanicButton();
    }
    
    renumberConnectivity();

    // if (rank == 0) std::cout << "Number of elements " << number_elements << " " 
                             // << numElem << " " << numBoundElems << std::endl;
    mirrorData << std::endl << "Element Connectivity" << std::endl;        
    
    for (int jel = 0; jel < numElem; jel++){
        VecInt connec = elements_[jel] -> getConnectivity();       
        for (int i=0; i < nElNodes; i++){
            mirrorData << connec[i] << " ";
        };
        mirrorData << std::endl;
    };

    setBoundaryConstrains();

    //Print nodal constrains
    for (int i=0; i<numNodes; i++){

        mirrorData<< "Constrains " << i
                  << " " << nodes_[i] -> getConstrains(0)
                  << " " << nodes_[i] -> getConstrainValue(0)
                  << " " << nodes_[i] -> getConstrains(1)
                  << " " << nodes_[i] -> getConstrainValue(1) << std::endl;
    }; 

    for (int i=0; i<numBoundElems; i++){

        mirrorData<< "Bound Elements " << i
                  << " " << boundary_[i] -> getBoundaryGroup() << std::endl;
    }; 

    setBoundarySides();
 
    domainDecompositionMETIS();

    iAux = 0;


    //Closing the file
    file.close();
    if (deleteFiles)
        system((rm + inputFile).c_str());

    // printResults(100);

return;
};

//------------------------------------------------------------------------------
//-------------------------SOLVE STEADY LAPLACE PROBLEM-------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
int Fluid<DIM,DEG>::solveSteadyLaplaceProblem(int iterNumber, double tolerance) {

    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione, iterations;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
   
    int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);        

    for (int inewton = 0; inewton < iterNumber; inewton++){

        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                            2*numNodes, 2*numNodes,
                            100,NULL,100,NULL,&A); CHKERRQ(ierr);
        
        ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
        
        //Create PETSc vectors
        ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
        ierr = VecSetSizes(b,PETSC_DECIDE,2*numNodes);CHKERRQ(ierr);
        ierr = VecSetFromOptions(b);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&All);CHKERRQ(ierr);
        
        //std::cout << "Istart = " << Istart << " Iend = " << Iend << std::endl;
        
        for (int jel = 0; jel < numElem; jel++){               
            //Compute Element matrix
            VecInt connec = elements_[jel] -> getConnectivity();

            MatrixDouble matrix(nLocDOF,nLocDOF);
            VecDouble rhs(nLocDOF);
            rhs.setZero();
            matrix.setZero();

            elements_[jel] -> getSolidProblem(matrix,rhs);
            
            //Disperse local contributions into the global matrix
            //Matrix K and C
            for (int i=0; i<nElNodes; i++){
                for (int j=0; j<nElNodes; j++){
                    if (fabs(matrix(2*i  ,2*j  )) >= 1.e-8){
                        int dof_i = 2*connec[i];
                        int dof_j = 2*connec[j];
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(2*i  ,2*j  ),ADD_VALUES);
                    };
                    if (fabs(matrix(2*i+1,2*j  )) >= 1.e-8){
                        int dof_i = 2*connec[i]+1;
                        int dof_j = 2*connec[j];
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(2*i+1,2*j  ),ADD_VALUES);
                    };
                    if (fabs(matrix(2*i  ,2*j+1)) >= 1.e-8){
                        int dof_i = 2*connec[i];
                        int dof_j = 2*connec[j]+1;
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(2*i  ,2*j+1),ADD_VALUES);
                    };
                    if (fabs(matrix(2*i+1,2*j+1)) >= 1.e-8){
                        int dof_i = 2*connec[i]+1;
                        int dof_j = 2*connec[j]+1;
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(2*i+1,2*j+1),ADD_VALUES);
                    };
                };
                                    
                //Rhs vector
                if (fabs(rhs[2*i  ]) >= 1.e-8){
                    int dof_i = 2*connec[i];
                    ierr = VecSetValues(b,1,&dof_i,&rhs[2*i  ],ADD_VALUES);
                };
                
                if (fabs(rhs[2*i+1]) >= 1.e-8){
                    int dof_i = 2*connec[i]+1;
                    ierr = VecSetValues(b,1,&dof_i,&rhs[2*i+1],ADD_VALUES);
                };
            };
        };
        
        //Assemble matrices and vectors
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        
        ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
        
        //MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        //ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        
        //Create KSP context to solve the linear system
        ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
        
        ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
        
#if defined(PETSC_HAVE_MUMPS)
        ierr = KSPSetType(ksp,KSPPREONLY);
        ierr = KSPGetPC(ksp,&pc);
        ierr = PCSetType(pc, PCLU);
#endif
        
        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        ierr = KSPSetUp(ksp);
        
        
        
        ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
        
        ierr = KSPGetTotalIterations(ksp, &iterations);

        //std::cout << "GMRES Iterations = " << iterations << std::endl;
        
        //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKERRQ(ierr);
        
        //Gathers the solution vector to the master process
        ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
        
        ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        
        ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        
        ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
                
        //Updates nodal values
        double u_ [2];
        Ione = 1;

        for (int i = 0; i < numNodes; ++i){
            Ii = 2*i;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            u_[0] = val;
            Ii = 2*i+1;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            u_[1] = val;
            if (nodes_[i] -> getConstrainsLaplace(0) != 1) nodes_[i] -> incrementCoordinate(0,u_[0]);
            if (nodes_[i] -> getConstrainsLaplace(1) != 1) nodes_[i] -> incrementCoordinate(1,u_[1]);




            // double x = nodes_[i] -> getCoordinateValue(0);
            // double y = nodes_[i] -> getCoordinateValue(1);
            // double xp = nodes_[i] -> getPreviousCoordinateValue(0);
            // double yp = nodes_[i] -> getPreviousCoordinateValue(1);
            // double vx = nodes_[i] -> getMeshVelocity(0);
            // double vy = nodes_[i] -> getMeshVelocity(1);
            // double ax = nodes_[i] -> getMeshAcceleration(0);
            // double ay = nodes_[i] -> getMeshAcceleration(1);

            // double accelx = (x - xp) / (0.25 * dTime * dTime) - vx / (0.25 * dTime) - ax * (0.5/0.25 - 1.0);
            // double accely = (y - yp) / (0.25 * dTime * dTime) - vy / (0.25 * dTime) - ay * (0.5/0.25 - 1.0);

            // nodes_[i] -> setMeshAccelerationComponent(0,accelx);
            // nodes_[i] -> setMeshAccelerationComponent(1,accely);

            // double velx = 0.5 * dTime * accelx + vx + dTime * (1.0 - 0.5) * ax;
            // double vely = 0.5 * dTime * accely + vy + dTime * (1.0 - 0.5) * ay;
            // // std::cout << "asd as " << vx << " " << velx << " " << accelx << " " << ax << std::endl;
            // nodes_[i] -> setMeshVelocityComponent(0,velx);
            // nodes_[i] -> setMeshVelocityComponent(1,vely);

        };
        
        //Computes the solution vector norm
        ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
        
        if(rank == 0){
            std::cout << "MESH MOVING - ERROR = " << val 
                      << std::scientific <<  std::endl;
        };

        ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
        ierr = VecDestroy(&b); CHKERRQ(ierr);
        ierr = VecDestroy(&u); CHKERRQ(ierr);
        ierr = VecDestroy(&All); CHKERRQ(ierr);
        ierr = MatDestroy(&A); CHKERRQ(ierr);

        if(val <= tolerance){
            break;
        };
    };
    
    // for (int i=0; i<numElem; i++){
    //     elements_[i] -> computeNodalGradient();            
    // };

    if (rank == 0) {
        //Computing velocity divergent
        //      printResults(1);
    };

    return 0;
};

//------------------------------------------------------------------------------
//-------------------------SOLVE TRANSIENT FLUID PROBLEM------------------------
//------------------------------------------------------------------------------
template<>
void Fluid<2,2>::readInitialValues(const std::string& inputPrev, const std::string& inputCurr) {

    int rank;

    // MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // hid_t filePrevious;
    // hid_t fileCurrent;
    // hid_t dataset;
    // herr_t status;

    // double *vecValues;
    // double *scaValues;
    // vecValues = new double[3*numNodes];
    // scaValues = new double[numNodes];

    // filePrevious = H5Fopen(inputPrev.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
    // fileCurrent = H5Fopen(inputCurr.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);

    // char datasetName[] = "/velocity";
    // dataset = H5Dopen( filePrevious, datasetName, H5P_DEFAULT );
    // status = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecValues[0]);
    // status = H5Dclose(dataset);
    // for (int i = 0; i < numNodes; ++i){
    //     nodes_[i] -> setPreviousVelocityComponent(0,vecValues[3*i  ]);
    //     nodes_[i] -> setPreviousVelocityComponent(1,vecValues[3*i+1]);
    // }
    // dataset = H5Dopen( fileCurrent, datasetName, H5P_DEFAULT );
    // status = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecValues[0]);
    // status = H5Dclose(dataset);
    // for (int i = 0; i < numNodes; ++i){
    //     nodes_[i] -> setVelocityComponent(0,vecValues[3*i  ]);
    //     nodes_[i] -> setVelocityComponent(1,vecValues[3*i+1]);
    // }

    // char datasetName2[] = "/acceleration";
    // dataset = H5Dopen( filePrevious, datasetName2, H5P_DEFAULT );
    // status = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecValues[0]);
    // status = H5Dclose(dataset);
    // for (int i = 0; i < numNodes; ++i){
    //     nodes_[i] -> setPreviousAccelerationComponent(0,vecValues[3*i  ]);
    //     nodes_[i] -> setPreviousAccelerationComponent(1,vecValues[3*i+1]);
    // }
    // dataset = H5Dopen( fileCurrent, datasetName2, H5P_DEFAULT );
    // status = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecValues[0]);
    // status = H5Dclose(dataset);
    // for (int i = 0; i < numNodes; ++i){
    //     nodes_[i] -> setAccelerationComponent(0,vecValues[3*i  ]);
    //     nodes_[i] -> setAccelerationComponent(1,vecValues[3*i+1]);
    // }

    // char datasetName3[] = "/lagrangeMultiplers";
    // dataset = H5Dopen( fileCurrent, datasetName3, H5P_DEFAULT );
    // status = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecValues[0]);
    // status = H5Dclose(dataset);
    // for (int i = 0; i < numNodes; ++i){
    //     nodes_[i] -> setLagrangeMultiplier(0,vecValues[3*i  ]);
    //     nodes_[i] -> setLagrangeMultiplier(1,vecValues[3*i+1]);
    // }

    // char datasetName4[] = "/pressure";
    // dataset = H5Dopen( fileCurrent, datasetName4, H5P_DEFAULT );
    // status = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &scaValues[0]);
    // status = H5Dclose(dataset);
    // for (int i = 0; i < numNodes; ++i){
    //     nodes_[i] -> setPressure(vecValues[i]);
    // }

    // delete [] vecValues;
    // delete [] scaValues;

    // //End HDF5 file
    // status = H5Fclose(filePrevious);

    return;
}


//------------------------------------------------------------------------------
//-------------------------SOLVE TRANSIENT FLUID PROBLEM------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
int Fluid<DIM,DEG>::solveFSIFluid(int iterNumber, double tolerance, int problem_type){

    Vec               b, u, All, Allu;
    PetscErrorCode    ierr;
    PetscInt          Ii, Ione, iterations;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    //    MatNullSpace      nullsp;

    int rank;

    iAux++;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    //Check if the problem type can be computed
    if ((problem_type > 2) || (problem_type < 1)){
        std::cout << "WRONG PROBLEM TYPE." << std::endl;
        return 0;
    };
            
    double duNorm=100.;

    double &alpha_f = fluidParameters.getAlphaF();
    double &alpha_m = fluidParameters.getAlphaM();
    double &gamma = fluidParameters.getGamma();
         
    for (int inewton = 0; inewton < iterNumber; inewton++){
        
        std::clock_t t1 = std::clock();
        
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                            2*numNodes+numNodes, 2*numNodes+numNodes,
                            100,NULL,300,NULL,&A); 
        CHKERRQ(ierr);
                
        //Create PETSc vectors
        ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
        ierr = VecSetSizes(b,PETSC_DECIDE,2*numNodes+numNodes);
        CHKERRQ(ierr);
        ierr = VecSetFromOptions(b);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&All);CHKERRQ(ierr);
        
        //std::cout << "Istart = " << Istart << " Iend = " << Iend << std::endl;

        for (int jel = 0; jel < numElem; jel++){   
            
            if (part_elem[jel] == rank) {
                //Compute Element matrix
                VecInt connec = elements_[jel] -> getConnectivity();

                MatrixDouble matrix(18,18);
                matrix.setZero();
                VecDouble rhs(18);
                rhs.setZero();

                elements_[jel] -> getTransientNavierStokes(matrix,rhs);
                
                //Disperse local contributions into the global matrix
                //Matrix K and C
                for (int i=0; i<6; i++){
                    for (int j=0; j<6; j++){
                        int dof_i = 2 * connec[i];
                        int dof_j = 2 * connec[j];
                        MatSetValues(A, 1, &dof_i,1, &dof_j, &matrix(2*i  ,2*j  ), ADD_VALUES);
                        
                        dof_i = 2 * connec[i] + 1;
                        dof_j = 2 * connec[j];
                        MatSetValues(A, 1, &dof_i, 1, &dof_j, &matrix(2*i+1,2*j  ), ADD_VALUES);
                        
                        dof_i = 2 * connec[i];
                        dof_j = 2 * connec[j] + 1;
                        MatSetValues(A, 1, &dof_i, 1, &dof_j, &matrix(2*i  ,2*j+1), ADD_VALUES);

                        dof_i = 2 * connec[i] + 1;
                        dof_j = 2 * connec[j] + 1;
                        MatSetValues(A, 1, &dof_i, 1, &dof_j, &matrix(2*i+1,2*j+1), ADD_VALUES);
                        
                        //Matrix Q and Qt
                        dof_i = 2 * connec[i];
                        dof_j = 2 * numNodes + connec[j];
                        MatSetValues(A, 1, &dof_i, 1, &dof_j, &matrix(2*i  ,12+j), ADD_VALUES);
                    
                        dof_i = 2 * connec[i];
                        dof_j = 2 * numNodes + connec[j];
                        MatSetValues(A, 1, &dof_j, 1, &dof_i, &matrix(12+j,2*i  ), ADD_VALUES);
                        
                        dof_i = 2 * connec[i] + 1;
                        dof_j = 2 * numNodes + connec[j];
                        MatSetValues(A, 1, &dof_i, 1, &dof_j, &matrix(2*i+1,12+j), ADD_VALUES);
                        
                        dof_i = 2 * connec[i] + 1;
                        dof_j = 2 * numNodes + connec[j];
                        MatSetValues(A, 1, &dof_j, 1, &dof_i, &matrix(12+j,2*i+1), ADD_VALUES);
                        
                        dof_i = 2 * numNodes + connec[i];
                        dof_j = 2 * numNodes + connec[j];
                        MatSetValues(A, 1, &dof_i, 1, &dof_j, &matrix(12+i,12+j), ADD_VALUES);
                    };
                    
                    //Rhs vector
                    int dof_i = 2 * connec[i];
                    VecSetValues(b, 1, &dof_i, &rhs[2*i  ], ADD_VALUES);
                    
                    dof_i = 2 * connec[i]+1;
                    VecSetValues(b, 1, &dof_i, &rhs[2*i+1], ADD_VALUES);
                    
                    dof_i = 2 * numNodes + connec[i];
                    VecSetValues(b, 1, &dof_i, &rhs[12+i], ADD_VALUES);
                };
            };
        }; //Elements
        
        //Assemble matrices and vectors
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        
        ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
        
        // MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        // ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        
        //Create KSP context to solve the linear system
        ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
        
        ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
        



    //     ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,
    //                             500);CHKERRQ(ierr);
        
    //     ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        
    //     ierr = KSPGetPC(ksp,&pc);
        
    //     ierr = PCSetType(pc,PCJACOBI);
        
    //     //ierr = KSPSetType(ksp,KSPBCGS); CHKERRQ(ierr);

    //     // ierr = KSPGMRESSetRestart(ksp, 10); CHKERRQ(ierr);
        
    //        //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
        

    //   //   ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,NULL,&nullsp);
    // // ierr = MatSetNullSpace(A, nullsp);
    // // ierr = MatNullSpaceDestroy(&nullsp);



#if defined(PETSC_HAVE_MUMPS)
        ierr = KSPSetType(ksp,KSPPREONLY);
        ierr = KSPGetPC(ksp,&pc);
        ierr = PCSetType(pc, PCLU);
#endif          
        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        ierr = KSPSetUp(ksp);



        ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);

        ierr = KSPGetTotalIterations(ksp, &iterations);            

        //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKERRQ(ierr);
        
        //Gathers the solution vector to the master process
        ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
        ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);

        ierr = VecScatterCreateToAll(b, &ctx, &Allu);CHKERRQ(ierr);
        ierr = VecScatterBegin(ctx, b, Allu, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        ierr = VecScatterEnd(ctx, b, Allu, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
        
        //Updates nodal values
        double p_;
        duNorm = 0.;
        double dpNorm = 0.;
        Ione = 1;
        
        for (int i = 0; i < numNodes; ++i){
            Ii = 2*i;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            nodes_[i] -> incrementAcceleration(0,val);
            nodes_[i] -> incrementVelocity(0,val*gamma*dTime);

            ierr = VecGetValues(Allu, Ione, &Ii, &val);CHKERRQ(ierr);
            duNorm += val*val;
        
            Ii = 2*i+1;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            nodes_[i] -> incrementAcceleration(1,val);
            nodes_[i] -> incrementVelocity(1,val*gamma*dTime);

            ierr = VecGetValues(Allu, Ione, &Ii, &val);CHKERRQ(ierr);
            duNorm += val*val;
        };
        
        for (int i = 0; i<numNodes; i++){
            Ii = 2*numNodes+i;
            ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
            p_ = val;
            nodes_[i] -> incrementPressure(p_);

            ierr = VecGetValues(Allu,Ione,&Ii,&val);CHKERRQ(ierr);
            dpNorm += val*val;
        };
        
        //Computes the solution vector norm
        //ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);

        std::clock_t t2 = std::clock();
     
        if(rank == 0){

            std::cout << "Iteration = " << inewton 
                      << " (" << iterations << ")"  
                      << "   Du Norm = " << std::scientific << sqrt(duNorm) 
                      << " " << sqrt(dpNorm)
                      << "  Time (s) = " << std::fixed
                      << 1000.*(t2-t1)/CLOCKS_PER_SEC/1000. << std::endl;
        };
                  
        ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
        ierr = VecDestroy(&b); CHKERRQ(ierr);
        ierr = VecDestroy(&u); CHKERRQ(ierr);
        ierr = VecDestroy(&All); CHKERRQ(ierr);
        ierr = VecDestroy(&Allu); CHKERRQ(ierr);
        ierr = MatDestroy(&A); CHKERRQ(ierr);

        if (sqrt(duNorm) <= tolerance) {
            break;
        };



    };//Newton-Raphson
    
    return 0;
};



//------------------------------------------------------------------------------
//-------------------------SOLVE TRANSIENT FLUID PROBLEM------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
int Fluid<DIM,DEG>::solvePoisson(){

    Vec               b, u, All, Allu;
    PetscErrorCode    ierr;
    PetscInt          Ii, Ione, iterations;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    //    MatNullSpace      nullsp;

    int rank;

    iAux++;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
            
    double duNorm=100.;
                 
        std::clock_t t1 = std::clock();
        
        if (fluidParameters.getSolverType() == SolverType::ESuiteSparse){
            ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, numNodes, numNodes, 100,NULL,&A);
        } else {
            ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                            numNodes, numNodes,100,NULL,300,NULL,&A); 
        }

        CHKERRQ(ierr);
                
        //Create PETSc vectors
        ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
        ierr = VecSetSizes(b,PETSC_DECIDE,numNodes);
        CHKERRQ(ierr);
        ierr = VecSetFromOptions(b);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&All);CHKERRQ(ierr);
        
        //std::cout << "Istart = " << Istart << " Iend = " << Iend << std::endl;

        for (int jel = 0; jel < numElem; jel++){   
            
            if (part_elem[jel] == rank) {
                //Compute Element matrix
                VecInt connec = elements_[jel] -> getConnectivity();

                MatrixDouble matrix(nElNodes,nElNodes);
                matrix.setZero();
                VecDouble rhs(nElNodes);
                rhs.setZero();

                elements_[jel] -> getPoisson(matrix,rhs);
                
                //Disperse local contributions into the global matrix
                //Matrix K and C
                for (int i=0; i<nElNodes; i++){
                    for (int j=0; j<nElNodes; j++){
                        int dof_i = connec[i];
                        int dof_j = connec[j];
                        MatSetValues(A, 1, &dof_i,1, &dof_j, &matrix(i,j), ADD_VALUES);
                    };
                    
                    //Rhs vector
                    int dof_i = connec[i];
                    VecSetValues(b, 1, &dof_i, &rhs[i], ADD_VALUES);
                };
            };
        }; //Elements
        
        //Assemble matrices and vectors
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        
        ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
        
        // MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        // ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        
        //Create KSP context to solve the linear system
        ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
        
        ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
        
        switch (fluidParameters.getSolverType())
        {
        case SolverType::ESuiteSparse:
            KSPGetPC(ksp, &pc);
            PCSetType(pc, PCLU);
            PCFactorSetMatSolverType(pc, MATSOLVERUMFPACK);
            break;
        case SolverType::EMumps:
            KSPGetPC(ksp, &pc);
            PCSetType(pc, PCLU);
            PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
            break;

        case SolverType::EIterative:
            KSPSetType(ksp,KSPFGMRES);
            KSPGetPC(ksp, &pc);
            PCSetType(pc,PCBJACOBI);
            KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,200);
            break;

        default:
            PanicButton();
            break;
        }

    //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
        
        ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);

        ierr = KSPGetTotalIterations(ksp, &iterations);            
// 
        // ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKERRQ(ierr);
        
        //Gathers the solution vector to the master process
        ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
        ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);

        ierr = VecScatterCreateToAll(b, &ctx, &Allu);CHKERRQ(ierr);
        ierr = VecScatterBegin(ctx, b, Allu, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        ierr = VecScatterEnd(ctx, b, Allu, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
        
        //Updates nodal values
        double p_;
        duNorm = 0.;
        double dpNorm = 0.;
        Ione = 1;
        
        for (int i = 0; i < numNodes; ++i){
            Ii = i;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            nodes_[i] -> incrementVelocity(0,val);
        };
        
        //Computes the solution vector norm
        //ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);

        std::clock_t t2 = std::clock();
     
        if(rank == 0){

            std::cout << "  Time (s) = " << std::fixed
                      << 1000.*(t2-t1)/CLOCKS_PER_SEC/1000. << std::endl;
        };
                  
        ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
        ierr = VecDestroy(&b); CHKERRQ(ierr);
        ierr = VecDestroy(&u); CHKERRQ(ierr);
        ierr = VecDestroy(&All); CHKERRQ(ierr);
        ierr = VecDestroy(&Allu); CHKERRQ(ierr);
        ierr = MatDestroy(&A); CHKERRQ(ierr);

        VecDouble errorsTotal(3);
    if (fluidParameters.getExactSolutionPoisson()) computeError(errorsTotal);

    printResultsPoisson();

    
    return 0;
};


//------------------------------------------------------------------------------
//-------------------------SOLVE TRANSIENT FLUID PROBLEM------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Fluid<DIM,DEG>::printResultsPoisson(){

     //    std::cout << "Printing Velocity Results" << std::endl;
    std::string s = "ResultsPoisson.vtu";
    
    std::fstream output_v(s.c_str(), std::ios_base::out);

    output_v << "<?xml version=\"1.0\"?>" << std::endl
             << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
             << "  <UnstructuredGrid>" << std::endl
             << "  <Piece NumberOfPoints=\"" << numNodes
             << "\"  NumberOfCells=\"" << numElem
             << "\">" << std::endl;

    //WRITE NODAL COORDINATES
    output_v << "    <Points>" << std::endl
             << "      <DataArray type=\"Float64\" "
             << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

    for (int i=0; i<numNodes; i++){
        auto x=nodes_[i]->getCoordinates();
        output_v << x[0] << " " << x[1] << " " << 0.0 << std::endl;        
    };
    output_v << "      </DataArray>" << std::endl
             << "    </Points>" << std::endl;
    
    //WRITE ELEMENT CONNECTIVITY
    output_v << "    <Cells>" << std::endl
             << "      <DataArray type=\"Int32\" "
             << "Name=\"connectivity\" format=\"ascii\">" << std::endl;
    
    for (int i=0; i<numElem; i++){
        auto connec=elements_[i]->getConnectivity();
        output_v << connec[0] << " " << connec[1] << " " << connec[2] ;
        if (DEG == 2){
            output_v << " " << connec[3] << " " << connec[4] << " " << connec[5] ;
        } else if (DEG == 3){
            output_v << " " << connec[3] << " " << connec[4] << " " << connec[5] << " " << connec[6] << " " << connec[7] << " " << connec[8] << " " << connec[9] ;
        }
        output_v << std::endl;
    };
    output_v << "      </DataArray>" << std::endl;
  
    //WRITE OFFSETS IN DATA ARRAY
    output_v << "      <DataArray type=\"Int32\""
             << " Name=\"offsets\" format=\"ascii\">" << std::endl;
    
    int aux = 0;
    for (int i=0; i<numElem; i++){
        output_v << aux + nElNodes << std::endl;
        aux += nElNodes;
    };
    output_v << "      </DataArray>" << std::endl;
  
    //WRITE ELEMENT TYPES
    output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
             << "format=\"ascii\">" << std::endl;
    
    int val = 0;
    if (DEG == 1)val = 5;
    if (DEG == 2)val = 22;
    if (DEG == 3)val = 69;

    for (int i=0; i<numElem; i++){
        output_v << val << std::endl;
    };

    output_v << "      </DataArray>" << std::endl
             << "    </Cells>" << std::endl;

    //WRITE NODAL RESULTS
    output_v << "    <PointData>" << std::endl;

    if (printVelocity){
        output_v<< "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numNodes; i++){
            output_v << nodes_[i] -> getVelocity(0) << " "             
                     << nodes_[i] -> getVelocity(1) << " " << 0. << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    };

    
   

    output_v << "    </PointData>" << std::endl; 

    //WRITE ELEMENT RESULTS
    output_v << "    <CellData>" << std::endl;
    
    if (printProcess){
        output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                 << "Name=\"Process\" format=\"ascii\">" << std::endl;
        for (int i=0; i<numElem; i++){
            output_v << part_elem[i] << std::endl;
        };
        output_v << "      </DataArray> " << std::endl;
    };
    
    output_v << "    </CellData>" << std::endl; 

    //FINALIZE OUTPUT FILE
    output_v << "  </Piece>" << std::endl;
    
    // output_v << "  <FieldData>" << std::endl;

    // output_v << "      <DataArray type=\"Float64\" Name=\"Time\" NumberOfTuples=\"1\" "
    //          << " format=\"ascii\">" << std::endl;
    // output_v << step << std::endl;
    // output_v << "      </DataArray> " << std::endl;

    // output_v << "      <DataArray type=\"Float64\" Name=\"LiftCoefficient\" NumberOfTuples=\"1\" "
    //          << " format=\"ascii\">" << std::endl;
    // output_v << liftCoefficient << std::endl;
    // output_v << "      </DataArray> " << std::endl;

    // // output_v << "      <DataSet type=\"Float64\" Name=\"Drag Coefficient\" NumberOfTuples=\"1\" "
    // //          << " format=\"ascii\">" << std::endl;
    // // output_v << dragCoefficient << std::endl;
    // // output_v << "      </DataSet> " << std::endl;

    // output_v << "  </FieldData>" << std::endl
    output_v << "  </UnstructuredGrid>" << std::endl
             << "</VTKFile>" << std::endl;



}

//------------------------------------------------------------------------------
//----------------COMPUTE AND ASSEMBLE THE GLOBAL MATRIX AND VECTOR-------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Fluid<DIM,DEG>::computeError(VecDouble &errorsTotal) {

    std::ofstream rprint("errors.txt",std::ios::app);

    VecDouble errorsProcess(3);
    errorsTotal.setZero();
    errorsProcess.setZero();

    // Loop over the elements
    for (int jel = numElem; jel--; ){

        VecDouble errors(3);

        if (fProbType == ProblemType::EPoisson) elements_[jel] -> computeErrorPoisson(errors);
        
        errorsProcess += errors;

    }; //Elements

    MPI_Allreduce(&errorsProcess[0],&errorsTotal[0],errorsTotal.size(),MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

    if (rank == 0){
        std::cout << "\n\nERROR REPORT:\n" << std::scientific << std::setprecision(10)
            << "L2 state var = " << sqrt(errorsTotal[0]) << "\n" 
            << "Semi H1 state var = " << sqrt(errorsTotal[1]) << "\n" 
            << "H1 state var = " << sqrt(errorsTotal[2]) << "\n"; 
        rprint << sqrt(errorsTotal[0]) << " " << sqrt(errorsTotal[1]) << " " << sqrt(errorsTotal[2]) << std::endl;
    }
    return;
}

template class Fluid<2,1>;
template class Fluid<2,2>;
template class Fluid<2,3>;
template class Fluid<3,1>;
template class Fluid<3,2>;
template class Fluid<3,3>;