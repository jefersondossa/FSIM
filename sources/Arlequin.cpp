#include "Arlequin.h"
#include "hdf5.h"
#include "ElCouplingLocal.h"
#include "ElCouplingGlobal.h"
#include <set>

#define  REAL double

// External includes
#ifndef TRILIBRARY
#define TRILIBRARY
#endif

#include "triangle.h"

extern "C"
{
	void triangulate(char *, struct triangulateio *, struct triangulateio *,struct triangulateio *);
	void trifree(void *);
}

struct PointTriangle {
    double x, y;
    PointTriangle()=default;
    PointTriangle(double x, double y) : x(x), y(y) {}
    PointTriangle(VecDouble &coord) : x(coord[0]), y(coord[1]) {}
};

struct LineTriangle {
    PointTriangle start, end;

    LineTriangle(const PointTriangle& start, const PointTriangle& end) : start(start), end(end) {}
};

// Helper function to calculate the cross product of two vectors.
double CrossProduct(const PointTriangle& v1, const PointTriangle& v2) {
    return v1.x * v2.y - v1.y * v2.x;
}

PointTriangle Sum(const PointTriangle& a, const PointTriangle& b){
    PointTriangle res;
    res.x = a.x + b.x;
    res.y = a.y + b.y;
    return res;
}

PointTriangle SumScalar(const PointTriangle& a, const PointTriangle& b, double scal){
    PointTriangle res;
    res.x = a.x + scal*b.x;
    res.y = a.y + scal*b.y;
    return res;
}

PointTriangle Dif(const PointTriangle& a, const PointTriangle& b){
    PointTriangle res;
    res.x = a.x - b.x;
    res.y = a.y - b.y;
    return res;
}

// Check if point p is inside the triangle defined by a, b, and c.
bool IsPointInsideTriangle(const PointTriangle& p, const PointTriangle& a, const PointTriangle& b, const PointTriangle& c) {
    PointTriangle v0 = Dif(c,a);
    PointTriangle v1 = Dif(b,a);
    PointTriangle v2 = Dif(p,a);

    double dot00 = v0.x * v0.x + v0.y * v0.y;
    double dot01 = v0.x * v1.x + v0.y * v1.y;
    double dot02 = v0.x * v2.x + v0.y * v2.y;
    double dot11 = v1.x * v1.x + v1.y * v1.y;
    double dot12 = v1.x * v2.x + v1.y * v2.y;

    // Compute barycentric coordinates
    double invDenom = 1. / (dot00 * dot11 - dot01 * dot01);
    double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    // Check if the point is inside the triangle
    return (u >= 0) && (v >= 0) && (u + v <= 1);
}

// Calculate the intersection point of two line segments.
PointTriangle LineIntersectionPoint(const LineTriangle& line1, const LineTriangle& line2) {
    PointTriangle dir1 = Dif(line1.end,line1.start);
    PointTriangle dir2 = Dif(line2.end,line2.start);
    PointTriangle start1_to_start2 = Dif(line2.start,line1.start);

    double denominator = CrossProduct(dir1, dir2);

    if (denominator == 0) {
        return PointTriangle(0, 0); // Lines are parallel or collinear, no intersection.
    }

    double t1 = CrossProduct(start1_to_start2, dir2) / denominator;
    double t2 = CrossProduct(start1_to_start2, dir1) / denominator;

    if (t1 >= 0.0 && t1 <= 1.0 && t2 >= 0.0 && t2 <= 1.0) {
        return SumScalar(line1.start,dir1,t1);
    }

    return PointTriangle(0, 0); // No intersection.
}

// Check if all nodes of one triangle are inside the other triangle.
bool AreTriangleNodesInside(const std::vector<PointTriangle>& triangle1, const std::vector<PointTriangle>& triangle2) {
    for (const PointTriangle& node : triangle1) {
        if (!IsPointInsideTriangle(node, triangle2[0], triangle2[1], triangle2[2])) {
            return false;
        }
    }
    return true;
}

// Find nodes of one triangle that are inside the other triangle.
std::vector<PointTriangle> NodesInsideTriangle(const std::vector<PointTriangle>& triangle1, const std::vector<PointTriangle>& triangle2) {
    std::vector<PointTriangle> nodesInside;

    for (const PointTriangle& node : triangle1) {
        if (IsPointInsideTriangle(node, triangle2[0], triangle2[1], triangle2[2])) {
            nodesInside.push_back(node);
        }
    }

    return nodesInside;
}

std::vector<PointTriangle> TriangleIntersectionPoints(const std::vector<PointTriangle>& triangle1, const std::vector<PointTriangle>& triangle2) {
    std::vector<PointTriangle> intersectionPoints;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            LineTriangle edge1(triangle1[i], triangle1[(i + 1) % 3]);
            LineTriangle edge2(triangle2[j], triangle2[(j + 1) % 3]);
            PointTriangle intersection = LineIntersectionPoint(edge1, edge2);

            if (intersection.x != 0 || intersection.y != 0) {
                if (IsPointInsideTriangle(intersection, triangle1[0], triangle1[1], triangle1[2]) &&
                    IsPointInsideTriangle(intersection, triangle2[0], triangle2[1], triangle2[2])) {
                    intersectionPoints.push_back(intersection);
                }
            }
        }
    }
    return intersectionPoints;
}







//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------
void Arlequin::ComputeTriangleIntersections(int64_t iEl, std::set<int64_t> &elIntersected){
    auto refElement = fMeshVector[2]->ElementVec()[iEl];
    auto connectref = fMeshVector[0]->ElementVec()[iEl]->getConnectivity();
    std::set<int> pointsInsideRefElement;
    
    for (auto interCoarse : elIntersected){
        auto intElement = fMeshVector[0]->ElementVec()[interCoarse];
        auto connectint=intElement->getConnectivity();

        PointTriangle ref1(fMeshVector[2]->NodeVec()[connectref[0]]->getCoordinates());
        PointTriangle ref2(fMeshVector[2]->NodeVec()[connectref[1]]->getCoordinates());
        PointTriangle ref3(fMeshVector[2]->NodeVec()[connectref[2]]->getCoordinates());
        PointTriangle int1(fMeshVector[2]->NodeVec()[connectint[0]]->getCoordinates());
        PointTriangle int2(fMeshVector[2]->NodeVec()[connectint[1]]->getCoordinates());
        PointTriangle int3(fMeshVector[2]->NodeVec()[connectint[2]]->getCoordinates());

        std::vector<PointTriangle> reftriangle = {ref1,ref2,ref3};
        std::vector<PointTriangle> inttriangle = {int1,int2,int3};

        // Find intersection points.
        std::vector<PointTriangle> intersectionPoints = TriangleIntersectionPoints(reftriangle, inttriangle);

        // Find nodes of Triangle 1 that are inside Triangle 2.
        std::vector<PointTriangle> nodesIn = NodesInsideTriangle(reftriangle, inttriangle);

        // if (nodesIn.size()>1) {
        //     std::cout << "The case of more than one point inside other element is\n";
        //     PanicButton();
        // }
        //Creating the containers for the input and output
        struct triangulateio in;
        struct triangulateio out;
        // clearTrianglesList(out);

        // buildInput(nodes, param, in);

        // int triangle_error = generateTesselation(in, out);
        
        std::string triflags = "znQP";
        char *triswitches = new char[triflags.size()+1];
        std::strcpy(triswitches, triflags.c_str());

        int triangle_error = 0;

        // try
        // {
        //     triangulate(triswitches, &in, &out, (struct triangulateio *)NULL);
        // }



        // setToContainer(out);

        // executePostMeshingProcesses(nodes, elements, param);

        // deleteInContainer(in);
        // deleteOutContainer(out);

    }
}




void Arlequin::CreateGlobalCouplingElements(){
    //The main idea of this code is: for each integration point in the local model,
    //we look for the interesected global element. Then, for each intersected global
    // element we create an ElCouplingGlobal associated with the corresponding local element.
    int numberIntPoints = fMeshVector[2]->ElementVec()[0] -> getNumberOfIntegrationPoints();
    int64_t nEl = fMeshVector[2]->NElements();
    int64_t index = nEl; 
    for (int iel = 0; iel<nEl; iel++){
        auto jel = fMeshVector[2]->ElementVec()[iel];
        std::set<int64_t> elIntersected;
        for (int i=0; i<numberIntPoints; i++){
            elIntersected.insert(fLocalIntPointToGlobalElement[iel][i]);
        };
        // int nElIntersected = elIntersected.size();
        if (elIntersected.size() > 1){
            ComputeTriangleIntersections(iel, elIntersected);
        } else {
        for (auto ielcoarse : elIntersected){
            ElCouplingGlobal *el = new ElCouplingGlobal(index++,ielcoarse,fMeshVector);
            fGlobalElToLocalEl[ielcoarse].insert(index-1);
            fMeshVector[2]->ElementVec().push_back(el);
            el->setConnectivity(fMeshVector[2]->ElementVec()[iel]->getConnectivity());
            el->SetGlobalXsi(fLocalIntPointToGlobalXsi[iel]);
            el->SetGlobalElemCorresp(fLocalIntPointToGlobalElement[iel]);
        }
        }
    }
    
    std::cout <<" Testsd " << std::endl;

}



//------------------------------------------------------------------------------
//-------------------------COMPUTE ELEMENT REGIONS/BOXES------------------------
//------------------------------------------------------------------------------

void Arlequin::SetElementBoxes() {
    
    int DIM = fMeshVector[0]->Dimension();
    VecInt connec;
    VecDouble xk(DIM), Xk(DIM);

    //Compute element boxes for coarse model
    //Only function for straight elements
    for (int jel = 0; jel < fMeshVector[0]->NElements(); jel++){
        connec = fMeshVector[0]->ElementVec()[jel] -> getConnectivity();
        VecDouble x1 = fMeshVector[0]->NodeVec()[connec[0]] -> getCoordinates();
        VecDouble x2 = fMeshVector[0]->NodeVec()[connec[1]] -> getCoordinates();
        VecDouble x3 = fMeshVector[0]->NodeVec()[connec[2]] -> getCoordinates();      

        xk[0] = std::min(x1[0],std::min(x2[0], x3[0]));
        xk[1] = std::min(x1[1],std::min(x2[1], x3[1]));

        Xk[0] = std::max(x1[0],std::max(x2[0], x3[0]));
        Xk[1] = std::max(x1[1],std::max(x2[1], x3[1]));        
        
        fMeshVector[0]->ElementVec()[jel] -> setIntersectionParameters(xk, Xk);
    };

    //Compute element boxes for fine model
    //Only function for straight elements
    for (int jel = 0; jel < fMeshVector[1]->NElements(); jel++){
        connec = fMeshVector[1]->ElementVec()[jel] -> getConnectivity();
        VecDouble x1 = fMeshVector[1]->NodeVec()[connec[0]] -> getCoordinates();
        VecDouble x2 = fMeshVector[1]->NodeVec()[connec[1]] -> getCoordinates();
        VecDouble x3 = fMeshVector[1]->NodeVec()[connec[2]] -> getCoordinates();      

        xk[0] = std::min(x1[0],std::min(x2[0], x3[0]));
        xk[1] = std::min(x1[1],std::min(x2[1], x3[1]));

        Xk[0] = std::max(x1[0],std::max(x2[0], x3[0]));
        Xk[1] = std::max(x1[1],std::max(x2[1], x3[1]));        
        
        fMeshVector[1]->ElementVec()[jel] -> setIntersectionParameters(xk, Xk);
    };

    return;
};


//------------------------------------------------------------------------------
//-------------------COMPUTE NODAL CORRESPONDECE WITH ELEMENTS------------------
//------------------------------------------------------------------------------

void Arlequin::searchNodeCorrespondence(VecDouble &x,CompMesh *cmesh, 
                                        int &elCorr, VecDouble &xsiCorr, int elSearch){
    
    int DIM = cmesh->Dimension();
    int DEG = cmesh->GetDefaultOrder();
    ShapeFunction shapeQuad(DIM,DEG);
    int nElNodes = cmesh->NElNodes();
    VecDouble phi_(nElNodes);    
    VecDouble xsiCC(3);
    std::pair<VecDouble,VecDouble> XK;

    elCorr = 150000;
    VecDouble xsi(DIM);
    VecDouble x_(DIM);
    VecDouble deltaX(DIM);
    VecDouble deltaXsi(DIM);
    xsi.setZero(); x_.setZero(); deltaX.setZero(); deltaXsi.setZero();
    bool flg = true;
    VecInt connec = cmesh->ElementVec()[elSearch] -> getConnectivity();
    
    xsiCC.fill(1.e10);
    xsiCorr.fill(1.e50);
    xsi.fill(1./3.);
    x_.fill(0.);
    
    shapeQuad.evaluate(xsi,phi_);

    for (int i = 0; i < nElNodes; i++){
        VecDouble xint = cmesh->NodeVec()[connec[i]] -> getCoordinates();
        for (int k = 0; k < DIM; k++){
            x_[k] += xint[k] * phi_[i];
        }        
    };

    double error = 1.e6;
    int iterations = 0;

    while ((error > 1.e-8) && (iterations < 4)) {
        
        iterations++;
        
        for (int k = 0; k < DIM; k++) deltaX[k] = x[k] - x_[k];
        deltaXsi.setZero();
        
        cmesh->ElementVec()[elSearch] -> ComputeJacobian(0);

        // for (int i = 0; i < DIM; i++)
        //     for (int j = 0; j < DIM; j++)
        //         deltaXsi[i] += ainv(j,i) * deltaX[j];
        deltaXsi = cmesh->ElementVec()[elSearch]->IntegrationData().fA0Inv.transpose()*deltaX;    

        xsi += deltaXsi;
        x_.setZero();
        
        shapeQuad.evaluate(xsi,phi_);
        
        for (int i=0; i<nElNodes; i++){
            VecDouble xint = cmesh->NodeVec()[connec[i]] -> getCoordinates();
            for (int k = 0; k < DIM; k++)x_[k] += xint[k] * phi_[i];
        };

        error = std::sqrt(deltaXsi[0]*deltaXsi[0] + deltaXsi[1]*deltaXsi[1]);
    };
    
    double t1 = -1.e-2;
    double t2 =  1. - t1;
    
    xsiCC[0] = xsi[0];
    xsiCC[1] = xsi[1];       
    xsiCC[2] = 1. - xsiCC[0] - xsiCC[1];

    if ((xsiCC[0] >= t1) && (xsiCC[1] >= t1) && (xsiCC[2] >= t1) &&
        (xsiCC[0] <= t2) && (xsiCC[1] <= t2) && (xsiCC[2] <= t2)){

        xsiCorr[0] = xsi[0]; xsiCorr[1] = xsi[1];
        elCorr = elSearch;
        // return;
    } else {

        int nEl;
        if (fMeshVector[1]->getProblemParameters().getTimeInstant() == 0){
            nEl = cmesh->NElements();
        } else {
            nEl = cmesh->ElementVec()[elSearch] -> getNumberOfNeighborElements(); 
        }   

        for (int jel = 0; jel < nEl; jel++){

            if (fMeshVector[1]->getProblemParameters().getTimeInstant() == 0){
                connec = cmesh->ElementVec()[jel] -> getConnectivity();
            } else {
                connec = cmesh->ElementVec()[cmesh->ElementVec()[elSearch] -> getNeighborElement(jel)] -> getConnectivity();
            }   

            //get boxes information        
            XK = cmesh->ElementVec()[jel] -> getXIntersectionParameter();

            //Chech if the node is inside the element box
            if ((x[0] < XK.first[0]*0.95) || (x[0] > XK.second[0]*1.05) ||
                (x[1] < XK.first[1]*0.95) || (x[1] > XK.second[1]*1.05)) continue;
            
            //Compute nodal correspondence
            for (int i = DIM+1; i--; ) xsiCC[i] = 1.e10;
    
            for (int i = DIM; i--; ){
                xsi[i] = 1. / 3.;
                x_[i] = 0.;
            }

            shapeQuad.evaluate(xsi,phi_);

            for (int i = 0; i < nElNodes; i++){
                VecDouble xint = cmesh->NodeVec()[connec[i]] -> getCoordinates();
                for (int k = DIM; k--; )
                    x_[k] += xint[k] * phi_[i];
            };
            
            double error = 1.e6;
            int iterations = 0;

            while ((error > 1.e-8) && (iterations < 4)) {
                
                iterations++;

                for (int k = DIM; k--; ){
                    deltaX[k] = x[k] - x_[k];
                    deltaXsi[k] = 0.;
                }
                
                cmesh->ElementVec()[jel] -> ComputeJacobian(0);
                auto ainv = cmesh->ElementVec()[jel]->IntegrationData().fA0Inv;
                for (int i = 0; i < DIM; i++)
                    for (int j = 0; j < DIM; j++)
                        deltaXsi[i] += ainv(j,i) * deltaX[j];
                
                for (int k = DIM; k--; ){
                    xsi[k] += deltaXsi[k];
                    x_[k] = 0.;
                }
            
                shapeQuad.evaluate(xsi,phi_);
                
                for (int i=0; i<nElNodes; i++){
                    VecDouble xint = cmesh->NodeVec()[connec[i]] -> getCoordinates();
                    for (int k = DIM; k--; ) x_[k] += xint[k] * phi_[i];
                };
                        
                error = std::sqrt(deltaXsi[0]*deltaXsi[0] + deltaXsi[1]*deltaXsi[1]);
            };
            
            double t1 = -1.e-2;
            double t2 =  1. - t1;
            
            xsiCC[0] = xsi[0];
            xsiCC[1] = xsi[1];       
            xsiCC[2] = 1. - xsiCC[0] - xsiCC[1];

            if ((xsiCC[0] >= t1) && (xsiCC[1] >= t1) && (xsiCC[2] >= t1) &&
                (xsiCC[0] <= t2) && (xsiCC[1] <= t2) && (xsiCC[2] <= t2)){

                xsiCorr[0] = xsi[0]; xsiCorr[1] = xsi[1];
                elCorr = jel;
                break;
            }
        }
    };

    if (fabs(xsi[0]) > 2.) {
        std::cout << "PROBLEM SEARCHING NODE CORRESPONDENCE " << std::endl;
        PanicButton();
    }
    return;
};

//------------------------------------------------------------------------------
//--------COMPUTE NODAL CORRESPONDECE OF FINE NODES WITH COARSE ELEMENTS--------
//------------------------------------------------------------------------------

void Arlequin::setNodalCorrespondenceFine() {

    int DIM = fMeshVector[0]->Dimension();
    int DEG = fMeshVector[0]->GetDefaultOrder();
    //FINE MESH
    double &alpha_f = fMeshVector[0]->getProblemParameters().getAlphaF();

    for (int inode = 0; inode < fMeshVector[2]->NNodes(); inode++) {
        
        VecDouble x = fMeshVector[2]->NodeVec()[inode] -> getCoordinates();

        int elCorr = 0;
        VecDouble xsiCorr(DIM);
        searchNodeCorrespondence(x, fMeshVector[0],elCorr,xsiCorr,fNodeLocalToElementGlobal[inode]);
        fNodeLocalToElementGlobal[inode]=elCorr;
        fNodeLocalToXsiGlobal[inode]=xsiCorr;
        // fMeshVector[2]->NodeVec()[inode] -> setNodalCorrespondence(elCorr,xsiCorr);             
        // }

        // std::cout << "INODE " << nodesGlueZoneFine_[inode] << " " << elCorr << " " << xsiCorr[0] << " " << xsiCorr[1] << std::endl;
   
        // std::cout << "CORRESP " << elCorr << " " << corresp.first << std::endl 
                  // << xsiCorr[0] << " " << xsiCorr[1] << " " << corresp.second[0] << " " << corresp.second[1] << std::endl;

        
            
        // std::cout << "corresp " << corresp.first 
        //           << " " << corresp.second(0) << " " << corresp.second(1)
        //           << std::endl;

    };

    // for (int i=0; i<fMeshVector[1]->NNodes(); i++){
    //     int elem = fMeshVector[1]->NodeVec()[i] -> getNodalElemCorrespondence();
    //     xsi = fMeshVector[1]->NodeVec()[i] -> getNodalXsiCorrespondence();
    //     std::cout << "node " << i << " elem " << elem 
    //               << " " << xsi(0) << " " << xsi(1) << std::endl;
    // };

    //Compute correspondence of integration points
    //int numberIntPoints = fMeshVector[1]->ElementVec()[0] -> getNumberOfIntegrationPoints();
    //if (rank == 0) std::cout << "Int Points " << numberIntPoints << std::endl;

    for (int ielem = 0; ielem < fMeshVector[2]->NElements(); ielem++) {
        int nElNodes = fMeshVector[2]->NElNodes();
        VecDouble x1(nElNodes), x2(nElNodes), x3(nElNodes);
        VecInt connec = fMeshVector[2]->ElementVec()[ielem] -> getConnectivity();
        
        for (int i = 0; i < nElNodes; i++){
            VecDouble x = fMeshVector[2]->NodeVec()[connec[i]] -> getCoordinates();
            VecDouble xp = fMeshVector[2]->NodeVec()[connec[i]] -> getPreviousCoordinates();
            
            x1[i] = alpha_f * x[0] + (1. - alpha_f) * xp[0];
            x2[i] = alpha_f * x[1] + (1. - alpha_f) * xp[1];
            if (DIM == 3) x3[i] = alpha_f * x[2] + (1. - alpha_f) * xp[2];
        };

        // std::cout << "XX1 " << x2 << " " << x22 << " " << x222 << std::endl;

        int numberIntPoints = fMeshVector[2]->ElementVec()[ielem] -> getNumberOfIntegrationPoints();

        IntegQuadratureSpecial quad(DIM,DEG);
        for (int i = 0; i < numberIntPoints; i++){

            VecDouble x(DIM);
            
            x[0] = quad.interpolateQuadraticVariable(x1,i,DIM);
            x[1] = quad.interpolateQuadraticVariable(x2,i,DIM);

            int elCorr = 0;
            VecDouble xsiCorr(DIM);

            if (fLocalIntPointToGlobalElement[ielem].size()==0) fLocalIntPointToGlobalElement[ielem].resize(numberIntPoints);
            if (fLocalIntPointToGlobalXsi[ielem].rows()==0) fLocalIntPointToGlobalXsi[ielem].resize(numberIntPoints,DIM);
            searchNodeCorrespondence(x,fMeshVector[0],elCorr,xsiCorr,
                                     fLocalIntPointToGlobalElement[ielem][i]);
                            
            fLocalIntPointToGlobalElement[ielem][i] = elCorr;
            fLocalIntPointToGlobalXsi[ielem](i,0)=xsiCorr[0];
            fLocalIntPointToGlobalXsi[ielem](i,1)=xsiCorr[1];
            if (DIM == 3)fLocalIntPointToGlobalXsi[ielem](i,2)=xsiCorr[2];
        };

    }; 
};

//------------------------------------------------------------------------------
//--------------------SETS THE COUPLING ZONE IN COARSE MODEL--------------------
//------------------------------------------------------------------------------

void Arlequin::setSignaledDistance(){
    
    std::map<int64_t,VecDouble> NodalNormalVector;

    int DIM = fMeshVector[0]->Dimension();
    int DEG = fMeshVector[0]->GetDefaultOrder();
    int nBdNodes = fMeshVector[1]->NBdNodes();
    int bconnec[nBdNodes];
    double dist;
    //approximate normal calculation

    for (int i = 0; i < fMeshVector[1]->NBoundElements(); i++){
        if (fMeshVector[1]->BoundaryVec()[i]->getConstrain(0) == 2){
            auto iel =fMeshVector[1]->BoundaryVec()[i]->getElement();
            auto nelements = fMeshVector[1]->NElements();
            VecInt connec = fMeshVector[1]->ElementVec()[fMeshVector[1]->BoundaryVec()[i]->getElement()] -> getConnectivity();
            
            fMeshVector[1]->ElementVec()[fMeshVector[1]->BoundaryVec()[i] -> getElement()] -> getBoundaryNodes(bconnec);
            
            //Loop over the 1D element segments
            for (int iSeg = 0; iSeg < DEG; iSeg++){
                int no1,no2;
                if (iSeg == 0){
                    no1 = bconnec[0];
                    if (DEG == 1) {
                        no2 = bconnec[1];
                    } else {
                        no2 = bconnec[2];
                    }
                } else {
                    no1 = bconnec[iSeg+1];
                    if (DEG == 2 || iSeg == 2){
                        no2 = bconnec[1];
                    } else {
                        no2 = bconnec[3];
                    }
                }
                
                VecDouble x1 = fMeshVector[1]->NodeVec()[no1] -> getCoordinates();
                VecDouble x2 = fMeshVector[1]->NodeVec()[no2] -> getCoordinates();

                double sLength = sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                                      (x1[0] - x2[0]) * (x1[0] - x2[0]));

                VecDouble n(2);
                n[0] = (x2[1] - x1[1]) / sLength;
                n[1] = (x1[0] - x2[0]) / sLength;

                NodalNormalVector[no1] = n;
                NodalNormalVector[no2] = n;
                // fMeshVector[1]->NodeVec()[no1] -> setInnerNormal(n);
                // fMeshVector[1]->NodeVec()[no2] -> setInnerNormal(n);
            }
        };
    };


    //Fine mesh nodes
    for (int ino = 0; ino < fMeshVector[1]->NNodes(); ino++){
        VecDouble x = fMeshVector[1]->NodeVec()[ino] -> getCoordinates();
        dist=10000000000000000000000000000.;
        
        for (int i = 0; i < fMeshVector[1]->NBoundElements(); i++){
            if (fMeshVector[1]->BoundaryVec()[i] -> getConstrain(0) == 2){
                VecInt connec = fMeshVector[1]->ElementVec()[fMeshVector[1]->BoundaryVec()[i] -> getElement()] -> getConnectivity();
                
                fMeshVector[1]->ElementVec()[fMeshVector[1]->BoundaryVec()[i] -> getElement()] -> getBoundaryNodes(bconnec);

                

                for (int iSeg = 0; iSeg < DEG ; ++iSeg){
                    
                
                    //first segment
                    int no1,no2;
                    if (iSeg == 0){
                        no1 = bconnec[0];
                        if (DEG == 1) {
                            no2 = bconnec[1];
                        } else {
                            no2 = bconnec[2];
                        }
                    } else {
                        no1 = bconnec[iSeg+1];
                        if (DEG == 2 || iSeg == 2){
                            no2 = bconnec[1];
                        } else {
                            no2 = bconnec[3];
                        }
                    }
                  
                    VecDouble x1 = fMeshVector[1]->NodeVec()[no1] -> getCoordinates();
                    VecDouble x2 = fMeshVector[1]->NodeVec()[no2] -> getCoordinates();
                    
                    double aux0 =  sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                                        (x2[0] - x1[0]) * (x2[0] - x1[0]));
                    double aux1 = ((x[0] - x1[0]) * (x2[0] - x1[0])+
                                   (x[1] - x1[1]) * (x2[1] - x1[1])) / aux0;
                    double dist2 =-((x2[1] - x1[1]) * x[0] - 
                                    (x2[0] - x1[0]) * x[1] +
                                    x2[0] * x1[1] - x2[1] * x1[0]) / aux0;
                    
                    if (aux1 > aux0){
                        dist2 = sqrt((x2[1] - x[1]) * (x2[1] - x[1]) +
                                     (x2[0] - x[0]) * (x2[0] - x[0]));
                        //find signal
                        //side normal vector
                        VecDouble n = NodalNormalVector[no2];
                        double test[2];

                        test[0] = x[0] - x2[0];
                        test[1] = x[1] - x2[1];
                        double signaltest = n[0]*test[0] + n[1]*test[1];
                        double signal = -1.;
                        
                        if (signaltest <= -0.001)signal = 1.;
                        
                        dist2 *= signal;
                    };

                    if (aux1 < 0.){
                        dist2 = sqrt((x[1] - x1[1]) * (x[1] - x1[1]) +
                                     (x[0] - x1[0]) * (x[0] - x1[0]));
                        //find signal
                        //side normal vector
                        VecDouble n = NodalNormalVector[no1];
                        double test[2];

                        test[0] = x[0] - x1[0];
                        test[1] = x[1] - x1[1];
                        double signaltest = n[0]*test[0] + n[1]*test[1];
                        double signal = -1.;
                        
                        if (signaltest <= -0.001) signal = 1.;
                        
                        dist2 *= signal;
                    };
                    
                    if (fabs(dist2) < fabs(dist)) dist = dist2;
                }
            }; //if bf is the glue boundary
        }; //
        if(dist < 0) dist = 0;
        dist = x[0] - 0.5;
        fLocalSignaledDistance[ino] = dist;
     };

    //Coarse mesh
     for (int ino = 0; ino < fMeshVector[0]->NNodes(); ino++){
        VecDouble x = fMeshVector[0]->NodeVec()[ino] -> getCoordinates();
        dist=10000000000000000000000000000.;
        
        for (int i = 0; i < fMeshVector[1]->NBoundElements(); i++){
            if (fMeshVector[1]->BoundaryVec()[i] -> getConstrain(0) == 2){

                VecInt connec = fMeshVector[1]->ElementVec()[fMeshVector[1]->BoundaryVec()[i] -> getElement()] -> getConnectivity();
                fMeshVector[1]->ElementVec()[fMeshVector[1]->BoundaryVec()[i] -> getElement()] -> getBoundaryNodes(bconnec);

                for (int iSeg = 0; iSeg < DEG ; ++iSeg){
                    //first segment
                    int no1,no2;
                    if (iSeg == 0){
                        no1 = bconnec[0];
                        if (DEG == 1) {
                            no2 = bconnec[1];
                        } else {
                            no2 = bconnec[2];
                        }
                    } else {
                        no1 = bconnec[iSeg+1];
                        if (DEG == 2 || iSeg == 2){
                            no2 = bconnec[1];
                        } else {
                            no2 = bconnec[3];
                        }
                    }
                    // std::cout<<no1<<" nos "<<no2<<std::endl;
                  
                    VecDouble x1 = fMeshVector[1]->NodeVec()[no1] -> getCoordinates();
                    VecDouble x2 = fMeshVector[1]->NodeVec()[no2] -> getCoordinates();
                    
                    double aux0 =  sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                                        (x2[0] - x1[0]) * (x2[0] - x1[0]));
                    double aux1 = ((x[0] - x1[0]) * (x2[0] - x1[0])+
                                   (x[1] - x1[1]) * (x2[1] - x1[1])) / aux0;
                    double dist2 =-((x2[1] - x1[1]) * x[0] - 
                                    (x2[0] - x1[0]) * x[1] +
                                    x2[0] * x1[1] - x2[1] * x1[0]) / aux0;
                    
                    if (aux1 > aux0){
                        dist2 = sqrt((x2[1] - x[1]) * (x2[1] - x[1]) +
                                     (x2[0] - x[0]) * (x2[0] - x[0]));
                        //find signal
                        //side normal vector
                        VecDouble n = NodalNormalVector[no2];
                        
                        double test[2];
                        test[0] = x[0] - x2[0];
                        test[1] = x[1] - x2[1];
                        double signaltest = n[0]*test[0] + n[1]*test[1];
                        double signal = -1.;
                        
                        if (signaltest <= -0.001)signal = 1.;
                        
                        dist2 *= signal;
                    };

                    if (aux1 < 0.){
                        dist2 = sqrt((x[1] - x1[1]) * (x[1] - x1[1]) +
                                     (x[0] - x1[0]) * (x[0] - x1[0]));
                        //find signal
                        //side normal vector
                        VecDouble n = NodalNormalVector[no1];
                        
                        double test[2];
                        test[0] = x[0] - x1[0];
                        test[1] = x[1] - x1[1];
                        double signaltest = n[0]*test[0] + n[1]*test[1];
                        double signal = -1.;
                        
                        if (signaltest <= -0.001) signal = 1.;
                        
                        dist2 *= signal;
                    };
                    
                    if (fabs(dist2) < fabs(dist)) dist = dist2;
                }
            }; //if bf is the glue boundary
        }; //
    
        if (fabs(fGlobalSignaledDistance[ino]) < 1.e-2){
            fGlobalSignaledDistance[ino] = dist;
        };
        dist = x[0] - 0.5;
        fGlobalSignaledDistance[ino] = dist;
     };

    for (int jel = 0; jel < fMeshVector[0]->NElements(); jel++){
        auto connect = fMeshVector[0]->ElementVec()[jel]->getConnectivity();
        VecDouble distfunction(connect.size());
        for (int i = 0; i < connect.size(); i++){
            distfunction[i] = fGlobalSignaledDistance[connect[i]];
        }
        fMeshVector[0]->ElementVec()[jel] -> ComputeIntPointDistFunction(distfunction);        
    };
    for (int jel = 0; jel < fMeshVector[1]->NElements(); jel++){
        auto connect = fMeshVector[1]->ElementVec()[jel]->getConnectivity();
        VecDouble distfunction(connect.size());
        for (int i = 0; i < connect.size(); i++){
            distfunction[i] = fLocalSignaledDistance[connect[i]];
        }
        fMeshVector[1]->ElementVec()[jel] -> ComputeIntPointDistFunction(distfunction);        
    };
    // for (int jel = 0; jel < fMeshVector[1]->NElements(); jel++){
    //     fMeshVector[1]->ElementVec()[jel] -> ComputeIntPointDistFunction();        
    // };

};

//------------------------------------------------------------------------------
//--------------------SETS THE COUPLING ZONE IN COARSE MODEL--------------------
//------------------------------------------------------------------------------

void Arlequin::setCouplingZone(){

    // double dist;
    int flag;
    int nodesCZ[fMeshVector[1]->NNodes()];
    int nodesCZ2[fMeshVector[0]->NNodes()];

    double lim1 = 0.06251;
    double lim2 = 0.93749;
    double tick = 0.01249;

    for (int i = 0; i < fMeshVector[1]->NNodes(); i++) nodesCZ[i] = 0;    
    int index = 0; 

    //Defines a criterion to select the elements that are in the glue zone
    for (int jel = 0; jel < fMeshVector[1]->NElements(); jel++){
        
        VecInt connec = fMeshVector[1]->ElementVec()[jel] -> getConnectivity();
        flag = 0;
        int nElNodes = fMeshVector[1]->NElNodes();
        for (int ino = 0; ino < nElNodes; ino++){
            VecDouble x = fMeshVector[1]->NodeVec()[connec[ino]] -> getCoordinates();
            double dist = fLocalSignaledDistance[connec[ino]];
            //  std::cout << "DIST " << dist << std::endl;
            if (dist <= fGlueZoneThickness + 0.001) flag++;
            
        };

        if (flag == nElNodes) {
            // elementsGlueZoneFine_.push_back(jel);
            // fMeshVector[1]->ElementVec()[jel] -> setGlueZone();

            ElCouplingLocal *el = new ElCouplingLocal(index++,jel,fMeshVector);
            fMeshVector[2]->ElementVec().push_back(el);            
        };        
    };

    //Defines which nodes are in the glue zone
    numElemGlueZoneFine = fMeshVector[2]->NElements();
    for (int i = 0; i < numElemGlueZoneFine; i++){
        auto *elcoup = dynamic_cast<ElCouplingLocal *> (fMeshVector[2]->ElementVec()[i]);
        auto locindex = elcoup->GetLocalIndex();
        VecInt connec = fMeshVector[1]->ElementVec()[locindex] -> getConnectivity();
        int nElNodes = fMeshVector[1]->NElNodes();
        for (int ino = 0; ino < nElNodes; ino++) nodesCZ[connec[ino]]++;
    };

    //Compute number of nodes in the glue zone
    // numNodesGlueZoneFine = 0;
    std::map<int,int> FineToGluing;
    int64_t nodeindex = 0;
    for (int i = 0; i < fMeshVector[1]->NNodes(); i++){
        if(nodesCZ[i] > 0) {
            // numNodesGlueZoneFine++;
            auto x = fMeshVector[1]->NodeVec()[i]->getCoordinates();
            Node *node = new Node(x,nodeindex,fMeshVector[1]->NState());
            fMeshVector[2]->NodeVec().push_back(node);
            FineToGluing[i] = nodeindex;
            nodeindex++;
        };
    };

    if (rank == 0) std::cout << "GLUE ZONE - Number of Nodes = " 
                             << fMeshVector[2]->NNodes() 
                             << " - Number of Elements = " 
                             << fMeshVector[2]->NElements() << std::endl;
    

    for (int i = 0; i < fMeshVector[2]->NElements(); i++){
        int nElNodes = fMeshVector[1]->NElNodes();
        VecInt connecAux(nElNodes);
        auto *elcoup = dynamic_cast<ElCouplingLocal *> (fMeshVector[2]->ElementVec()[i]);
        auto locindex = elcoup->GetLocalIndex();

        VecInt connec = fMeshVector[1]->ElementVec()[locindex] -> getConnectivity();
        
        // for (int ino = 0; ino < numNodesGlueZoneFine; ino++)
        for (int k = 0; k < nElNodes; k++) connecAux[k] = FineToGluing[connec[k]];
            // if (nodesGlueZoneFine_[ino] == connec[k]) connecAux[k] = ino;

        fMeshVector[2]->ElementVec()[i] -> setConnectivity(connecAux);
        // glueZoneFine_[i] -> setNodes(nodesLagrangeFine_);

    };

    fMeshVector[2]->NGlobalDOF() = fMeshVector[2]->NNodes() * fMeshVector[2]->NState();  

    // // Glue Zone in coarse mesh
    // for (int i = 0; i < fMeshVector[0]->NNodes(); i++) nodesCZ2[i] = 0;    

    // elementsGlueZoneCoarse_.reserve(fMeshVector[0]->NElements() / 3);
    // nodesGlueZoneCoarse_.reserve(fMeshVector[0]->NNodes() / 3);

    // //Defines a criterion to select the elements that are in the glue zone
    // for (int jel = 0; jel < fMeshVector[0]->NElements(); jel++){
    //     int nElNodes = fMeshVector[0]->NElNodes();
    //     VecInt connec = fMeshVector[0]->ElementVec()[jel] -> getConnectivity();
    //     flag = 0;

    //     for (int ino = 0; ino < nElNodes; ino++){
    //         VecDouble x = fMeshVector[0]->NodeVec()[connec[ino]] -> getCoordinates();

    //         if ((x[0] < lim1) || (x[0] > lim2) || 
    //             (x[1] < lim1) || (x[1] > lim2)){

    //         }else{
    //             if ((x[0] > lim1 + tick) && (x[0] < lim2 - tick) &&
    //                 (x[1] > lim1 + tick) && (x[1] < lim2 - tick)){

    //             }else{
    //                 flag = 1;
    //                 break;
    //             };
    //         };
    //     };
    //     if (flag > 0) {
    //         elementsGlueZoneCoarse_.push_back(jel);
    //         fMeshVector[0]->ElementVec()[jel] -> setGlueZone();

    //         Glue *el = new Glue(index++,jel);
    //         glueZoneCoarse_.push_back(el);
            
    //         for (int i=0; i < fMeshVector[0]->ElementVec()[jel] -> getNumberOfIntegrationPoints(); i++){
    //             VecDouble x = fMeshVector[0]->ElementVec()[jel] -> getIntegPointCoordinatesValue(i);

    //             if ((x[0] < lim1) || (x[0] > lim2) || 
    //                 (x[1] < lim1) || (x[1] > lim2)){

    //             }else{
    //                 if ((x[0] > lim1 + tick) && (x[0] < lim2 - tick) &&
    //                     (x[1] > lim1 + tick) && (x[1] < lim2 - tick)){

    //                 }else{
    //                 fMeshVector[0]->ElementVec()[jel] -> setIntegPointInGlueZone(i);
    //                 };
    //             };
    //         };
    //     };        
    // };

    // //Defines which nodes are in the glue zone
    // numElemGlueZoneCoarse = elementsGlueZoneCoarse_.size();
    // for (int i = 0; i < numElemGlueZoneCoarse; i++){
    //     VecInt connec = fMeshVector[0]->ElementVec()[elementsGlueZoneCoarse_[i]] -> getConnectivity();
    //     int nElNodes = fMeshVector[0]->NElNodes();
    //     for (int ino = 0; ino < nElNodes; ino++){
    //         nodesCZ2[connec[ino]] += 1;
    //     };
    // };

    // //Compute number of nodes in the glue zone
    // numNodesGlueZoneCoarse = 0;
    // for (int i = 0; i < fMeshVector[0]->NNodes(); i++){
    //     if(nodesCZ2[i] > 0) {
    //         numNodesGlueZoneCoarse += 1;
    //         nodesGlueZoneCoarse_.push_back(i);
    //     };
    // };

    // for (int i = 0; i < numNodesGlueZoneCoarse; i++){
    //     VecDouble x = fMeshVector[0]->NodeVec()[nodesGlueZoneCoarse_[i]] -> getCoordinates();
        
    //     Node *no = new Node(x,i);
    //     nodesLagrangeCoarse_.push_back(no);
    // };

    // for (int i = 0; i < numElemGlueZoneCoarse; i++){
    //     int nElNodes = fMeshVector[0]->NElNodes();
    //     VecInt connecAux(nElNodes);

    //     VecInt connec = fMeshVector[0]->ElementVec()[elementsGlueZoneCoarse_[i]] -> getConnectivity();
        
    //     for (int ino = 0; ino < numNodesGlueZoneCoarse; ino++)
    //         for (int k = 0; k < nElNodes; k++)
    //             if (nodesGlueZoneCoarse_[ino] == connec[k]) connecAux[k] = ino;
            
        
    //     glueZoneCoarse_[i] -> setConnectivity(connecAux);
    // };
     
};


//------------------------------------------------------------------------------
//---------------------------SETS THE WEIGHT FUNCTION---------------------------
//------------------------------------------------------------------------------
double Arlequin::GlobalWeightFunction(double dist){

    if (dist < 0){
        return 1.;
        //if (fabs(r) < 1.e-5) wFuncValue = 0.5;
    } else {
        if (dist >= fGlueZoneThickness*1.01){
            return fArlequinEpsilon;
        } else {
            // wFuncValue = 1. - (1. - epsilon) / lambda * r;
            double wFuncValue = 0.5;
            //wFuncValue = epsilon;
            // wFuncValue = 1+3.*(epsilon-1.)/(lambda*lambda) * r * r
            //     -2.*(epsilon-1.)/(lambda*lambda*lambda) * r * r * r;
            
            if (wFuncValue < fArlequinEpsilon) wFuncValue = fArlequinEpsilon;
            return wFuncValue;
        };
    };


}

double Arlequin::LocalWeightFunction(double dist){
    if (dist >= fGlueZoneThickness*1.01){
        return 1. - fArlequinEpsilon;
    } else {
        // wFuncValue = (1. - epsilon) / lambda * r;
        double wFuncValue = 0.5;
        //wFuncValue = 1. - epsilon;
        // wFuncValue = -3.*(epsilon-1.)/(lambda*lambda) * r * r
        //         +2.*(epsilon-1.)/(lambda*lambda*lambda) * r * r * r;

        if (wFuncValue > (1. - fArlequinEpsilon)) wFuncValue = 1. - fArlequinEpsilon;
        return wFuncValue;
    };  
}


void Arlequin::setWeightFunction(double val){
    
    double epsilon = fArlequinEpsilon;
    double lambda = fGlueZoneThickness*1.01;
 
    for (int i = 0; i < fMeshVector[0]->NNodes(); i++){
        
        double r = fGlobalSignaledDistance[i];
        double wFuncValue = GlobalWeightFunction(r);
        fMeshVector[0]->NodeVec()[i] -> setWeightFunction(wFuncValue);
    };         

    for (int jel = 0; jel < fMeshVector[0]->NElements(); jel++){
        for (int i = 0; i < fMeshVector[0]->ElementVec()[jel]->getNumberOfIntegrationPoints(); i++){
            double dist = fMeshVector[0]->ElementVec()[jel]->GetIntPointDistFunction(i);
            double wFuncValue = GlobalWeightFunction(dist);
            fMeshVector[0]->ElementVec()[jel] -> setIntegPointWeightFunction(i,wFuncValue); 
        }
    };


    for (int i=0; i<fMeshVector[1]->NNodes(); i++){
        double r = fLocalSignaledDistance[i];
        double wFuncValue = LocalWeightFunction(r);
        fMeshVector[1]->NodeVec()[i] -> setWeightFunction(wFuncValue);
    };
    
    for (int jel = 0; jel < fMeshVector[1]->NElements(); jel++){
        for (int i = 0; i < fMeshVector[1]->ElementVec()[jel]->getNumberOfIntegrationPoints(); i++){
            double dist = fMeshVector[1]->ElementVec()[jel]->GetIntPointDistFunction(i);
            double wFuncValue = LocalWeightFunction(dist);
            fMeshVector[1]->ElementVec()[jel] -> setIntegPointWeightFunction(i,wFuncValue);    
        }
            
    };

    return;
};


void Arlequin::DeleteCoarseEls(){

    std::set<int64_t> KeptNodes;
    std::set<int64_t> deletedElements;
    int64_t nElements = fMeshVector[0]->NElements();
    for (auto it=fMeshVector[0]->ElementVec().begin(); it!=fMeshVector[0]->ElementVec().end();){
        // auto elem = fMeshVector[0]->ElementVec()[iel];
        
        int counter = 0;
        for (int iintPoint = 0; iintPoint < (*it)->getNumberOfIntegrationPoints(); iintPoint++){
            if ((*it)->getIntegPointWeightFunction(iintPoint) > fArlequinEpsilon*1.1) continue;
            counter++;
        }

        if (counter == (*it)->getNumberOfIntegrationPoints()){
            //These are the elements to be deleted
            // delete fMeshVector[0]->ElementVec()[iel];
            deletedElements.insert((*it)->Index());
            fMeshVector[0]->ElementVec().erase(it);
        } else {
            //These are the elements to be kept, so store the nodes that will not be deleted.
            auto connec = (*it)->getConnectivity(); 
            for (int i = 0; i < connec.size(); i++) KeptNodes.insert(connec[i]);
            ++it;
        }
        
    }
    
    //Now delete the nodes;
    int64_t nNodes = fMeshVector[0]->NNodes();
    // for (int64_t inode = 0; inode < nNodes; inode++){
    int count = 0;
    std::set<int64_t> deletedNodes;
    for (auto it=fMeshVector[0]->NodeVec().begin(); it!=fMeshVector[0]->NodeVec().end();){
        // std::cout << "Node index " << count << " " << (*it)->Index() << std::endl;
        // count++;
        if (KeptNodes.find((*it)->Index()) == KeptNodes.end()) {
            deletedNodes.insert((*it)->Index());
            fMeshVector[0]->NodeVec().erase(it);
        } else {
            ++it;
        }
    }   

    fMeshVector[0]->ElementVec().shrink_to_fit();
    fMeshVector[0]->NodeVec().shrink_to_fit();
    //Finally, update the connects
    std::map<int64_t,int64_t> prevToNewNodeIndex;
    
    nNodes = fMeshVector[0]->NNodes();
    for (int64_t inode = 0; inode < nNodes; inode++){
        prevToNewNodeIndex[fMeshVector[0]->NodeVec()[inode]->Index()]=inode;
        fMeshVector[0]->NodeVec()[inode]->Index()=inode;
    }
    
    nElements = fMeshVector[0]->NElements();
    for (int64_t iel = 0; iel < nElements; iel++){
        auto connec = fMeshVector[0]->ElementVec()[iel]->getConnectivity();
        for (int i = 0; i < connec.size(); i++) connec[i]=prevToNewNodeIndex[connec[i]];
        fMeshVector[0]->ElementVec()[iel]->setConnectivity(connec);
    }
    
}

//------------------------------------------------------------------------------
//-----------------------------PRINT COARSE RESULTS-----------------------------
//------------------------------------------------------------------------------

void Arlequin::printResultsCoarse(int step) {

    // // PRINT COARSE MODEL RESULTS
    // std::string result;
    // std::ostringstream convert;

    // convert << step+100000;
    // result = convert.str();

    // std::string s = "saidaVelCoarse"+result+".xdmf";
    // std::fstream output_v(s.c_str(), std::ios_base::out);

    // std::string s1 = "resultCoarse"+result+".h5";
    // std::fstream filename(s1.c_str(), std::ios_base::out);

    // std::string s2 = "geometryCoarse.h5";
    // if (step == 0){
    //     std::fstream filename(s2.c_str(), std::ios_base::out);
    // }

    // //Auxiliary vectors to write HDF5 file
    // double *pointVector;
    // int *connec2;
    // double *pointScalar;
    // int *intCellScalar;
    // double *douCellScalar;
    // int nElNodes = fMeshVector[0]->NElNodes();

    // pointVector = new double[3*fMeshVector[0]->NNodes()]();
    // connec2 = new int[nElNodes*fMeshVector[0]->NElements()]();
    // pointScalar = new double[fMeshVector[0]->NNodes()]();
    // intCellScalar = new int[fMeshVector[0]->NElements()]();
    // douCellScalar = new double[fMeshVector[0]->NElements()]();

    // for (int i = 0; i < fMeshVector[0]->NNodes(); i++) {
    //     VecDouble x = fMeshVector[0]->NodeVec()[i] -> getCoordinates();
    //     pointVector[3*i  ] = x[0];
    //     pointVector[3*i+1] = x[1];
    //     pointVector[3*i+2] = 0.0;
    // }
    // for (int iElem = 0; iElem < fMeshVector[0]->NElements(); iElem++){
    //     VecInt con = fMeshVector[0]->ElementVec()[iElem] -> getConnectivity();

    //     if (DIM == 2){
    //         switch (DEG)
    //         {
    //         case 1:
    //             for (int i = 0; i < nElNodes; ++i) connec2[nElNodes*iElem+i] = con[i];
    //             break;
    //         case 2:
    //             {
    //                 connec2[nElNodes*iElem+0] = con[4];
    //                 connec2[nElNodes*iElem+1] = con[1];
    //                 connec2[nElNodes*iElem+2] = con[3];
    //                 connec2[nElNodes*iElem+3] = con[0];
    //                 connec2[nElNodes*iElem+4] = con[5];
    //                 connec2[nElNodes*iElem+5] = con[2];
    //                 break;
    //             }
    //         case 3:
    //             {
    //                 connec2[nElNodes*iElem+0] = con[6];
    //                 connec2[nElNodes*iElem+1] = con[9];
    //                 connec2[nElNodes*iElem+2] = con[5];
    //                 connec2[nElNodes*iElem+3] = con[1];
    //                 connec2[nElNodes*iElem+4] = con[4];
    //                 connec2[nElNodes*iElem+5] = con[3];
    //                 connec2[nElNodes*iElem+6] = con[0];
    //                 connec2[nElNodes*iElem+7] = con[8];
    //                 connec2[nElNodes*iElem+8] = con[7];
    //                 connec2[nElNodes*iElem+9] = con[2];
    //                 break;
    //             }
            
    //         default:
    //             PanicButton();
    //             break;
    //         }
    //     } else {
    //         for (int i = 0; i < nElNodes; ++i) connec2[nElNodes*iElem+i] = con[i];
    //     }
        
        
    // }
   
    // hid_t file, file2; 
    // hid_t dataset; 
    // hid_t dataspace;

    // herr_t status;
    // hsize_t xdim = fMeshVector[0]->NNodes();
    // hsize_t eldim = fMeshVector[0]->NElements();
    // hsize_t elnod = nElNodes;
    // hsize_t pointVectorDims[2] = { xdim, 3 };
    // hsize_t pointScalarDims[2] = { xdim, 1 };
    // hsize_t pointCellScalarDims[2] = { eldim, 1 };
    // hsize_t connec2Dims[2] = { eldim, elnod };

    // std::string elType;
    // if (DIM == 2) elType = "TRIANGLE";
    // if (DIM == 3) elType = "TETRAHEDRON";

    // //Create HDF5 file
    // file = H5Fcreate(s1.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    // if (step == 0) file2 = H5Fcreate(s2.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // //Write xdmf file
    // output_v << "<?xml version=\"1.0\"?>" << std::endl
    //          << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl
    //          << "<Xdmf Version=\"2.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\" >" << std::endl
    //          << "<Domain>" << std::endl
    //          << "  <Grid>"  << std::endl
    //          << "    <Topology TopologyType=\"" << elType << "\" NumberOfElements=\"" << fMeshVector[0]->NElements() << "\" >" << std::endl
    //          << "      <DataItem Format=    \"HDF\" NumberType=\"int\" Dimensions=\"" << fMeshVector[0]->NElements() << " " << nElNodes << "\" >" << std::endl;
    
    // //Connectivity
    // output_v << "        " << s2 << ":/connec" << std::endl;
    // //Start connectivity
    // if (step == 0){
    //     dataspace = H5Screate_simple(2, connec2Dims, NULL);
    //     dataset = H5Dcreate2(file2, "/connec", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &connec2[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    // }
    // //End connectivity

    // output_v << "      </DataItem>" << std::endl
    //          << "    </Topology>" << std::endl
    //          << "    <Geometry GeometryType=\"XYZ\">" << std::endl
    //          << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[0]->NNodes() <<" 3\">" << std::endl;

    // //Coordinates   
    // output_v << "        " << s2 << ":/coords" << std::endl;
    // //Start coordinates
    // if (step == 0){
    //     dataspace = H5Screate_simple(2, pointVectorDims, NULL);
    //     dataset = H5Dcreate2(file2, "/coords", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    // }
    // //End coordinates
    // output_v << "      </DataItem>" << std::endl
    //          << "    </Geometry>" << std::endl;

    // //LISTS
    // //Velocity
    // if (fMeshVector[0]->printVelocity){
    //     output_v << "    <Attribute Name=\"Velocity\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[0]->NNodes() <<" 3\">" << std::endl;  
    //     output_v << "        " << s1 << ":/velocity" << std::endl;     
    //     //Start Velocity
    //     for (int i = 0; i < fMeshVector[0]->NNodes(); i++) {
    //         pointVector[3*i  ] = fMeshVector[0]->NodeVec()[i] -> getVelocity(0);
    //         pointVector[3*i+1] = fMeshVector[0]->NodeVec()[i] -> getVelocity(1);
    //         pointVector[3*i+2] = 0.0;
    //     };
    //     dataspace = H5Screate_simple(2, pointVectorDims, NULL);
    //     dataset = H5Dcreate2(file, "/velocity", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Velocity
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //Acceleration
    // if (fMeshVector[0]->printAcceleration){
    //     output_v << "    <Attribute Name=\"Acceleration\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[0]->NNodes() <<" 3\">" << std::endl;  
    //     output_v << "        " << s1 << ":/acceleration" << std::endl;     
    //     //Start Acceleration
    //     for (int i = 0; i < fMeshVector[0]->NNodes(); i++) {
    //         pointVector[3*i  ] = fMeshVector[0]->NodeVec()[i] -> getAcceleration(0);
    //         pointVector[3*i+1] = fMeshVector[0]->NodeVec()[i] -> getAcceleration(1);
    //         pointVector[3*i+2] = 0.0;
    //     };
    //     dataspace = H5Screate_simple(2, pointVectorDims, NULL);
    //     dataset = H5Dcreate2(file, "/acceleration", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Acceleration
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //Real Velocity
    // if (fMeshVector[0]->printRealVelocity){
    //     output_v << "    <Attribute Name=\"Real Velocity\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[0]->NNodes() <<" 3\">" << std::endl;  
    //     output_v << "        " << s1 << ":/realVelocity" << std::endl;     
    //     //Start Real Velocity
    //     for (int i = 0; i < fMeshVector[0]->NNodes(); i++) {
    //         pointVector[3*i  ] = fMeshVector[0]->NodeVec()[i] -> getVelocityArlequin(0);
    //         pointVector[3*i+1] = fMeshVector[0]->NodeVec()[i] -> getVelocityArlequin(1);
    //         pointVector[3*i+2] = 0.0;
    //     };
    //     dataspace = H5Screate_simple(2, pointVectorDims, NULL);
    //     dataset = H5Dcreate2(file, "/realVelocity", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Real Velocity
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //Lagrange Multipliers
    // if (fMeshVector[0]->printLagrangeMultipliers){
    //     output_v << "    <Attribute Name=\"Lagrange Multipliers\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[0]->NNodes() <<" 3\">" << std::endl;  
    //     output_v << "        " << s1 << ":/lagrangeMultiplers" << std::endl;     
    //     //Start Lagrange Multipliers
    //     for (int i = 0; i < fMeshVector[0]->NNodes(); i++) {
    //         pointVector[3*i  ] = fMeshVector[0]->NodeVec()[i] -> getLagrangeMultiplier(0);
    //         pointVector[3*i+1] = fMeshVector[0]->NodeVec()[i] -> getLagrangeMultiplier(1);
    //         pointVector[3*i+2] = 0.0;
    //     };
    //     dataspace = H5Screate_simple(2, pointVectorDims, NULL);
    //     dataset = H5Dcreate2(file, "/lagrangeMultiplers", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Real Velocity
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //Dist Function
    // if (fMeshVector[0]->printDistFunction){
    //     output_v << "    <Attribute Name=\"Dist Function\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\"1 "<< fMeshVector[0]->NNodes() <<"\">" << std::endl;  
    //     output_v << "        " << s1 << ":/distfunction" << std::endl;     
    //     //Start Dist Function
    //     for (int i = 0; i < fMeshVector[0]->NNodes(); i++) pointScalar[i] = fMeshVector[0]->NodeVec()[i] -> getDistFunction();
    //     dataspace = H5Screate_simple(2, pointScalarDims, NULL);
    //     dataset = H5Dcreate2(file, "/distfunction", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Dist Function
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //Weight Function
    // if (fMeshVector[0]->printEnergyWeightFunction){
    //     output_v << "    <Attribute Name=\"Weight Function\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\"1 "<< fMeshVector[0]->NNodes() <<"\">" << std::endl;  
    //     output_v << "        " << s1 << ":/weightFunction" << std::endl;     
    //     //Start Weigth Function
    //     for (int i = 0; i < fMeshVector[0]->NNodes(); i++) pointScalar[i] = fMeshVector[0]->NodeVec()[i] -> getWeightFunction();
    //     dataspace = H5Screate_simple(2, pointScalarDims, NULL);
    //     dataset = H5Dcreate2(file, "/weightFunction", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Weight Function
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //Pressure
    // if (fMeshVector[0]->printPressure){
    //     output_v << "    <Attribute Name=\"Pressure\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[0]->NNodes() <<" 1\">" << std::endl;  
    //     output_v << "        " << s1 << ":/pressure" << std::endl;     
    //     //Start Pressure
    //     for (int i = 0; i < fMeshVector[0]->NNodes(); i++) pointScalar[i] = fMeshVector[0]->NodeVec()[i] -> getPressure();
    //     dataspace = H5Screate_simple(2, pointScalarDims, NULL);
    //     dataset = H5Dcreate2(file, "/pressure", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Pressure
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //Real Pressure
    // if (fMeshVector[0]->printRealPressure){
    //     output_v << "    <Attribute Name=\"Real Pressure\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[0]->NNodes() <<" 1\">" << std::endl;  
    //     output_v << "        " << s1 << ":/realPressure" << std::endl;     
    //     //Start Real Pressure
    //     for (int i = 0; i < fMeshVector[0]->NNodes(); i++) pointScalar[i] = fMeshVector[0]->NodeVec()[i] -> getPressureArlequin();
    //     dataspace = H5Screate_simple(2, pointScalarDims, NULL);
    //     dataset = H5Dcreate2(file, "/realPressure", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Real Pressure
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //Process
    // if (fMeshVector[0]->printProcess){
    //     output_v << "    <Attribute Name=\"Process\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"int\" Dimensions=\""<< fMeshVector[0]->NElements() <<" 1\">" << std::endl;  
    //     output_v << "        " << s1 << ":/process" << std::endl;
    //     //Start Process
    //     dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
    //     dataset = H5Dcreate2(file, "/process", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &domDecompCoarse.first[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Process
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //WeightFunction2
    // if (fMeshVector[0]->printEnergyWeightFunction){
    //     output_v << "    <Attribute Name=\"Weight Function\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[0]->NElements() <<" 1\">" << std::endl;  
    //     output_v << "        " << s1 << ":/weightFunction2" << std::endl;     
    //     //Start Weight Function
    //     for (int i = 0; i < fMeshVector[0]->NElements(); i++) douCellScalar[i] = fMeshVector[0]->ElementVec()[i] -> getIntegPointWeightFunction(0);
    //     dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
    //     dataset = H5Dcreate2(file, "/weightFunction2", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &douCellScalar[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Weight Function
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //Lines
    // if (fMeshVector[0]->printLines){
    //     output_v << "    <Attribute Name=\"Lines\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"int\" Dimensions=\""<< fMeshVector[0]->NElements() <<" 1\">" << std::endl;  
    //     output_v << "        " << s1 << ":/lines" << std::endl;     
    //     //Start Lines
    //     for (int i = 0; i < fMeshVector[0]->NElements(); i++){
    //         int res = 0;
    //         for (int j = 0; j < numBoundElemCoarse; j++)
    //            if (boundaryCoarse_[j] -> getElement() == i) res = boundaryCoarse_[j] -> getBoundaryGroup();
    //         intCellScalar[i] = res;
    //     };
    //     dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
    //     dataset = H5Dcreate2(file, "/lines", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &intCellScalar[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Lines
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // // Jacobian
    // if (fMeshVector[0]->printJacobian){
    //     output_v << "    <Attribute Name=\"Jacobian\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[0]->NElements() <<" 1\">" << std::endl;  
    //     output_v << "        " << s1 << ":/jacobian" << std::endl;     
    //     //Start Weight Function
    //     for (int i = 0; i < fMeshVector[0]->NElements(); i++) douCellScalar[i] = fMeshVector[0]->ElementVec()[i] -> getJacobian();
    //     dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
    //     dataset = H5Dcreate2(file, "/jacobian", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &douCellScalar[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Weight Function
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // // END FILE
    // output_v << "  </Grid>" << std::endl
    //          << "</Domain>" << std::endl
    //          << "</Xdmf>" << std::endl;

    // //Delete auxiliary vectors
    // delete [] pointVector;
    // delete [] connec2;
    // delete [] pointScalar;
    // delete [] intCellScalar;
    // delete [] douCellScalar;

    // //End HDF5 file
    // status = H5Fclose(file);
    // if (step == 0) status = H5Fclose(file2);

    return;
};

//------------------------------------------------------------------------------
//------------------------------PRINT FINE RESULTS------------------------------
//------------------------------------------------------------------------------

void Arlequin::printResultsFine(int step) {

    // //PRINT FINE MODEL RESULTS
    // std::string result;
    // std::ostringstream convert;

    // convert << step+100000;
    // result = convert.str();

    // std::string s = "saidaVelFine"+result+".xdmf";
    // std::fstream output_v(s.c_str(), std::ios_base::out);

    // std::string s1 = "resultFine"+result+".h5";
    // std::fstream filename(s1.c_str(), std::ios_base::out);

    // //Auxiliary vectors to write HDF5 file
    // double *pointVector;
    // int *connec2;
    // double *pointScalar;
    // int *intCellScalar;
    // double *douCellScalar;
    // int nElNodes = fMeshVector[1]->NElNodes();

    // pointVector = new double[3*fMeshVector[1]->NNodes()]();
    // connec2 = new int[nElNodes*fMeshVector[1]->NElements()]();
    // pointScalar = new double[fMeshVector[1]->NNodes()]();
    // intCellScalar = new int[fMeshVector[1]->NElements()]();
    // douCellScalar = new double[fMeshVector[1]->NElements()]();

    // for (int i = 0; i < fMeshVector[1]->NNodes(); i++) {
    //     VecDouble x = fMeshVector[1]->NodeVec()[i] -> getCoordinates();
    //     pointVector[3*i  ] = x[0];
    //     pointVector[3*i+1] = x[1];
    //     pointVector[3*i+2] = 0.0;
    // }

    // for (int iElem = 0; iElem < fMeshVector[1]->NElements(); iElem++){
    //     VecInt con = fMeshVector[1]->ElementVec()[iElem] -> getConnectivity();
    //     if (DIM == 2){
    //         switch (DEG)
    //         {
    //         case 1:
    //             for (int i = 0; i < nElNodes; ++i) connec2[nElNodes*iElem+i] = con[i];
    //             break;
    //         case 2:
    //             {
    //                 connec2[nElNodes*iElem+0] = con[4];
    //                 connec2[nElNodes*iElem+1] = con[1];
    //                 connec2[nElNodes*iElem+2] = con[3];
    //                 connec2[nElNodes*iElem+3] = con[0];
    //                 connec2[nElNodes*iElem+4] = con[5];
    //                 connec2[nElNodes*iElem+5] = con[2];
    //                 break;
    //             }
    //         case 3:
    //             {
    //                 connec2[nElNodes*iElem+0] = con[6];
    //                 connec2[nElNodes*iElem+1] = con[9];
    //                 connec2[nElNodes*iElem+2] = con[5];
    //                 connec2[nElNodes*iElem+3] = con[1];
    //                 connec2[nElNodes*iElem+4] = con[4];
    //                 connec2[nElNodes*iElem+5] = con[3];
    //                 connec2[nElNodes*iElem+6] = con[0];
    //                 connec2[nElNodes*iElem+7] = con[8];
    //                 connec2[nElNodes*iElem+8] = con[7];
    //                 connec2[nElNodes*iElem+9] = con[2];
    //                 break;
    //             }      
            
    //         default:
    //             PanicButton();
    //             break;
    //         }
    //     } else {
    //         for (int i = 0; i < nElNodes; ++i) connec2[nElNodes*iElem+i] = con[i];
    //     }
    // } 

    // hid_t file; 
    // hid_t dataset; 
    // hid_t dataspace;

    // herr_t status;
    // hsize_t xdim = fMeshVector[1]->NNodes();
    // hsize_t eldim = fMeshVector[1]->NElements();
    // hsize_t elnod = nElNodes;
    // hsize_t pointVectorDims[2] = { xdim, 3 };
    // hsize_t pointScalarDims[2] = { xdim, 1 };
    // hsize_t pointCellScalarDims[2] = { eldim, 1 };
    // hsize_t connec2Dims[2] = { eldim, elnod };

    // std::string elType;
    // if (DIM == 2) elType = "TRIANGLE";
    // if (DIM == 3) elType = "TETRAHEDRON";

    // //Create HDF5 file
    // file = H5Fcreate(s1.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // //Write xdmf file
    // output_v << "<?xml version=\"1.0\"?>" << std::endl
    //          << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl
    //          << "<Xdmf Version=\"2.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\" >" << std::endl
    //          << "<Domain>" << std::endl
    //          << "  <Grid>"  << std::endl
    //          << "    <Topology TopologyType=\"" << elType << "\" NumberOfElements=\"" << fMeshVector[1]->NElements() << "\" >" << std::endl
    //          << "      <DataItem Format=    \"HDF\" NumberType=\"int\" Dimensions=\"" << fMeshVector[1]->NElements() << " " << nElNodes << "\" >" << std::endl;
    
    // //Connectivity
    // output_v << "        " << s1 << ":/connec" << std::endl;
    // //Start connectivity
    // dataspace = H5Screate_simple(2, connec2Dims, NULL);
    // dataset = H5Dcreate2(file, "/connec", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &connec2[0]);
    // status = H5Dclose(dataset);
    // status = H5Sclose(dataspace);
    // //End connectivity

    // output_v << "      </DataItem>" << std::endl
    //          << "    </Topology>" << std::endl
    //          << "    <Geometry GeometryType=\"XYZ\">" << std::endl
    //          << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[1]->NNodes() <<" 3\">" << std::endl;

    // //Coordinates   
    // output_v << "        " << s1 << ":/coords" << std::endl;
    // //Start coordinates
    // dataspace = H5Screate_simple(2, pointVectorDims, NULL);
    // dataset = H5Dcreate2(file, "/coords", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
    // status = H5Dclose(dataset);
    // status = H5Sclose(dataspace);
    // //End coordinates
    // output_v << "      </DataItem>" << std::endl
    //          << "    </Geometry>" << std::endl;

    // //LISTS
    // //Velocity
    // if (fMeshVector[1]->printVelocity){
    //     output_v << "    <Attribute Name=\"Velocity\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[1]->NNodes() <<" 3\">" << std::endl;  
    //     output_v << "        " << s1 << ":/velocity" << std::endl;     
    //     //Start Velocity
    //     for (int i = 0; i < fMeshVector[1]->NNodes(); i++) {
    //         pointVector[3*i  ] = fMeshVector[1]->NodeVec()[i] -> getVelocity(0);
    //         pointVector[3*i+1] = fMeshVector[1]->NodeVec()[i] -> getVelocity(1);
    //         pointVector[3*i+2] = 0.0;
    //     };
    //     dataspace = H5Screate_simple(2, pointVectorDims, NULL);
    //     dataset = H5Dcreate2(file, "/velocity", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Velocity
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //Acceleration
    // if (fMeshVector[1]->printAcceleration){
    //     output_v << "    <Attribute Name=\"Acceleration\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[1]->NNodes() <<" 3\">" << std::endl;  
    //     output_v << "        " << s1 << ":/acceleration" << std::endl;     
    //     //Start Acceleration
    //     for (int i = 0; i < fMeshVector[1]->NNodes(); i++) {
    //         pointVector[3*i  ] = fMeshVector[1]->NodeVec()[i] -> getAcceleration(0);
    //         pointVector[3*i+1] = fMeshVector[1]->NodeVec()[i] -> getAcceleration(1);
    //         pointVector[3*i+2] = 0.0;
    //     };
    //     dataspace = H5Screate_simple(2, pointVectorDims, NULL);
    //     dataset = H5Dcreate2(file, "/acceleration", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Acceleration
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //Real Velocity
    // if (fMeshVector[1]->printRealVelocity){
    //     output_v << "    <Attribute Name=\"Real Velocity\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[1]->NNodes() <<" 3\">" << std::endl;  
    //     output_v << "        " << s1 << ":/realVelocity" << std::endl;     
    //     //Start Real Velocity
    //     for (int i = 0; i < fMeshVector[1]->NNodes(); i++) {
    //         pointVector[3*i  ] = fMeshVector[1]->NodeVec()[i] -> getVelocityArlequin(0);
    //         pointVector[3*i+1] = fMeshVector[1]->NodeVec()[i] -> getVelocityArlequin(1);
    //         pointVector[3*i+2] = 0.0;
    //     };
    //     dataspace = H5Screate_simple(2, pointVectorDims, NULL);
    //     dataset = H5Dcreate2(file, "/realVelocity", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Real Velocity
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //Lagrange Multipliers
    // if (fMeshVector[1]->printLagrangeMultipliers){
    //     output_v << "    <Attribute Name=\"Lagrange Multipliers\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[1]->NNodes() <<" 3\">" << std::endl;  
    //     output_v << "        " << s1 << ":/lagrangeMultiplers" << std::endl;     
    //     //Start Lagrange Multipliers
    //     for (int i = 0; i < fMeshVector[1]->NNodes(); i++) {
    //         pointVector[3*i  ] = fMeshVector[1]->NodeVec()[i] -> getLagrangeMultiplier(0);
    //         pointVector[3*i+1] = fMeshVector[1]->NodeVec()[i] -> getLagrangeMultiplier(1);
    //         pointVector[3*i+2] = 0.0;
    //     };
    //     dataspace = H5Screate_simple(2, pointVectorDims, NULL);
    //     dataset = H5Dcreate2(file, "/lagrangeMultiplers", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Lagrange Multipliers
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //Normal vector
    // if (fMeshVector[1]->printInnerNormal){
    //     output_v << "    <Attribute Name=\"Normal Vector\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[1]->NNodes() <<" 3\">" << std::endl;  
    //     output_v << "        " << s1 << ":/normalVector" << std::endl;     
    //     //Start Normal vector
    //     for (int i = 0; i < fMeshVector[1]->NNodes(); i++) {
    //         VecDouble n = fMeshVector[1]->NodeVec()[i] -> getInnerNormal();
    //         pointVector[3*i  ] = n[0];
    //         pointVector[3*i+1] = n[1];
    //         pointVector[3*i+2] = 0.0;
    //     };
    //     dataspace = H5Screate_simple(2, pointVectorDims, NULL);
    //     dataset = H5Dcreate2(file, "/normalVector", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Normal Vector
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //Mesh Velocity
    // if (fMeshVector[1]->printMeshVelocity){
    //     output_v << "    <Attribute Name=\"Mesh Velocity\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[1]->NNodes() <<" 3\">" << std::endl;  
    //     output_v << "        " << s1 << ":/meshVelocity" << std::endl;     
    //     //Start mesh Velocity
    //     for (int i = 0; i < fMeshVector[1]->NNodes(); i++) {
    //         pointVector[3*i  ] = fMeshVector[1]->NodeVec()[i] -> getMeshVelocity(0);
    //         pointVector[3*i+1] = fMeshVector[1]->NodeVec()[i] -> getMeshVelocity(1);
    //         pointVector[3*i+2] = 0.0;
    //     };
    //     dataspace = H5Screate_simple(2, pointVectorDims, NULL);
    //     dataset = H5Dcreate2(file, "/meshVelocity", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Mesh Velocity
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // };

    // //Mesh Displacement
    // if (fMeshVector[1]->printMeshDisplacement){
    //     output_v << "    <Attribute Name=\"Mesh Displacement\" Center=\"Node\" AttributeType=\"Vector\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[1]->NNodes() <<" 3\">" << std::endl;  
    //     output_v << "        " << s1 << ":/meshDisplacement" << std::endl;     
    //     //Start Mesh displacement
    //     for (int i = 0; i < fMeshVector[1]->NNodes(); i++) {
    //         VecDouble x = fMeshVector[1]->NodeVec()[i] -> getCoordinates();
    //         VecDouble xi = fMeshVector[1]->NodeVec()[i] -> getInitialCoordinates();
    //         pointVector[3*i  ] = x[0] - xi[0];
    //         pointVector[3*i+1] = x[1] - xi[1];
    //         pointVector[3*i+2] = 0.0;
    //     };
    //     dataspace = H5Screate_simple(2, pointVectorDims, NULL);
    //     dataset = H5Dcreate2(file, "/meshDisplacement", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointVector[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Mesh displacement
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // };

    // //Element correspondence
    // if (fMeshVector[1]->printElementCorrespondence){
    //     output_v << "    <Attribute Name=\"Element\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\"1 "<< fMeshVector[1]->NNodes() <<"\">" << std::endl;  
    //     output_v << "        " << s1 << ":/elementCorresp" << std::endl;     
    //     //Start Element correspondence
    //     for (int i = 0; i < fMeshVector[1]->NNodes(); i++) pointScalar[i] = fMeshVector[1]->NodeVec()[i] -> getNodalElemCorrespondence();
    //     dataspace = H5Screate_simple(2, pointScalarDims, NULL);
    //     dataset = H5Dcreate2(file, "/elementCorresp", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Element correspondence
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //Dist Function
    // if (fMeshVector[1]->printDistFunction){
    //     output_v << "    <Attribute Name=\"Dist Function\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\"1 "<< fMeshVector[1]->NNodes() <<"\">" << std::endl;  
    //     output_v << "        " << s1 << ":/distfunction" << std::endl;     
    //     //Start Dist Function
    //     for (int i = 0; i < fMeshVector[1]->NNodes(); i++) pointScalar[i] = fMeshVector[1]->NodeVec()[i] -> getDistFunction();
    //     dataspace = H5Screate_simple(2, pointScalarDims, NULL);
    //     dataset = H5Dcreate2(file, "/distfunction", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Dist Function
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //Weight Function
    // if (fMeshVector[1]->printEnergyWeightFunction){
    //     output_v << "    <Attribute Name=\"Weight Function\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\"1 "<< fMeshVector[1]->NNodes() <<"\">" << std::endl;  
    //     output_v << "        " << s1 << ":/weightFunction" << std::endl;     
    //     //Start Weigth Function
    //     for (int i = 0; i < fMeshVector[1]->NNodes(); i++) pointScalar[i] = fMeshVector[1]->NodeVec()[i] -> getWeightFunction();
    //     dataspace = H5Screate_simple(2, pointScalarDims, NULL);
    //     dataset = H5Dcreate2(file, "/weightFunction", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Weight Function
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //Pressure
    // if (fMeshVector[1]->printPressure){
    //     output_v << "    <Attribute Name=\"Pressure\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[1]->NNodes() <<" 1\">" << std::endl;  
    //     output_v << "        " << s1 << ":/pressure" << std::endl;     
    //     //Start Pressure
    //     for (int i = 0; i < fMeshVector[1]->NNodes(); i++) pointScalar[i] = fMeshVector[1]->NodeVec()[i] -> getPressure();
    //     dataspace = H5Screate_simple(2, pointScalarDims, NULL);
    //     dataset = H5Dcreate2(file, "/pressure", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Pressure
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //Real Pressure
    // if (fMeshVector[1]->printRealPressure){
    //     output_v << "    <Attribute Name=\"Real Pressure\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[1]->NNodes() <<" 1\">" << std::endl;  
    //     output_v << "        " << s1 << ":/realPressure" << std::endl;     
    //     //Start Real Pressure
    //     for (int i = 0; i < fMeshVector[1]->NNodes(); i++) pointScalar[i] = fMeshVector[1]->NodeVec()[i] -> getPressureArlequin();
    //     dataspace = H5Screate_simple(2, pointScalarDims, NULL);
    //     dataset = H5Dcreate2(file, "/realPressure", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Real Pressure
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // if (fMeshVector[1]->printVorticity){
    //     output_v << "    <Attribute Name=\"Vorticity\" Center=\"Node\" AttributeType=\"Scalar\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[1]->NNodes() <<" 1\">" << std::endl;  
    //     output_v << "        " << s1 << ":/vorticity" << std::endl;     
    //     //Start Real Pressure
    //     for (int i = 0; i < fMeshVector[1]->NNodes(); i++) pointScalar[i] = fMeshVector[1]->NodeVec()[i] -> getVorticity();
    //     dataspace = H5Screate_simple(2, pointScalarDims, NULL);
    //     dataset = H5Dcreate2(file, "/vorticity", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pointScalar[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Real Pressure
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // };

    // //Process
    // if (fMeshVector[1]->printProcess){
    //     output_v << "    <Attribute Name=\"Process\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"int\" Dimensions=\""<< fMeshVector[1]->NElements() <<" 1\">" << std::endl;  
    //     output_v << "        " << s1 << ":/process" << std::endl;
    //     //Start Process
    //     dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
    //     dataset = H5Dcreate2(file, "/process", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &domDecompFine.first[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Process
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //WeightFunction2
    // if (fMeshVector[1]->printEnergyWeightFunction){
    //     output_v << "    <Attribute Name=\"Weight Function\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[1]->NElements() <<" 1\">" << std::endl;  
    //     output_v << "        " << s1 << ":/weightFunction2" << std::endl;     
    //     //Start Weight Function
    //     for (int i = 0; i < fMeshVector[1]->NElements(); i++) douCellScalar[i] = fMeshVector[1]->ElementVec()[i] -> getIntegPointWeightFunction(0);
    //     dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
    //     dataset = H5Dcreate2(file, "/weightFunction2", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &douCellScalar[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Weight Function
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //Lines
    // if (fMeshVector[1]->printLines){
    //     output_v << "    <Attribute Name=\"Lines\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"int\" Dimensions=\""<< fMeshVector[1]->NElements() <<" 1\">" << std::endl;  
    //     output_v << "        " << s1 << ":/lines" << std::endl;     
    //     //Start Lines
    //     for (int i = 0; i < fMeshVector[1]->NElements(); i++){
    //         int res = 0;
    //         for (int j = 0; j < fMeshVector[1]->NBoundElements(); j++)
    //            if (boundaryFine_[j] -> getElement() == i) res = boundaryFine_[j] -> getBoundaryGroup();
    //         intCellScalar[i] = res;
    //     };
    //     dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
    //     dataset = H5Dcreate2(file, "/lines", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &intCellScalar[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Lines
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // // Jacobian
    // if (fMeshVector[1]->printJacobian){
    //     output_v << "    <Attribute Name=\"Jacobian\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\""<< fMeshVector[1]->NElements() <<" 1\">" << std::endl;  
    //     output_v << "        " << s1 << ":/jacobian" << std::endl;     
    //     //Start Weight Function
    //     for (int i = 0; i < fMeshVector[1]->NElements(); i++) douCellScalar[i] = fMeshVector[1]->ElementVec()[i] -> getJacobian();
    //     dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
    //     dataset = H5Dcreate2(file, "/jacobian", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &douCellScalar[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Weight Function
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // //Glue Zone
    // if (fMeshVector[1]->printGlueZone){
    //     output_v << "    <Attribute Name=\"Glue Zone\" Center=\"Cell\" AttributeType=\"Scalar\" >" << std::endl
    //              << "      <DataItem Format=\"HDF\" NumberType=\"int\" Dimensions=\""<< fMeshVector[1]->NElements() <<" 1\">" << std::endl;  
    //     output_v << "        " << s1 << ":/glueZone" << std::endl;     
    //     //Start Lines
    //     for (int i = 0; i < fMeshVector[1]->NElements(); i++) intCellScalar[i] = fMeshVector[1]->ElementVec()[i] -> getGlueZoneInt();
    //     dataspace = H5Screate_simple(2, pointCellScalarDims, NULL);
    //     dataset = H5Dcreate2(file, "/glueZone", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &intCellScalar[0]);
    //     status = H5Dclose(dataset);
    //     status = H5Sclose(dataspace);
    //     //End Lines
    //     output_v << "      </DataItem>" << std::endl
    //              << "    </Attribute>" << std::endl;
    // }

    // // END FILE
    // output_v << "  </Grid>" << std::endl
    //          << "</Domain>" << std::endl
    //          << "</Xdmf>" << std::endl;

    // //Delete auxiliary vectors
    // delete [] pointVector;
    // delete [] connec2;
    // delete [] pointScalar;
    // delete [] intCellScalar;
    // delete [] douCellScalar;

    // //End HDF5 file
    // status = H5Fclose(file);

    return;
};


//------------------------------------------------------------------------------
//------------SETS COARSE/FINE MESHES AND GETS ITS BASIC INFORMATIONS-----------
//------------------------------------------------------------------------------
// 
// void Arlequin::setFluidModels(FluidMesh& coarse, FluidMesh& fine){

//     fMeshVector[0] = &coarse;
//     fMeshVector[1] = &fine;

//     //Gets Fine and Coarse models basic information
//     // fMeshVector[0]->NElements() = fMeshVector[0]->NElements();
//     // fMeshVector[1]->NElements()   = fMeshVector[1]->NElements();
//     // fMeshVector[0]->NNodes() = fMeshVector[0]->NNodes();
//     // fMeshVector[1]->NNodes()   = fMeshVector[1]->NNodes();
 
//     nodesCoarse_  = &fMeshVector[0]->NodeVec();
//     nodesFine_    = &fMeshVector[1]->NodeVec();

//     // fMeshVector[0]->ElementVec() = fMeshVector[0]->ElementVec();
//     // fMeshVector[1]->ElementVec()   = fMeshVector[1]->ElementVec();
 
//     boundaryCoarse_ = fMeshVector[0]->BoundaryVec();
//     boundaryFine_   = fMeshVector[1]->BoundaryVec();
//     fProbType = fMeshVector[1]->ProbType();

//     // fMeshVector[1]->NBoundElements() = boundaryFine_.size();
//     numBoundElemCoarse = boundaryCoarse_.size();

//     // domDecompCoarse = fMeshVector[0]->getDomainDecomposition();
//     // domDecompFine = fMeshVector[1]->getDomainDecomposition();
    
//     // numTimeSteps = fMeshVector[1]->getNumberOfTimeSteps();
//     dTime = fMeshVector[1]->getProblemParameters().GetTimeStep();

//     for (int i=0; i < fMeshVector[1]->NElements(); i++){
//         fMeshVector[1]->ElementVec()[i] -> setModel(true);
//     };
//     for (int i=0; i < fMeshVector[0]->NElements(); i++){
//         fMeshVector[0]->ElementVec()[i] -> setModel(false);
//     };

//     //Sets the element boxes for all elements in both coarse and fine models
//     SetElementBoxes();
//     MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

//     // // Rotating fine mesh
//     // for (int i = 0; i < fMeshVector[1]->NNodes(); ++i){
//     //     typename Nodes::VecLocD x, xn;
//     //     x = fMeshVector[1]->NodeVec()[i] -> getCoordinates();       
    
//     //     double a = -90 * pi / 180;

//     //     xn(0) = 0.0 + (x(0)-0.0) * cos(a) - (x(1)-0.0) * sin(a);
//     //     xn(1) = 0.0 + (x(0)-0.0) * sin(a) + (x(1)-0.0) * cos(a);

//     //     fMeshVector[1]->NodeVec()[i] -> setCoordinates(xn);
//     //     fMeshVector[1]->NodeVec()[i] -> setPreviousCoordinates(0,xn(0));
//     //     fMeshVector[1]->NodeVec()[i] -> setPreviousCoordinates(1,xn(1));
//     // }

//     parametersFine = &fMeshVector[1]->fProbParameters;
//     parametersCoarse = &fMeshVector[0]->fProbParameters;

//     // std::cout << "AA1 " << rank << std::endl;
//     // MPI_Barrier(PETSC_COMM_WORLD);

//     setSignaledDistance();

//     // std::cout << "AA2 " << rank << std::endl;
//     // MPI_Barrier(PETSC_COMM_WORLD);

//     //Construct the glue zone based on some defined criterion
//     setCouplingZone();
//     // std::cout << "AA3 " << rank << std::endl;
//     // MPI_Barrier(PETSC_COMM_WORLD);
//     //Computes the Weight function for all the finite elements
//     setWeightFunction(16.); 
//     // std::cout << "AA4 " << rank << std::endl;
//     // MPI_Barrier(PETSC_COMM_WORLD);
//     //Update domain decomposition - Start
//     int size;
//     MPI_Comm_size(PETSC_COMM_WORLD, &size);

//     Vec  b;
//     PetscErrorCode    ierr;
//     int start[size], end[size];
//     int numDOF = (DIM+1)*fMeshVector[0]->NNodes() + (DIM+1)*fMeshVector[1]->NNodes() + DIM*numNodesGlueZoneFine;

//     ierr = VecCreate(PETSC_COMM_WORLD,&b);
//     ierr = VecSetSizes(b,PETSC_DECIDE,numDOF);
//     VecSetFromOptions(b);
//     ierr = VecGetOwnershipRange(b,&start[rank],&end[rank]);

//     VecDestroy(&b);
    
//     for (int i = 0; i < size; ++i){
//         MPI_Bcast(&start[i],1,MPI_INT,i,PETSC_COMM_WORLD);
//         MPI_Bcast(&end[i],1,MPI_INT,i,PETSC_COMM_WORLD);
//     }
//     // std::cout << "AA5 " << rank << std::endl;
//     // MPI_Barrier(PETSC_COMM_WORLD);
//     for (int i = 0; i < fMeshVector[0]->NElements(); i++){
//         VecInt connec = fMeshVector[0]->ElementVec()[i] -> getConnectivity();
//         for (int j = 0; j < size; j++){
//             if (((DIM+1)*connec[0] >= start[j]) && ((DIM+1)*connec[0] <= end[j])) {
//                 domDecompCoarse.first[i] = j;
//                 break;
//             }
//         }
//     }    
//     for (int i = 0; i < fMeshVector[1]->NElements(); i++){
//         VecInt connec = fMeshVector[1]->ElementVec()[i] -> getConnectivity();
//         for (int j = 0; j < size; j++){
//             if (((DIM+1)*(fMeshVector[0]->NNodes() + connec[0]) >= start[j]) && ((DIM+1)*(fMeshVector[0]->NNodes() + connec[0]) <= end[j])) {
//                 domDecompFine.first[i] = j;
//                 break;
//             }
//         }
//     }    
//     // std::cout << "AA6 " << rank << std::endl;
//     // MPI_Barrier(PETSC_COMM_WORLD);
//     //Update domain decomposition - End

//     if(rank == 0){
//         std::cout << "---------------------ARLEQUIN DATA---------------------" 
//                   << std::endl;
//         std::cout << "Coarse Model: " << fMeshVector[0]->NNodes() << " nodes, " << fMeshVector[0]->NElements() << " elements." << std::endl;
//         std::cout << "Fine Model: " << fMeshVector[1]->NNodes() << " nodes, " << fMeshVector[1]->NElements() << " elements." << std::endl;
//         std::cout << "Lagrange Multipliers: " << numNodesGlueZoneFine << " nodes, " << numElemGlueZoneFine << " elements." << std::endl;
//         std::cout << "Number of Degrees of Freedom: " << numDOF << std::endl; 
//     }

//     if(rank == 0){
//         printResultsCoarse(10);
//         printResultsFine(10);
//     }


//     // fMeshVector[1]->NodeVec()[269] -> setWeightFunction(0.);
//     // fMeshVector[1]->NodeVec()[80] -> setWeightFunction(0.);





// };

//------------------------------------------------------------------------------
//----------------------COMPUTES DRAG AND LIFT COEFFICIENTS---------------------
//------------------------------------------------------------------------------

void Arlequin::dragAndLiftCoefficients(std::ofstream& dragLift){
    int DIM = fMeshVector[0]->Dimension();
    int DEG = fMeshVector[0]->GetDefaultOrder();
    double dragCoefficient = 0.;
    double liftCoefficient = 0.;
    double pressureDragCoefficient = 0.;
    double pressureLiftCoefficient = 0.;
    double frictionDragCoefficient = 0.;
    double frictionLiftCoefficient = 0.;
    double theta = 0.;
    double pitchingMomentCoefficient = 0.;
    double pMom = 0.;
    double per = 0.;

    double rhoInf = 1.0;
    double velocityInf[DIM];
    velocityInf[0] = 1;
    velocityInf[1] = 0.;

    for (int jel = 0; jel < fMeshVector[1]->NBoundElements(); jel++){   
        

        
        double dForce = 0.;
        double lForce = 0.;
        double pDForce = 0.;
        double pLForce = 0.;
        double fDForce = 0.;
        double fLForce = 0.;
        double aux_Mom = 0.;
        double aux_Per = 0.;



        // for (int i=0; i<fMeshVector[1]->numberOfLines; i++){
        //     if (boundaryFine_[jel] -> getBoundaryGroup() == fMeshVector[1]->dragAndLiftBoundary[i]){
                
        //         int iel = boundaryFine_[jel] -> getElement();
        //         // fMeshVector[1]->ElementVec()[iel] -> computeDragAndLiftForces(pDForce, pLForce, fDForce, fLForce, dForce, lForce, aux_Mom, aux_Per);

        //         pMom += aux_Mom;
        //         per += aux_Per;
        //     };
        // };
        
        pressureDragCoefficient += pDForce ;/// 
            //(0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        pressureLiftCoefficient += pLForce ;/// 
            //(0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        
        frictionDragCoefficient += fDForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0] * 1);
        frictionLiftCoefficient += fLForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0] * 1);
        
        dragCoefficient += dForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0] * 1);
        liftCoefficient += lForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0] * 1);
        
    };
    
    pitchingMomentCoefficient = pMom / (rhoInf * velocityInf[0] * velocityInf[0] * per);

    if (rank == 0) {
        const int timeWidth = 13;
        const int numWidth = 13;
        dragLift << std::setprecision(5) << std::scientific;
        dragLift << std::left << std::setw(timeWidth) << fMeshVector[1]->getProblemParameters().getTimeInstant() * fMeshVector[1]->getProblemParameters().GetTimeStep();
        dragLift << std::setw(numWidth) << pressureDragCoefficient;
        dragLift << std::setw(numWidth) << pressureLiftCoefficient;
        dragLift << std::setw(numWidth) << frictionDragCoefficient;
        dragLift << std::setw(numWidth) << frictionLiftCoefficient;
        dragLift << std::setw(numWidth) << dragCoefficient;
        dragLift << std::setw(numWidth) << liftCoefficient;
        dragLift << std::setw(numWidth) << pitchingMomentCoefficient;
        //dragLift << std::setw(numWidth) << theta;
        dragLift << std::endl;
    }

    return;
}

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------

void Arlequin::setMatVecValuesCoarseModel(MatrixDouble &matrix, VecDouble &rhs, VecInt &connec){
    int DIM = fMeshVector[0]->Dimension();
    int DEG = fMeshVector[0]->GetDefaultOrder();
    int nElNodes = fMeshVector[0]->NElNodes();
    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){
            for (int k = 0; k < DIM; k++){
                for (int l = 0; l < DIM; l++){
                    //Matrix K and C            
                    int dof_i = (DIM+1) * connec[i] + k;
                    int dof_j = (DIM+1) * connec[j] + l;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(DIM*i+k,DIM*j+l),ADD_VALUES);
                }
                //Matrix G and Gt
                int dof_i = (DIM+1) * connec[i] + k;
                int dof_j = (DIM+1) * connec[j] + DIM;
                MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(DIM*i+k,DIM*nElNodes+j),ADD_VALUES);
                MatSetValues(A,1,&dof_j,1,&dof_i,&matrix(DIM*nElNodes+j,DIM*i+k),ADD_VALUES);
            }
            // Matrix Q
            int dof_i = (DIM+1) * connec[i] + DIM;
            int dof_j = (DIM+1) * connec[j] + DIM;
            MatSetValues(A,1,&dof_j,1,&dof_i,&matrix(DIM*nElNodes+i,DIM*nElNodes+j),ADD_VALUES);
        };
        for (int k = 0; k < DIM; k++){          
            //Rhs vector
            int dof_i = (DIM+1) * connec[i] + k;
            VecSetValues(b,1,&dof_i,&rhs[DIM*i+k],ADD_VALUES);
        }
        int dof_i = (DIM+1) * connec[i] + DIM;
        VecSetValues(b,1,&dof_i,&rhs[DIM*nElNodes+i],ADD_VALUES);
    };

    return;
}

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------

void Arlequin::setMatVecValuesCoarseModelElasticity(MatrixDouble &matrix, VecDouble &rhs, VecInt &connec){
    int DIM = fMeshVector[0]->Dimension();
    int DEG = fMeshVector[0]->GetDefaultOrder();
    int nElNodes = fMeshVector[0]->NElNodes();
    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){
            for (int k = 0; k < DIM; k++){
                for (int l = 0; l < DIM; l++){
                    //Matrix K 
                    int dof_i = DIM * connec[i] + k;
                    int dof_j = DIM * connec[j] + l;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(DIM*i+k,DIM*j+l),ADD_VALUES);
                }
            }
        };
        for (int k = 0; k < DIM; k++){          
            //Rhs vector
            int dof_i = DIM * connec[i] + k;
            VecSetValues(b,1,&dof_i,&rhs[DIM*i+k],ADD_VALUES);
        }
    };

    return;
}

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------

void Arlequin::setMatVecValuesCoarseModelPoisson(MatrixDouble &matrix, VecDouble &rhs, VecInt &connec){
    int DIM = fMeshVector[0]->Dimension();
    int DEG = fMeshVector[0]->GetDefaultOrder();
    int nElNodes = fMeshVector[0]->NElNodes();
    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        int dof_i = connec[i];
        for (int j = 0; j < nElNodes; j++){
            //Matrix K
            int dof_j = connec[j];
            MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(i,j),ADD_VALUES);
        };
        //Rhs vector
        VecSetValues(b,1,&dof_i,&rhs[i],ADD_VALUES);
    };

    return;
}

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------

void Arlequin::setMatVecValuesFineModel(MatrixDouble &matrix, VecDouble &rhs, VecInt &connec){
    int DIM = fMeshVector[0]->Dimension();
    int DEG = fMeshVector[0]->GetDefaultOrder();
    int nElNodes = fMeshVector[1]->NElNodes();
    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){
            for (int k = 0; k < DIM; k++){
                for (int l = 0; l < DIM; l++){
                    //Matrix K and C
                    int dof_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * connec[i] + k;
                    int dof_j = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * connec[j] + l;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(DIM*i+k,DIM*j+l),ADD_VALUES);
                }          
                //Matrix G and Gt
                int dof_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * connec[i] + k;
                int dof_j = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * connec[j] + DIM;
                MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(DIM*i+k,DIM*nElNodes+j),ADD_VALUES);
                MatSetValues(A,1,&dof_j,1,&dof_i,&matrix(DIM*nElNodes+j,DIM*i+k),ADD_VALUES);
            }
            //Matrix Q
            int dof_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * connec[i] + DIM;
            int dof_j = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * connec[j] + DIM;
            MatSetValues(A,1,&dof_j,1,&dof_i,&matrix(DIM*nElNodes+i,DIM*nElNodes+j),ADD_VALUES);
        };
        for (int k = 0; k < DIM; k++){  
            ///Rhs vector
            int dof_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * connec[i] + k;
            VecSetValues(b,1,&dof_i,&rhs[DIM*i+k],ADD_VALUES);
        }
        int dof_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * connec[i] + DIM;
        VecSetValues(b,1,&dof_i,&rhs[DIM*nElNodes+i],ADD_VALUES);
    }; 

    return;
};

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------

void Arlequin::setMatVecValuesFineModelElasticity(MatrixDouble &matrix, VecDouble &rhs, VecInt &connec){
    int DIM = fMeshVector[0]->Dimension();
    int DEG = fMeshVector[0]->GetDefaultOrder();
    int nElNodes = fMeshVector[1]->NElNodes();
    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){
            for (int k = 0; k < DIM; k++){
                for (int l = 0; l < DIM; l++){
                    //Matrix K and C
                    int dof_i = DIM * fMeshVector[0]->NNodes() + DIM * connec[i] + k;
                    int dof_j = DIM * fMeshVector[0]->NNodes() + DIM * connec[j] + l;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(DIM*i+k,DIM*j+l),ADD_VALUES);
                }
            }
        };
        for (int k = 0; k < DIM; k++){  
            ///Rhs vector
            int dof_i = DIM * fMeshVector[0]->NNodes() + DIM * connec[i] + k;
            VecSetValues(b,1,&dof_i,&rhs[DIM*i+k],ADD_VALUES);
        }
    }; 

    return;
};

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------

void Arlequin::setMatVecValuesFineModelPoisson(MatrixDouble &matrix, VecDouble &rhs, VecInt &connec){
    int DIM = fMeshVector[0]->Dimension();
    int DEG = fMeshVector[0]->GetDefaultOrder();
    int nElNodes = fMeshVector[1]->NElNodes();
    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        int dof_i = fMeshVector[0]->NNodes() + connec[i];
        for (int j = 0; j < nElNodes; j++){
            //Matrix K
            int dof_j = fMeshVector[0]->NNodes() + connec[j];
            MatSetValues(A,1,&dof_i,1,&dof_j,&matrix(i,j),ADD_VALUES);
        };
        ///Rhs vector
        VecSetValues(b,1,&dof_i,&rhs[i],ADD_VALUES);
    }; 

    return;
};

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------

void Arlequin::setMatVecValuesLagMultFineFine(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                                       MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                                       VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                                       VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                                       VecInt &connec, VecInt &connecL){

    double &alpha_f = fMeshVector[1]->getProblemParameters().getAlphaF();
    double &alpha_m = fMeshVector[1]->getProblemParameters().getAlphaM();
    double &gamma = fMeshVector[1]->getProblemParameters().getGamma();
    int nElNodes = fMeshVector[1]->NElNodes();
    int DIM = fMeshVector[0]->Dimension();
    int DEG = fMeshVector[0]->GetDefaultOrder();
    double integ = alpha_f * gamma * fMeshVector[1]->getProblemParameters().GetTimeStep();
    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){
            for (int k = 0; k < DIM; k++){
                for (int l = 0; l < DIM; l++){
                    //COUPLING OPERATOR
                    if (fabs(Ajac2(DIM*i+k,DIM*j+l)) >= 1.e-15){
                        int d_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
                        int d_j = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * connec[j] + l;
                        double value = Ajac2(DIM*i+k,DIM*j+l)*integ;
                        MatSetValues(A,1,&d_i,1,&d_j,&value,ADD_VALUES);
                        MatSetValues(A,1,&d_j,1,&d_i,&Ajac2(DIM*i+k,DIM*j+l),ADD_VALUES);
                    };
                }
                //SUPG STABILIZATION
                if (fabs(localMV_mat(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    int d_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * fMeshVector[1]->NNodes() + DIM*connecL[i] + k;
                    int d_j = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * connec[j] + k;
                    MatSetValues(A,1,&d_j,1,&d_i,&localMV_mat(2*i  ,2*j  ),ADD_VALUES);
                };
                //PSPG STABILIZATION
                if (fabs(localMV_mat(DIM*i+k,DIM*nElNodes+j)) >= 1.e-15){
                    int dof_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * connec[i] + DIM;
                    int dof_j = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * fMeshVector[1]->NNodes() + connecL[j] + k;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&localMV_mat(DIM*i+k,DIM*nElNodes+j),ADD_VALUES);
                };
                //ARLEQUIN STABILIZATION
                if (fabs(ArlequinA2(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    int d_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
                    int d_j = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * connec[j] + k;
                    MatSetValues(A,1,&d_i,1,&d_j,&ArlequinA2(DIM*i+k,DIM*j+k),ADD_VALUES);
                };
                if (fabs(ArlequinA2(DIM*nElNodes+i,DIM*j+k)) >= 1.e-15){
                    int dof_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
                    int dof_j = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * connec[j] + DIM;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA2(DIM*nElNodes+i,DIM*j+k),ADD_VALUES);
                };
                if (fabs(ArlequinA1(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    int dof_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
                    int dof_j = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * fMeshVector[1]->NNodes() + DIM * connecL[j] + k;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA1(DIM*i+k,DIM*j+k),ADD_VALUES);
                };      
            }                      
        };
        //RHS VECTOR
        for (int k = 0; k < DIM; k++){
            //COUPLING OPERATOR
            int dof_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
            VecSetValues(b,1,&dof_i,&Rhs2[DIM*i+k],ADD_VALUES);

            dof_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * connec[i] + k;
            VecSetValues(b,1,&dof_i,&rhsLagMult2[DIM*i+k],ADD_VALUES);

            //SUPG STABILIZATION
            dof_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * connec[i] + k;
            VecSetValues(b,1,&dof_i,&localMV_vec[DIM*i+k],ADD_VALUES);

            //ARLEQUIN STABILIZATION
            dof_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * fMeshVector[1]->NNodes() + DIM*connecL[i] + k;
            VecSetValues(b,1,&dof_i,&RhsArlequin2[DIM*i+k],ADD_VALUES);
        }
        //PSPG STABILIZATION
        int dof_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * connec[i] + DIM;
        VecSetValues(b,1,&dof_i,&localMV_vec[DIM*nElNodes+i],ADD_VALUES);
    };


    return;
};

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------

void Arlequin::setMatVecValuesLagMultFineFineElasticity(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                                       MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                                       VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                                       VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                                       VecInt &connec, VecInt &connecL){

    int nElNodes = fMeshVector[1]->NElNodes();
    int DIM = fMeshVector[0]->Dimension();
    int DEG = fMeshVector[0]->GetDefaultOrder();
    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){
            for (int k = 0; k < DIM; k++){
                for (int l = 0; l < DIM; l++){
                    //COUPLING OPERATOR
                    if (fabs(Ajac2(DIM*i+k,DIM*j+l)) >= 1.e-15){
                        int d_i = DIM * fMeshVector[0]->NNodes() + DIM * fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
                        int d_j = DIM * fMeshVector[0]->NNodes() + DIM * connec[j] + l;
                        MatSetValues(A,1,&d_i,1,&d_j,&Ajac2(DIM*i+k,DIM*j+l),ADD_VALUES);
                        MatSetValues(A,1,&d_j,1,&d_i,&Ajac2(DIM*i+k,DIM*j+l),ADD_VALUES);
                    };
                }
                //ARLEQUIN STABILIZATION
                if (fabs(ArlequinA2(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    int d_i = DIM * fMeshVector[0]->NNodes() + DIM * fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
                    int d_j = DIM * fMeshVector[0]->NNodes() + DIM * connec[j] + k;
                    MatSetValues(A,1,&d_i,1,&d_j,&ArlequinA2(DIM*i+k,DIM*j+k),ADD_VALUES);
                };
                // if (fabs(ArlequinA2(DIM*nElNodes+i,DIM*j+k)) >= 1.e-15){
                //     int dof_i = DIM * fMeshVector[0]->NNodes() + DIM * fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
                //     int dof_j = DIM * fMeshVector[0]->NNodes() + DIM * connec[j] + DIM;
                //     MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA2(DIM*nElNodes+i,DIM*j+k),ADD_VALUES);
                // };
                if (fabs(ArlequinA1(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    int dof_i = DIM * fMeshVector[0]->NNodes() + DIM * fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
                    int dof_j = DIM * fMeshVector[0]->NNodes() + DIM * fMeshVector[1]->NNodes() + DIM * connecL[j] + k;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA1(DIM*i+k,DIM*j+k),ADD_VALUES);
                };      
            }                      
        };
        //RHS VECTOR
        for (int k = 0; k < DIM; k++){
            //COUPLING OPERATOR
            int dof_i = DIM * fMeshVector[0]->NNodes() + DIM * fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
            VecSetValues(b,1,&dof_i,&Rhs2[DIM*i+k],ADD_VALUES);

            dof_i = DIM * fMeshVector[0]->NNodes() + DIM * connec[i] + k;
            VecSetValues(b,1,&dof_i,&rhsLagMult2[DIM*i+k],ADD_VALUES);

            //ARLEQUIN STABILIZATION
            dof_i = DIM * fMeshVector[0]->NNodes() + DIM * fMeshVector[1]->NNodes() + DIM*connecL[i] + k;
            VecSetValues(b,1,&dof_i,&RhsArlequin2[DIM*i+k],ADD_VALUES);
        }
    };


    return;
};
//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------

void Arlequin::setMatVecValuesLagMultFineFinePoisson(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                                              MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                                              VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                                              VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                                              VecInt &connec, VecInt &connecL){

    int nElNodes = fMeshVector[1]->NElNodes();
    int DIM = fMeshVector[0]->Dimension();
    int DEG = fMeshVector[0]->GetDefaultOrder();                                                 
    //Disperse local contributions into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){
            //COUPLING OPERATOR
            if (fabs(Ajac2(i,j)) >= 1.e-15){
                int d_i = fMeshVector[0]->NNodes() + fMeshVector[1]->NNodes() + connecL[i];
                int d_j = fMeshVector[0]->NNodes() + connec[j];
                MatSetValues(A,1,&d_i,1,&d_j,&Ajac2(i,j),ADD_VALUES);
                MatSetValues(A,1,&d_j,1,&d_i,&Ajac2(i,j),ADD_VALUES);
            };
            // ARLEQUIN STABILIZATION
            if (fabs(ArlequinA2(i,j)) >= 1.e-15){
                int dof_i = fMeshVector[0]->NNodes() + fMeshVector[1]->NNodes() + connecL[i];
                int dof_j = fMeshVector[0]->NNodes() + connec[j];
                MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA2(i,j),ADD_VALUES);
            };
            if (fabs(ArlequinA1(i,j)) >= 1.e-15){
                int dof_i = fMeshVector[0]->NNodes() + fMeshVector[1]->NNodes() + connecL[i];
                int dof_j = fMeshVector[0]->NNodes() + fMeshVector[1]->NNodes() + connecL[j];
                MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA1(i,j),ADD_VALUES);
            };                   
        };
        //RHS VECTOR
        //COUPLING OPERATOR
        int dof_i = fMeshVector[0]->NNodes() + fMeshVector[1]->NNodes() + connecL[i];
        VecSetValues(b,1,&dof_i,&Rhs2[i],ADD_VALUES);

        dof_i = fMeshVector[0]->NNodes() + connec[i];
        VecSetValues(b,1,&dof_i,&rhsLagMult2[i],ADD_VALUES);

        //ARLEQUIN STABILIZATION
        dof_i = fMeshVector[0]->NNodes() + fMeshVector[1]->NNodes() + connecL[i];
        VecSetValues(b,1,&dof_i,&RhsArlequin2[i],ADD_VALUES);
    };


    return;
};

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------

void Arlequin::setMatVecValuesLagMultFineCoarse(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                                   MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                                   VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                                   VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                                   VecInt &connecC, VecInt &connecL){

    double &alpha_f = fMeshVector[1]->getProblemParameters().getAlphaF();
    double &alpha_m = fMeshVector[1]->getProblemParameters().getAlphaM();
    double &gamma = fMeshVector[1]->getProblemParameters().getGamma();
    double integ = alpha_f * gamma * fMeshVector[1]->getProblemParameters().GetTimeStep();
    int nElNodes = fMeshVector[1]->NElNodes(); 
    int DIM = fMeshVector[0]->Dimension();
    int DEG = fMeshVector[0]->GetDefaultOrder();

    //Disperse local contribution into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){
            for (int k = 0; k < DIM; k++){
                for (int l = 0; l < DIM; l++){
                    //COUPLING OPERATOR
                    if (fabs(Ajac2(DIM*i+k,DIM*j+l)) >= 1.e-15){
                        int dof_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
                        int dof_j = (DIM+1) * connecC[j] + l;
                        double value = Ajac2(DIM*i+k,DIM*j+l) * integ;
                        MatSetValues(A,1,&dof_i,1,&dof_j,&value,ADD_VALUES);
                        MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac2(DIM*i+k,DIM*j+l),ADD_VALUES);
                    }
                }
                //SUPG STABILIZATION
                if (fabs(localMV_mat(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    int dof_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
                    int dof_j = (DIM+1) * connecC[j] + k;
                    MatSetValues(A,1,&dof_j,1,&dof_i,&localMV_mat(DIM*i+k,DIM*j+k),ADD_VALUES);
                };
                //PSPG STABILIZATION
                if (fabs(localMV_mat(DIM*i+k,DIM*nElNodes+j)) >= 1.e-15){
                    int dof_i = (DIM+1) * connecC[i] + DIM;
                    int dof_j = (DIM+1) * fMeshVector[0]->NNodes() + DIM * fMeshVector[1]->NNodes() + connecL[j] + k;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&localMV_mat(DIM*i+k,DIM*nElNodes+j),ADD_VALUES);
                };
                //ARLEQUIN STABILIZATION
                if (fabs(ArlequinA2(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    int d_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
                    int d_j = (DIM+1) * connecC[j] + k;
                    MatSetValues(A,1,&d_i,1,&d_j,&ArlequinA2(2*i  ,2*j  ),ADD_VALUES);
                };
                if (fabs(ArlequinA2(DIM*nElNodes+i,DIM*j+k)) >= 1.e-15){
                    int dof_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * fMeshVector[1]->NNodes() + DIM * connecL[i] + 1;
                    int dof_j = (DIM+1) * connecC[j] + DIM;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA2(DIM*nElNodes+i,DIM*j+k),ADD_VALUES);
                };
                if (fabs(ArlequinA1(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    int dof_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1)*fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
                    int dof_j = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1)*fMeshVector[1]->NNodes() + DIM * connecL[j] + k;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA1(DIM*i+k,DIM*j+k),ADD_VALUES);
                };
            }
            // if (fabs(Ajac2(12+i,12+j)) >= 1.e-15){
            //     int dof_i = 3*fMeshVector[0]->NNodes() + 3*fMeshVector[1]->NNodes() + 2*connecL[i];
            //     int dof_j = 3*fMeshVector[0]->NNodes() + 3*fMeshVector[1]->NNodes() + 2*connecL[j];
            //     dof_i++; dof_j++;
            //     ierr = MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac2(12+i,12+j),ADD_VALUES);
            // };
        };
        for (int k = 0; k < DIM; k++){
            //RHS VECTOR
            //COUPLING OPERATOR
            int d_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
            VecSetValues(b,1,&d_i,&Rhs2[DIM*i+k],ADD_VALUES);

            int dof_i = (DIM+1) * connecC[i] + k;
            VecSetValues(b,1,&dof_i,&rhsLagMult2[DIM*i+k],ADD_VALUES);

            //SUPG STABILIZATION
            dof_i = (DIM+1) * connecC[i] + k;
            VecSetValues(b,1,&dof_i,&localMV_vec[DIM*i+k],ADD_VALUES);

            dof_i = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
            VecSetValues(b,1,&dof_i,&RhsArlequin2[DIM*i+k],ADD_VALUES);
        }
        //PSPG STABILIZATION
        int dof_i = (DIM+1) * connecC[i] + DIM;
        VecSetValues(b,1,&dof_i,&localMV_vec[DIM*nElNodes+i],ADD_VALUES);
    };
    return;
}

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------

void Arlequin::setMatVecValuesLagMultFineCoarseElasticity(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                                   MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                                   VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                                   VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                                   VecInt &connecC, VecInt &connecL){

    int nElNodes = fMeshVector[1]->NElNodes(); 
    int DIM = fMeshVector[0]->Dimension();
    int DEG = fMeshVector[0]->GetDefaultOrder();
    //Disperse local contribution into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){
            for (int k = 0; k < DIM; k++){
                for (int l = 0; l < DIM; l++){
                    //COUPLING OPERATOR
                    if (fabs(Ajac2(DIM*i+k,DIM*j+l)) >= 1.e-15){
                        int dof_i = DIM * fMeshVector[0]->NNodes() + DIM * fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
                        int dof_j = DIM * connecC[j] + l;
                        MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac2(DIM*i+k,DIM*j+l),ADD_VALUES);
                        MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac2(DIM*i+k,DIM*j+l),ADD_VALUES);
                    }
                }
                //ARLEQUIN STABILIZATION
                if (fabs(ArlequinA2(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    int d_i = DIM * fMeshVector[0]->NNodes() + DIM * fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
                    int d_j = DIM * connecC[j] + k;
                    MatSetValues(A,1,&d_i,1,&d_j,&ArlequinA2(2*i  ,2*j  ),ADD_VALUES);
                };
                // if (fabs(ArlequinA2(DIM*nElNodes+i,DIM*j+k)) >= 1.e-15){
                //     int dof_i = DIM * fMeshVector[0]->NNodes() + DIM * fMeshVector[1]->NNodes() + DIM * connecL[i] + 1;
                //     int dof_j = DIM * connecC[j] + DIM;
                //     MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA2(DIM*nElNodes+i,DIM*j+k),ADD_VALUES);
                // };
                if (fabs(ArlequinA1(DIM*i+k,DIM*j+k)) >= 1.e-15){
                    int dof_i = DIM * fMeshVector[0]->NNodes() + DIM*fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
                    int dof_j = DIM * fMeshVector[0]->NNodes() + DIM*fMeshVector[1]->NNodes() + DIM * connecL[j] + k;
                    MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA1(DIM*i+k,DIM*j+k),ADD_VALUES);
                };
            }
        };
        for (int k = 0; k < DIM; k++){
            //RHS VECTOR
            //COUPLING OPERATOR
            int d_i = DIM * fMeshVector[0]->NNodes() + DIM * fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
            VecSetValues(b,1,&d_i,&Rhs2[DIM*i+k],ADD_VALUES);

            int dof_i = DIM * connecC[i] + k;
            VecSetValues(b,1,&dof_i,&rhsLagMult2[DIM*i+k],ADD_VALUES);
            //Arlequin stabilization
            dof_i = DIM * fMeshVector[0]->NNodes() + DIM * fMeshVector[1]->NNodes() + DIM * connecL[i] + k;
            VecSetValues(b,1,&dof_i,&RhsArlequin2[DIM*i+k],ADD_VALUES);
        }
    };
    return;
}

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------

void Arlequin::setMatVecValuesLagMultFineCoarsePoisson(MatrixDouble &Ajac2, MatrixDouble &localMV_mat, 
                                                                MatrixDouble &ArlequinA1, MatrixDouble &ArlequinA2, 
                                                                VecDouble &Rhs2, VecDouble &rhsLagMult2,
                                                                VecDouble &localMV_vec, VecDouble &RhsArlequin2,
                                                                VecInt &connecC, VecInt &connecL){
                                                                    
    int nElNodes = fMeshVector[1]->NElNodes(); 
    int DIM = fMeshVector[0]->Dimension();
    int DEG = fMeshVector[0]->GetDefaultOrder();
    //Disperse local contribution into the global matrix
    for (int i = 0; i < nElNodes; i++){
        for (int j = 0; j < nElNodes; j++){
            //COUPLING OPERATOR
            if (fabs(Ajac2(i,j)) >= 1.e-15){
                int dof_i = fMeshVector[0]->NNodes() + fMeshVector[1]->NNodes() + connecL[i];
                int dof_j = connecC[j];
                MatSetValues(A,1,&dof_i,1,&dof_j,&Ajac2(i,j),ADD_VALUES);
                MatSetValues(A,1,&dof_j,1,&dof_i,&Ajac2(i,j),ADD_VALUES);
            }
             
            //ARLEQUIN STABILIZATION
            if (fabs(ArlequinA2(i,j)) >= 1.e-15){
                int dof_i = fMeshVector[0]->NNodes() + fMeshVector[1]->NNodes() + connecL[i];
                int dof_j = connecC[j];
                MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA2(i,j),ADD_VALUES);
            };
            if (fabs(ArlequinA1(i,j)) >= 1.e-15){
                int dof_i = fMeshVector[0]->NNodes() + fMeshVector[1]->NNodes() + connecL[i];
                int dof_j = fMeshVector[0]->NNodes() + fMeshVector[1]->NNodes() + connecL[j];
                MatSetValues(A,1,&dof_i,1,&dof_j,&ArlequinA1(i,j),ADD_VALUES);
            };
        };
        
        //RHS VECTOR
        //COUPLING OPERATOR
        int d_i = fMeshVector[0]->NNodes() + fMeshVector[1]->NNodes() + connecL[i];
        VecSetValues(b,1,&d_i,&Rhs2[i],ADD_VALUES);

        int dof_i = connecC[i];
        VecSetValues(b,1,&dof_i,&rhsLagMult2[i],ADD_VALUES);


        dof_i = fMeshVector[0]->NNodes() + fMeshVector[1]->NNodes() + connecL[i];
        VecSetValues(b,1,&dof_i,&RhsArlequin2[i],ADD_VALUES);
    };
    return;
}

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------

void Arlequin::assembleCoarseModel(){

    //Coarse mesh
    for (int jel = 0; jel < fMeshVector[0]->NElements(); jel++){   
        if (fMeshVector[0]->part_elem[jel] == rank) {
            int nLocDOF = fMeshVector[0]->NLocDOF();
            //Compute Element matrix
            MatrixDouble matrix(nLocDOF,nLocDOF);
            matrix.setZero();
            VecDouble rhs(nLocDOF);
            rhs.setZero();
            
            fMeshVector[0]->ElementVec()[jel] -> ComputeElContribution(matrix,rhs);

            switch (fProbType){
            case ENavierStokes:
                setMatVecValuesCoarseModel(matrix,rhs,fMeshVector[0]->ElementVec()[jel] -> getConnectivity());
                break;
            case EPoisson:
                setMatVecValuesCoarseModelPoisson(matrix,rhs,fMeshVector[0]->ElementVec()[jel] -> getConnectivity());
                break;
            case EElastic:
                setMatVecValuesCoarseModelElasticity(matrix,rhs,fMeshVector[0]->ElementVec()[jel] -> getConnectivity());
                break;
            
            default:
                PanicButton();
                break;
            }
        };
    };
}


void Arlequin::assembleFineModel(){

    //Fine mesh
    for (int jel = 0; jel < fMeshVector[1]->NElements(); jel++){           
        if (fMeshVector[1]->part_elem[jel] == rank) {
            int nLocDOF = fMeshVector[1]->NLocDOF();
            //Compute Element matrix                    
            MatrixDouble matrix(nLocDOF,nLocDOF);
            matrix.setZero();
            VecDouble rhs(nLocDOF);
            rhs.setZero();
            fMeshVector[1]->ElementVec()[jel] -> ComputeElContribution(matrix,rhs);

            switch (fProbType){
            case ENavierStokes:
                setMatVecValuesFineModel(matrix,rhs,fMeshVector[1]->ElementVec()[jel] -> getConnectivity());
                break;
            case EPoisson:
                setMatVecValuesFineModelPoisson(matrix,rhs,fMeshVector[1]->ElementVec()[jel] -> getConnectivity());
                break;
            case EElastic:
                setMatVecValuesFineModelElasticity(matrix,rhs,fMeshVector[1]->ElementVec()[jel] -> getConnectivity());
                break;
            
            default:
                PanicButton();
                break;
            }
            

        };                
    };

}


void Arlequin::assembleCouplingOperator(){
    // int DIM = fMeshVector[0]->Dimension();
    // int DEG = fMeshVector[0]->GetDefaultOrder();
    // //Lagrange Multipliers
    // for (int l=0; l< numElemGlueZoneFine; l++){
    //     int jel = elementsGlueZoneFine_[l];
    //     if (fMeshVector[1]->part_elem[jel] == rank) {
            
    //         int nElNodes = fMeshVector[1]->NElNodes();
    //         int dofMatrices = 0;
    //         if (fProbType == EPoisson){
    //             dofMatrices = nElNodes;
    //         } else if (fProbType == EElastic) {
    //             dofMatrices = nElNodes * DIM;
    //         } else {
    //             dofMatrices = nElNodes * (DIM+1);
    //         }

    //         VecInt connecC;
    //         VecInt connec = fMeshVector[1]->ElementVec()[jel] -> getConnectivity();
    //         VecInt connecL = glueZoneFine_[l] -> getConnectivity();
            
    //         //FINE MESH
    //         //Matrices
    //         MatrixDouble matC0(dofMatrices,dofMatrices),matC1(dofMatrices,dofMatrices);
    //         MatrixDouble matA0(dofMatrices,dofMatrices),matA1(dofMatrices,dofMatrices);
    //         MatrixDouble localMV_mat(dofMatrices,dofMatrices);
    //         MatrixDouble matE0(dofMatrices,dofMatrices),matE1(dofMatrices,dofMatrices);
    //         matC1.setZero();
    //         matA1.setZero();
    //         localMV_mat.setZero();
    //         matE1.setZero();

    //         double tARLQ0_, tARLQ1_;
            
    //         //Vectors
    //         VecDouble vecC0(dofMatrices),vecC1(dofMatrices);
    //         VecDouble vecU0(dofMatrices),vecU1(dofMatrices);
    //         VecDouble RhsA0(dofMatrices),RhsA1(dofMatrices);
    //         VecDouble localMV_vec(dofMatrices);
    //         VecDouble vecE0(dofMatrices),vecE1(dofMatrices);
    //         vecC1.setZero();
    //         vecU1.setZero();
    //         localMV_vec.setZero();
    //         vecE1.setZero();

    //         // FINE MESH
    //         //Computes element matrix
    //         fMeshVector[1]->ElementVec()[jel] -> getLagrangeMultipliersSameMesh(matC1, vecC1, vecU1);
            
    //         if (fArlequinStab != ArlequinStabType::ENoStab){
    //             if (fProbType == ENavierStokes || fProbType == EStokes){
    //                 //PSPG and SUPG stabilizations
    //                 fMeshVector[1]->ElementVec()[jel] -> getLagrangeMultipliersSUPG_PSPG_SameMesh(localMV_mat,localMV_vec);
    //             }
    //             //Arlequin Stabilization
    //             fMeshVector[1]->ElementVec()[jel] -> getLagrangeMultipliersArlequinSameMesh(matE1, matA1, vecE1);
    //         }

    //         //COAESE MESH
    //         //Counts number of coarse mesh intersecting the fine element
    //         int numberIntPoints = fMeshVector[1]->ElementVec()[jel] -> 
    //             getNumberOfIntegrationPoints();
    //         int aux;
            
    //         std::vector<int> ele, diffElem;
    //         ele.clear();
    //         diffElem.clear();

    //         ele.reserve(3);
    //         for (int i=0; i<numberIntPoints; i++){
    //             aux = fMeshVector[1]->ElementVec()[jel] -> 
    //                 getIntegPointCorrespondenceElement(i);
    //             ele.push_back(aux);
    //             //std::cout << "Num elem inters " << jel << " " << aux << std::endl;
    //         };
            
    //         int numElemIntersect = 1;
    //         int flag = 0;
    //         diffElem.push_back(ele[0]);
            
    //         for (int i = 1; i<numberIntPoints; i++){
    //             flag = 0;
    //             for (int j = 0; j<numElemIntersect; j++){
    //                 if (ele[i] == diffElem[j]) {
    //                     break;
    //                 }else{
    //                     flag++;
    //                 };
    //                 if(flag == numElemIntersect){
    //                     numElemIntersect++;
    //                     diffElem.push_back(ele[i]);
    //                 };
    //             };
    //         };
    //         //Compute the Lagrange Multiplier element matrix
    //         for (int ielem = 0; ielem < numElemIntersect; ielem++){
    //             int nElNodes = fMeshVector[1]->NElNodes();
    //             int iElemCoarse = diffElem[ielem];
    //             double pspg = 0;//(*fMeshVector[0]->ElementVec()[iElemCoarse]) -> getPSPG();
    //             VecDouble press_(nElNodes), velX_(nElNodes), velY_(nElNodes), velXPrev_(nElNodes), velYPrev_(nElNodes);

    //             connecC = fMeshVector[0]->ElementVec()[iElemCoarse] -> getConnectivity();

    //             for (int k = 0; k < nElNodes; k++){
    //                 press_[k] = fMeshVector[0]->NodeVec()[connecC[k]] -> getPressure();
    //                 velX_[k] = fMeshVector[0]->NodeVec()[connecC[k]] -> getVelocity(0);
    //                 velY_[k] = fMeshVector[0]->NodeVec()[connecC[k]] -> getVelocity(1);
    //                 velXPrev_[k] = fMeshVector[0]->NodeVec()[connecC[k]] -> getPreviousVelocity(0);
    //                 velYPrev_[k] = fMeshVector[0]->NodeVec()[connecC[k]] -> getPreviousVelocity(1);
    //             }
                
    //             matC0.setZero();
    //             localMV_mat.setZero();
    //             matE0.setZero();
                
    //             //Vectors
    //             vecC0.setZero();
    //             vecU0.setZero();
    //             localMV_vec.setZero();
    //             vecE0.setZero();
    //             matA0.setZero();
                
    //             fMeshVector[1]->ElementVec()[jel] -> getLagrangeMultipliersDifferentMesh(iElemCoarse,pspg,press_,velX_,velY_,velXPrev_,velYPrev_,matC0,vecC0,vecU0);

    //             if (fArlequinStab != ArlequinStabType::ENoStab){
    //                 if (fProbType == ENavierStokes || fProbType == EStokes){
    //                     fMeshVector[1]->ElementVec()[jel] -> getLagrangeMultipliersSUPG_PSPG_DifferentMesh(iElemCoarse,pspg,press_,velX_,velY_,localMV_mat,localMV_vec);
    //                 }

    //                 fMeshVector[1]->ElementVec()[jel] -> getLagrangeMultipliersArlequinDifferentMesh(iElemCoarse,pspg,press_,velX_,velY_,matE0,matA0,vecE0);
    //             }
    //             MatrixDouble matE = matE1;
    //             stabilizeArlequin(matA0,matA1,matC0,matC1,matE0,vecE0,vecE1,tARLQ0_,tARLQ1_);

    //             matE0 *= tARLQ0_;
    //             matA0 *= tARLQ0_;
    //             vecE0 *= tARLQ0_;
                
    //             switch (fProbType)
    //             {
    //             case ENavierStokes:
    //             case EStokes:
    //                 setMatVecValuesLagMultFineCoarse(matC0, localMV_mat, matE0, matA0, 
    //                                              vecU0, vecC0, localMV_vec, vecE0,
    //                                              fMeshVector[0]->ElementVec()[iElemCoarse] -> getConnectivity(), 
    //                                              glueZoneFine_[l] -> getConnectivity());
    //                 break;
    //             case EPoisson:
    //                 setMatVecValuesLagMultFineCoarsePoisson(matC0, localMV_mat, matE0, matA0, 
    //                                                     vecU0, vecC0, localMV_vec, vecE0,
    //                                                     fMeshVector[0]->ElementVec()[iElemCoarse] -> getConnectivity(), 
    //                                                     glueZoneFine_[l] -> getConnectivity());
    //                 break;
    //             case EElastic:
    //                 setMatVecValuesLagMultFineCoarseElasticity(matC0, localMV_mat, matE0, matA0, 
    //                                              vecU0, vecC0, localMV_vec, vecE0,
    //                                              fMeshVector[0]->ElementVec()[iElemCoarse] -> getConnectivity(), 
    //                                              glueZoneFine_[l] -> getConnectivity());
    //                 break;
                
    //             default:
    //                 PanicButton();
    //                 break;
    //             }
    //         }; //Number of intersections

    //         matE1 *= tARLQ1_;
    //         matA1 *= tARLQ1_;
    //         vecE1 *= tARLQ1_;

    //         switch (fProbType)
    //         {
    //         case ENavierStokes:
    //         case EStokes:
    //             setMatVecValuesLagMultFineFine(matC1,localMV_mat,matE1,matA1, 
    //                                     vecU1,vecC1,localMV_vec,vecE1,
    //                                     fMeshVector[1]->ElementVec()[jel] -> getConnectivity(),
    //                                     glueZoneFine_[l] -> getConnectivity());
    //             break;
    //         case EPoisson:
    //             setMatVecValuesLagMultFineFinePoisson(matC1,localMV_mat,matE1,matA1, 
    //                                                   vecU1,vecC0,localMV_vec,vecE1,
    //                                                   fMeshVector[1]->ElementVec()[jel] -> getConnectivity(),
    //                                                   glueZoneFine_[l] -> getConnectivity());
    //             break;
    //         case EElastic:
    //             setMatVecValuesLagMultFineFineElasticity(matC1,localMV_mat,matE1,matA1, 
    //                                                   vecU1,vecC0,localMV_vec,vecE1,
    //                                                   fMeshVector[1]->ElementVec()[jel] -> getConnectivity(),
    //                                                   glueZoneFine_[l] -> getConnectivity());
    //             break;
            
    //         default:
    //             PanicButton();
    //             break;
    //         }

            
           
    //     }; // if element belongs to the glue zone
    // }; // Glue zone
}


void Arlequin::assembleArlequinSystem(){

    assembleCoarseModel();
    assembleFineModel();
    assembleCouplingOperator();
   
}

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------

int Arlequin::solveArlequinProblem(int iterNumber, double tolerance,
                                            int problem_type, int time_dependency){


//     std::ofstream dragLift;
//     int DIM = fMeshVector[0]->Dimension();
//     int DEG = fMeshVector[0]->GetDefaultOrder();
//     dragLift.open("dragLift.dat", std::ofstream::out | std::ofstream::app);
//     if (rank == 0) {
//         dragLift << "Time   Pressure Drag   Pressure Lift " 
//                  << "Friction Drag  Friction Lift Drag    Lift " 
//                  << std::endl;
//     };

//     if ((problem_type < 1) || (problem_type > 2)){
//         std::cout << "WRONG PROBLEM TYPE." << std::endl;
//         return 0;
//     };

//     if (time_dependency == 0) fMeshVector[1]->getProblemParameters().GetNTimeSteps() = 1;
    
//     fMeshVector[1]->getProblemParameters().setTimeInstant(0);
   
//     // //Computes the Weight function for all the finite elements
//     setWeightFunction(16.);

//     //Computes the Nodal correspondence between fine nodes and coarse elements
//     setNodalCorrespondenceFine();
//     if (rank == 0) {
//         printResultsCoarse(100);
//         printResultsFine(100);
//     }
//     // Computes the system size
//     int sysSize;
//     if (fProbType == ProblemType::ENavierStokes || fProbType == ProblemType::EStokes){
//         sysSize = (DIM+1)*fMeshVector[0]->NNodes() + (DIM+1)*fMeshVector[1]->NNodes() + DIM*numNodesGlueZoneFine;
//     } else if (fProbType == ProblemType::EPoisson){
//         sysSize = fMeshVector[0]->NNodes() + fMeshVector[1]->NNodes() + numNodesGlueZoneFine;
//     } else if (fProbType == ProblemType::EElastic){
//         sysSize = (fMeshVector[0]->NNodes() + fMeshVector[1]->NNodes() + numNodesGlueZoneFine)*DIM;
//     }
//     double integScheme = fMeshVector[1]->getProblemParameters().GetSpectralRadius();

//     double alpha_f = 1. / (1. + integScheme);
//     double alpha_m = 0.5 * (3. - integScheme) / (1. + integScheme);
//     double gamma = 0.5 + alpha_m - alpha_f;

//     for (int iTimeStep = 0; iTimeStep < fMeshVector[1]->getProblemParameters().GetNTimeSteps(); iTimeStep++){
        
//         fMeshVector[0]->getProblemParameters().setTimeInstant(iTimeStep);
//         fMeshVector[1]->getProblemParameters().setTimeInstant(iTimeStep);

//         if (rank == 0) {std::cout << "----------------------------" 
//                                   << " TIME STEP = "
//                                   << iTimeStep 
//                                   << " ---------------------------"
//                                   << std::endl;}
//         PetscMemoryGetCurrentUsage(&bytes);
//         PetscPrintf(PETSC_COMM_WORLD,"Memory used %g M\n",bytes/(1024*1024));

//         //Updates velocity and acceleration
//         for (int i = 0; i < fMeshVector[0]->NNodes(); i++){
//             VecDouble accel(DIM), u(DIM), uprev(DIM);
            
//             //Compute acceleration
//             u[0] = fMeshVector[0]->NodeVec()[i] -> getVelocity(0);
//             u[1] = fMeshVector[0]->NodeVec()[i] -> getVelocity(1);

//             fMeshVector[0]->NodeVec()[i] -> setPreviousVelocity(u);
            
//             accel[0] = fMeshVector[0]->NodeVec()[i] -> getAcceleration(0);
//             accel[1] = fMeshVector[0]->NodeVec()[i] -> getAcceleration(1);
            
//             fMeshVector[0]->NodeVec()[i] -> setPreviousAcceleration(accel);

//             accel[0] *= (gamma - 1.) / gamma;
//             accel[1] *= (gamma - 1.) / gamma;
            
//             fMeshVector[0]->NodeVec()[i] -> setAcceleration(accel);            

//         };

//         for (int i = 0; i < fMeshVector[1]->NNodes(); i++){
//             VecDouble accel(DIM), u(DIM), uprev(DIM), lag(DIM);
            
//             //Compute acceleration
//             u[0] = fMeshVector[1]->NodeVec()[i] -> getVelocity(0);
//             u[1] = fMeshVector[1]->NodeVec()[i] -> getVelocity(1);

//             fMeshVector[1]->NodeVec()[i] -> setPreviousVelocity(u);
            
//             accel[0] = fMeshVector[1]->NodeVec()[i] -> getAcceleration(0);
//             accel[1] = fMeshVector[1]->NodeVec()[i] -> getAcceleration(1);
            
//             fMeshVector[1]->NodeVec()[i] -> setPreviousAcceleration(accel);

//             accel[0] *= (gamma - 1.) / gamma;
//             accel[1] *= (gamma - 1.) / gamma;
            
//             fMeshVector[1]->NodeVec()[i] -> setAcceleration(accel);
//         };

//         //STARTS NEWTON-RAPHSON
//         for (int inewton = 0; inewton < iterNumber; inewton++){
            
//             std::clock_t t1 = std::clock();
            
//             // Preallocates the matrix
//             if (fMeshVector[1]->getProblemParameters().getSolverType() == SolverType::ESuiteSparse){
//                 ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, sysSize, sysSize, 100,NULL,&A);
//             } else {
//                 ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
//                                 sysSize, sysSize,400,NULL,600,NULL,&A); 
//             }

//             // for (int i=0; i<sysSize; i++){
//             //     double val = 1.e-10;
//             //     ierr = MatSetValues(A,1,&i,1,&i,&val,ADD_VALUES);   
//             // }  
//             CHKERRQ(ierr);
            
//             // Divides the matrix between the processes
//             ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
            
//             //Create PETSc vectors
//             ierr = VecCreate(PETSC_COMM_WORLD, &b); CHKERRQ(ierr);
//             ierr = VecSetSizes(b, PETSC_DECIDE, sysSize); CHKERRQ(ierr);
//             ierr = VecSetFromOptions(b); CHKERRQ(ierr); 
//             ierr = VecDuplicate(b, &u); CHKERRQ(ierr);
                        
//             for (int i=0; i<sysSize; i++){
//                 double val = 1.e-20;
//                 ierr = MatSetValues(A,1,&i,1,&i,&val,ADD_VALUES);
                
//             }
            

//             std::clock_t t3 = std::clock();
//             assembleArlequinSystem();
//             // if (fProbType == ProblemType::ENavierStokes){
//             //     assembleArlequinSystem();
//             // } else if (fProbType == ProblemType::EPoisson){
//             //     assembleArlequinSystemPoisson();
//             // }
            
//             std::clock_t t4 = std::clock();

//             //std::cout << "Enter PETSc " << rank << std::endl;
            
//             //Assemble matrices and vectors
//             ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//             ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            
//             ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
//             ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
//             std::clock_t t5 = std::clock();

          
            
//  // PetscViewer    viewer;

//  // PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,300,300,&viewer);
//  // PetscObjectSetName((PetscObject)viewer,"Line graph Plot");
//  //  PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);

//             // ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//             // ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            
//             //Create KSP context to solve the linear system
//             ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
            
//             ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);


//             switch (fMeshVector[1]->getProblemParameters().getSolverType())
//             {
//             case SolverType::ESuiteSparse:
//                 KSPGetPC(ksp, &pc);
//                 PCSetType(pc, PCLU);
//                 PCFactorSetMatSolverType(pc, MATSOLVERUMFPACK);
//                 break;
//             case SolverType::EMumps:
//                 KSPGetPC(ksp, &pc);
//                 PCSetType(pc, PCLU);
//                 PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
//                 break;

//             case SolverType::EIterative:
//                 KSPSetType(ksp,KSPFGMRES);
//                 KSPGetPC(ksp, &pc);
//                 PCSetType(pc,PCBJACOBI);
//                 KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,200);
//                 break;

//             default:
//                 PanicButton();
//                 break;
//             }


            
//             ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
            
//             ierr = KSPGetTotalIterations(ksp, &iterations);

//             std::clock_t t6 = std::clock();

//             if (rank == 0) std::cout << "TIME " << 1000.*(t4-t3)/CLOCKS_PER_SEC/1000. << " " 
//                                                 << 1000.*(t5-t4)/CLOCKS_PER_SEC/1000. << " " 
//                                                 << 1000.*(t6-t5)/CLOCKS_PER_SEC/1000. << std::endl;
//             //if (rank == 0)std::cout << "GMRES Iterations = " << iterations << std::endl;
        
//             //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            
//             //Gathers the solution vector to the master process
//             ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
            
//             ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
//             CHKERRQ(ierr);
            
//             ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
//             CHKERRQ(ierr);
            
//             ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
            
//             //Updates nodal values
//             double u_[DIM];
//             double normU = 0.;
//             double normP = 0.;
//             double normL = 0.;
//             double p_;
//             Ione = 1;

//             if (fProbType == ENavierStokes || fProbType == EStokes){
//                 for (int i = 0; i < fMeshVector[0]->NNodes(); ++i){
//                     double w_ = fMeshVector[0]->NodeVec()[i] -> getWeightFunction();
//                     for (int k = 0; k < DIM; k++){
//                         Ii = (DIM+1) * i + k;
//                         ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                         // if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
//                         u_[k] = val;
//                         normU += val*w_*val*w_;
//                         fMeshVector[0]->NodeVec()[i] -> incrementAcceleration(k,u_[k]);
//                         fMeshVector[0]->NodeVec()[i] -> incrementVelocity(k,u_[k]*gamma*fMeshVector[1]->getProblemParameters().GetTimeStep());
//                     }
//                     Ii = (DIM+1) * i + DIM;
//                     ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
//                     p_ = val;
//                     normP += val*w_*val*w_;
//                     fMeshVector[0]->NodeVec()[i] -> incrementPressure(p_);
//                 };
                
//                 for (int i = 0; i < fMeshVector[1]->NNodes(); ++i){
//                     double w_ = fMeshVector[1]->NodeVec()[i] -> getWeightFunction();
//                     for (int k = 0; k < DIM; k++){
//                         Ii = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * i + k;
//                         ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                         u_[k] = val;
//                         normU += val*w_*val*w_;
//                         fMeshVector[1]->NodeVec()[i] -> incrementAcceleration(k,u_[k]);
//                         fMeshVector[1]->NodeVec()[i] -> incrementVelocity(k,u_[k]*gamma*fMeshVector[1]->getProblemParameters().GetTimeStep());
//                     }        
//                     Ii = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * i + DIM;
//                     ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
//                     p_ = val;
//                     normP += val*w_*val*w_;
//                     fMeshVector[1]->NodeVec()[i] -> incrementPressure(p_);
//                 };
                
//                 for (int i = 0; i < numNodesGlueZoneFine; ++i){
//                     for (int k = 0; k < DIM; k++){
//                         Ii = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * fMeshVector[1]->NNodes() + DIM * i + k;
//                         ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                         u_[k] = val;
//                         fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> incrementLagrangeMultiplier(k,u_[k]);
//                         normL += val*val;
//                     }
//                 };
//             } else if (fProbType == EPoisson){
//                 for (int i = 0; i < fMeshVector[0]->NNodes(); ++i){
//                     double w_ = fMeshVector[0]->NodeVec()[i] -> getWeightFunction();
//                     Ii = i;
//                     ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                     // if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
//                     u_[0] = val;
//                     normU += val*w_*val*w_;
//                     fMeshVector[0]->NodeVec()[i] -> incrementVelocity(0,u_[0]);
//                 };
                
//                 for (int i = 0; i < fMeshVector[1]->NNodes(); ++i){
//                     double w_ = fMeshVector[1]->NodeVec()[i] -> getWeightFunction();
//                     Ii = fMeshVector[0]->NNodes() + i;
//                     ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                     u_[0] = val;
//                     normU += val*w_*val*w_;
//                     fMeshVector[1]->NodeVec()[i] -> incrementVelocity(0,u_[0]);
//                 };
                
//                 for (int i = 0; i < numNodesGlueZoneFine; ++i){
//                     Ii = fMeshVector[0]->NNodes() + fMeshVector[1]->NNodes() + i;
//                     ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                     u_[0] = val;
//                     fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> incrementLagrangeMultiplier(0,u_[0]);
//                     normL += val*val;

//                 };
//             } else if (fProbType == EElastic){
//                 for (int i = 0; i < fMeshVector[0]->NNodes(); ++i){
//                     double w_ = fMeshVector[0]->NodeVec()[i] -> getWeightFunction();
//                     for (int k = 0; k < DIM; k++){
//                         Ii = DIM * i + k;
//                         ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                         // if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
//                         u_[k] = val;
//                         normU += val*w_*val*w_;
//                         fMeshVector[0]->NodeVec()[i] -> incrementAcceleration(k,u_[k]);
//                         fMeshVector[0]->NodeVec()[i] -> incrementVelocity(k,u_[k]*gamma*fMeshVector[1]->getProblemParameters().GetTimeStep());
//                     }
//                 };
                
//                 for (int i = 0; i < fMeshVector[1]->NNodes(); ++i){
//                     double w_ = fMeshVector[1]->NodeVec()[i] -> getWeightFunction();
//                     for (int k = 0; k < DIM; k++){
//                         Ii = DIM * fMeshVector[0]->NNodes() + DIM * i + k;
//                         ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                         u_[k] = val;
//                         normU += val*w_*val*w_;
//                         fMeshVector[1]->NodeVec()[i] -> incrementAcceleration(k,u_[k]);
//                         fMeshVector[1]->NodeVec()[i] -> incrementVelocity(k,u_[k]*gamma*fMeshVector[1]->getProblemParameters().GetTimeStep());
//                     }        
//                 };
                
//                 for (int i = 0; i < numNodesGlueZoneFine; ++i){
//                     for (int k = 0; k < DIM; k++){
//                         Ii = DIM * fMeshVector[0]->NNodes() + DIM * fMeshVector[1]->NNodes() + DIM * i + k;
//                         ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                         u_[k] = val;
//                         fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> incrementLagrangeMultiplier(k,u_[k]);
//                         normL += val*val;
//                     }
//                 };

//             } else{
//                 PanicButton();
//             }
            
//             //Computes the solution vector norm
//             ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
            
//             std::clock_t t2 = std::clock();
            
//             if(rank == 0){
//                 std::cout<<"Iteration = " << inewton << " (" << iterations <<  
//                     ")  Du Norm = " << std::scientific << sqrt(normU) 
//                          << " " << sqrt(normP) 
//                          << " " << sqrt(normL) 
//                          << " " << val << 
//                     "  Time (s) = " << std::fixed << 
//                     1000.*(t2-t1)/CLOCKS_PER_SEC/1000. << std::endl;
//             };
            
//             ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
//             ierr = VecDestroy(&b); CHKERRQ(ierr);
//             ierr = VecDestroy(&u); CHKERRQ(ierr);
//             ierr = VecDestroy(&All); CHKERRQ(ierr);
//             ierr = MatDestroy(&A); CHKERRQ(ierr);
//            // ierr = MatDestroy(&F); CHKERRQ(ierr);
//             //ierr = MatDestroy(&C); CHKERRQ(ierr);
            
//             if(val <= tolerance){
//                 break;            
//             }; 
//         };
        
//         if (rank == 0){
//             double normUUprev, normU;
//             normUUprev = 0.;
//             normU = 0.;
//             for (int i = 0; i < fMeshVector[0]->NNodes(); ++i){
//                 for (int k = 0; k < DIM; k++){
//                     double u, uPr, weight;
//                     weight = fMeshVector[0]->NodeVec()[i] -> getWeightFunction();
//                     u = fMeshVector[0]->NodeVec()[i] -> getVelocity(k) * weight;
//                     uPr = fMeshVector[0]->NodeVec()[i] -> getPreviousVelocity(k) * weight;
//                     normUUprev += (u-uPr) * (u-uPr);
//                     normU += u*u;
//                 }
//             }
//             for (int i = 0; i < fMeshVector[1]->NNodes(); ++i){
//                 for (int k = 0; k < DIM; k++){
//                     double u, uPr, weight;
//                     weight = fMeshVector[1]->NodeVec()[i] -> getWeightFunction();
//                     u = fMeshVector[1]->NodeVec()[i] -> getVelocity(k) * weight;
//                     uPr = fMeshVector[1]->NodeVec()[i] -> getPreviousVelocity(k) * weight;
//                     normUUprev += (u-uPr) * (u-uPr);
//                     normU += u*u;
//                 }
//             }
//             std::cout << "NORM U  " << std::scientific <<  sqrt(normUUprev/normU) << std::endl;
//         }

//         //Compute real velocity
//         ShapeFunction shapeQuad(DIM,DEG);
//         int nElNodes = fMeshVector[1]->NElNodes(); 
//         VecDouble phi_(nElNodes);
        
//         for (int i = 0; i<fMeshVector[1]->NNodes(); i++){
//             for (int k = 0; k < DIM; k++) 
//                 fMeshVector[1]->NodeVec()[i] -> setVelocityArlequin(k,fMeshVector[1]->NodeVec()[i] -> getVelocity(k));
//             fMeshVector[1]->NodeVec()[i] -> setPressureArlequin(fMeshVector[1]->NodeVec()[i] ->getPressure());
//         };
        
//         for (int i = 0; i<numNodesGlueZoneFine; i++){
//             double u_coarse[nElNodes], v_coarse[nElNodes], p_coarse[nElNodes];
            
//             double u = 0.;
//             double v = 0.;
//             double p = 0.;
            
//             int elCoarse = fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> getNodalElemCorrespondence();
//             VecDouble xsi = fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> getNodalXsiCorrespondence();
                        
//             VecInt connecCoarse = fMeshVector[0]->ElementVec()[elCoarse] -> getConnectivity();
            
//             for (int j=0; j<nElNodes; j++){
//                 u_coarse[j] = fMeshVector[0]->NodeVec()[connecCoarse[j]] -> getVelocity(0);
//                 v_coarse[j] = fMeshVector[0]->NodeVec()[connecCoarse[j]] -> getVelocity(1);
//                 p_coarse[j] = fMeshVector[0]->NodeVec()[connecCoarse[j]] -> getPressure();
//             };
            
//             shapeQuad.evaluate(xsi,phi_);
            
//             for (int j=0; j<nElNodes; j++){
//                 u += u_coarse[j] * phi_[j];
//                 v += v_coarse[j] * phi_[j];
//                 p += p_coarse[j] * phi_[j];
//             };
            
//             double wFunc = fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> 
//                 getWeightFunction();
            
//             double u_int = fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> getVelocity(0) * wFunc + u * (1. - wFunc);
//             double v_int = fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> getVelocity(1) * wFunc + v * (1. - wFunc);
//             double p_int = fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> getPressure() * wFunc + p * (1. - wFunc);
            
//             fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> setVelocityArlequin(0,u_int);
//             fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> setVelocityArlequin(1,v_int);
//             fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> setPressureArlequin(p_int);
            
//         };
        
//         for (int i=0; i<fMeshVector[0]->NNodes(); i++){
//             fMeshVector[0]->NodeVec()[i] -> setVelocityArlequin(0,fMeshVector[0]->NodeVec()[i] -> getVelocity(0));
//             fMeshVector[0]->NodeVec()[i] -> setVelocityArlequin(1,fMeshVector[0]->NodeVec()[i] -> getVelocity(1));
//             fMeshVector[0]->NodeVec()[i] -> setPressureArlequin(fMeshVector[0]->NodeVec()[i] -> getPressure());
//         };

//         if (fMeshVector[0]->getProblemParameters().getExactSolution() && fMeshVector[1]->getProblemParameters().getExactSolution()) computeErrorPoisson();
        
//         // Compute and print drag and lift coefficients
//         // if (fMeshVector[1]->getComputeDragAndLift()){
//         //     dragAndLiftCoefficients(dragLift);
//         // };

//         if (rank == 0) {
//             //Printing results
//             printResultsCoarse(iTimeStep);
//             printResultsFine(iTimeStep);
//         };        
//      };
        
    return 0;
};

//------------------------------------------------------------------------------
//----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
//------------------------------------------------------------------------------

int Arlequin::solveArlequinProblemMoving(int iterNumber, double tolerance,
                                                  int problem_type, 
                                                  int time_dependency){

                                                    
//     std::ofstream dragLift;
//     int DIM = fMeshVector[0]->Dimension();
//     int DEG = fMeshVector[0]->GetDefaultOrder();
//     dragLift.open("dragLift.dat", std::ofstream::out | std::ofstream::app);
//     if (rank == 0) {
//         dragLift << "Time   Pressure Drag   Pressure Lift " 
//                  << "Friction Drag  Friction Lift Drag    Lift " 
//                  << std::endl;
//     };   

//     if ((problem_type < 1) || (problem_type > 2)){
//         std::cout << "WRONG PROBLEM TYPE." << std::endl;
//         return 0;
//     };

//     if (time_dependency == 0) fMeshVector[1]->getProblemParameters().GetNTimeSteps() = 1;
    
//     fMeshVector[1]->getProblemParameters().setTimeInstant(0);
    
//     // //Construct the glue zone based on some defined criterion
//     // setCouplingZone();

//     // //Computes the Weight function for all the finite elements
//     // setWeightFunction(16.);
   
//     //Computes the Nodal correspondence between fine nodes and coarse elements
//     // setNodalCorrespondenceFine();

//     double integScheme = fMeshVector[1]->getProblemParameters().GetSpectralRadius();

//     double alpha_f = 1. / (1. + integScheme);
//     double alpha_m = 0.5 * (3. - integScheme) / (1. + integScheme);
//     double gamma = 0.5 + alpha_m - alpha_f;
//     if (rank == 0) std::cout << "Time integ parameters: " << alpha_f << " " << alpha_m << " " << gamma << std::endl;

//     // Computes the system size
//     int sysSize = 3 * fMeshVector[0]->NNodes() + 3 * fMeshVector[1]->NNodes() + 2 * numNodesGlueZoneFine;
    
//     fMeshVector[1]->getProblemParameters().SetNTimeSteps(2000);

//     for (int iTimeStep = 0; iTimeStep < fMeshVector[1]->getProblemParameters().GetNTimeSteps(); iTimeStep++){
        
//         if (rank == 0) {std::cout << "------------------------- TIME STEP = "
//                                   << iTimeStep << " -------------------------"
//                                   << std::endl;}
//         PetscMemoryGetCurrentUsage(&bytes);
//         PetscPrintf(PETSC_COMM_WORLD,"Memory used %g M\n",bytes/(1024*1024));
        
//         //Updates velocity and acceleration
//         for (int i = 0; i < fMeshVector[0]->NNodes(); i++){
//             VecDouble accel(DIM), u(DIM), uprev(DIM);
            
//             //Compute acceleration
//             u[0] = fMeshVector[0]->NodeVec()[i] -> getVelocity(0);
//             u[1] = fMeshVector[0]->NodeVec()[i] -> getVelocity(1);

//             fMeshVector[0]->NodeVec()[i] -> setPreviousVelocity(u);
            
//             accel[0] = fMeshVector[0]->NodeVec()[i] -> getAcceleration(0);
//             accel[1] = fMeshVector[0]->NodeVec()[i] -> getAcceleration(1);
            
//             fMeshVector[0]->NodeVec()[i] -> setPreviousAcceleration(accel);

//             accel[0] *= (gamma - 1.) / gamma;
//             accel[1] *= (gamma - 1.) / gamma;
            
//             fMeshVector[0]->NodeVec()[i] -> setAcceleration(accel);            

//         };

//         // double f = .35;
//         // double w = 2 * pi * f;

//         for (int i = 0; i < fMeshVector[1]->NNodes(); i++){
//             VecDouble accel(DIM), u(DIM), uprev(DIM);
            
//             //Compute acceleration
//             u[0] = fMeshVector[1]->NodeVec()[i] -> getVelocity(0);
//             u[1] = fMeshVector[1]->NodeVec()[i] -> getVelocity(1);

//             fMeshVector[1]->NodeVec()[i] -> setPreviousVelocity(u);
            
//             accel[0] = fMeshVector[1]->NodeVec()[i] -> getAcceleration(0);
//             accel[1] = fMeshVector[1]->NodeVec()[i] -> getAcceleration(1);
            
//             fMeshVector[1]->NodeVec()[i] -> setPreviousAcceleration(accel);

//             accel[0] *= (gamma - 1.) / gamma;
//             accel[1] *= (gamma - 1.) / gamma;
            
//             fMeshVector[1]->NodeVec()[i] -> setAcceleration(accel);



//             VecDouble xn(DIM);
//             VecDouble xi = fMeshVector[1]->NodeVec()[i] -> getInitialCoordinates();       
//             VecDouble x = fMeshVector[1]->NodeVec()[i] -> getCoordinates();       
//             double pi = fMeshVector[0]->getProblemParameters().getPi();
//             double a = -20 * pi / 180 + 10 * pi / 180 * std::cos(2.*pi*iTimeStep*fMeshVector[1]->getProblemParameters().GetTimeStep());// + 10 * pi / 180;

//             // std::cout << " AAA " << a << std::endl;

//             xn[0] = 0.5 + (xi[0]-0.5) * std::cos(a) - (xi[1]-0.0) * std::sin(a);
//             xn[1] = 0.0 + (xi[0]-0.5) * std::sin(a) + (xi[1]-0.0) * std::cos(a);

//             u[0] = (xn[0] - x[0]) / fMeshVector[1]->getProblemParameters().GetTimeStep();
//             u[1] = (xn[1] - x[1]) / fMeshVector[1]->getProblemParameters().GetTimeStep();


//             fMeshVector[1]->NodeVec()[i] -> setMeshVelocity(u);
      
//             fMeshVector[1]->NodeVec()[i] -> setPreviousCoordinates(0,x[0]);
//             fMeshVector[1]->NodeVec()[i] -> setPreviousCoordinates(1,x[1]);

//             fMeshVector[1]->NodeVec()[i] -> setCoordinates(xn);
//         };
        
//         setNodalCorrespondenceFine();
//         setSignaledDistance();
//         setWeightFunction(1.);
        
//         //STARTS NEWTON-RAPHSON
//         for (int inewton = 0; inewton < iterNumber; inewton++){
            
//             std::clock_t t1 = std::clock();
            
//             // Preallocates the matrix
//             ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
//                                 sysSize, sysSize, 400, NULL, 600, NULL, &A); 
            
//             CHKERRQ(ierr);
            
//             // Divides the matrix between the processes
//             ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
            
//             //Create PETSc vectors
//             ierr = VecCreate(PETSC_COMM_WORLD, &b); CHKERRQ(ierr);
//             ierr = VecSetSizes(b, PETSC_DECIDE, sysSize); CHKERRQ(ierr);
//             ierr = VecSetFromOptions(b); CHKERRQ(ierr); 
//             ierr = VecDuplicate(b, &u); CHKERRQ(ierr);
                        
//             assembleArlequinSystem();
            
//             //Assemble matrices and vectors
//             ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//             ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            
//             ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
//             ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
            
//             //ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//             // ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            
//             //Create KSP context to solve the linear system
//             ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
            
//             ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
            



//             // ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE,0, NULL, &nullsp);
//             // ierr = MatSetNullSpace(A, nullsp);
//             // ierr = MatNullSpaceDestroy(&nullsp);




//             // ierr = KSPSetTolerances(ksp,1.e-7,1.e-10,PETSC_DEFAULT,
//             //                         1000);CHKERRQ(ierr);
            
//             // //ierr = KSPGMRESSetRestart(ksp, 10); CHKERRQ(ierr);
            
//             // ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
            
//             // ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);
            
//             // //ierr = KSPSetPCSide(ksp, PC_RIGHT);
//             // ierr = KSPSetType(ksp,KSPLSQR); CHKERRQ(ierr);
            
//             // ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
//             // // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);


// #if defined(PETSC_HAVE_MUMPS)
//         ierr = KSPSetType(ksp,KSPPREONLY);
//         ierr = KSPGetPC(ksp,&pc);
//         ierr = PCSetType(pc, PCLU);

//         ierr = PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
//         PCFactorSetUpMatSolverType(pc);
//         PCFactorGetMatrix(pc,&F);

//         PetscInt ival,icntl;
//         icntl = 14; ival = 60;
//         MatMumpsSetIcntl(F,icntl,ival);
        
// #endif
        
        
//         ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
//         ierr = KSPSetUp(ksp);
        
        
// #if defined(PETSC_HAVE_MUMPS)
//         PetscInt  info1,info2,icntl14;

//         MatMumpsGetInfo(F,1,&info1);
//         MatMumpsGetInfo(F,2,&info2);
//         MatMumpsGetIcntl(F,14,&icntl14);
//         if((rank == 0) && (info1 != 0)) std::cout << " INFO(1) = " << info1
//                                                       << " " << info2 << " " 
//                                                       << icntl14 << std::endl;
// #endif
            
// // #if defined(PETSC_HAVE_MUMPS)
// //             ierr = KSPSetType(ksp,KSPPREONLY);
// //             ierr = KSPGetPC(ksp,&pc);
// //             ierr = PCSetType(pc, PCLU);
// // #endif
            
// //             ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
// //             ierr = KSPSetUp(ksp);
            




//            // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);

//             ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
            
//             ierr = KSPGetTotalIterations(ksp, &iterations);
            
//             //if (rank == 0)std::cout << "GMRES Iterations = " << iterations << std::endl;
            
//             //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            
//             //Gathers the solution vector to the master process
//             ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
            
//             ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
//             CHKERRQ(ierr);
            
//             ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
//             CHKERRQ(ierr);
            
//             ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
            
//             //Updates nodal values
//             double u_[DIM];
//             double normU = 0.;
//             double normP = 0.;
//             double normL = 0.;
//             double normT = 0.;
//             double p_;
//             Ione = 1;

//             for (int i = 0; i < fMeshVector[0]->NNodes(); ++i){
//                 for (int k = 0; k < DIM; k++){
//                     Ii = (DIM+1) * i + k;
//                     ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                     // if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
//                     u_[k] = val;
//                     normU += val*val;
//                     fMeshVector[0]->NodeVec()[i] -> incrementAcceleration(k,u_[k]);
//                     fMeshVector[0]->NodeVec()[i] -> incrementVelocity(k,u_[k]*gamma*fMeshVector[1]->getProblemParameters().GetTimeStep());
//                 }
//                 Ii = (DIM+1) * i + DIM;
//                 ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
//                 p_ = val;
//                 normP += val*val;
//                 fMeshVector[0]->NodeVec()[i] -> incrementPressure(p_);
//             };
//             for (int i = 0; i < fMeshVector[1]->NNodes(); ++i){
//                 for (int k = 0; k < DIM; k++){
//                     Ii = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * i + k;
//                     ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                     u_[k] = val;
//                     normU += val*val;
//                     fMeshVector[1]->NodeVec()[i] -> incrementAcceleration(k,u_[k]);
//                     fMeshVector[1]->NodeVec()[i] -> incrementVelocity(k,u_[k]*gamma*fMeshVector[1]->getProblemParameters().GetTimeStep());
//                 }        
//                 Ii = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * i + DIM;
//                 ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
//                 p_ = val;
//                 normP += val*val;
//                 fMeshVector[1]->NodeVec()[i] -> incrementPressure(p_);
//             };
//             for (int i = 0; i < numNodesGlueZoneFine; ++i){
//                 for (int k = 0; k < DIM; k++){
//                     Ii = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * fMeshVector[1]->NNodes() + DIM * i + k;
//                     ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                     u_[k] = val;
//                     fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> incrementLagrangeMultiplier(k,u_[k]);
//                     normL += val*val;
//                 }
//             };

//             //Computes the solution vector norm
//             ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
            
//             std::clock_t t2 = std::clock();
            
//             if(rank == 0){                
//                 std::cout<<"Iteration = " << inewton << " (" << iterations <<  
//                     ")  Du Norm = " << std::scientific << sqrt(normU) 
//                          << " " << sqrt(normP) 
//                          << " " << sqrt(normL) 
//                          << " " << val << 
//                     "  Time (s) = " << std::fixed << 
//                     1000.*(t2-t1)/CLOCKS_PER_SEC/1000. << std::endl;
//             };
            
//             ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
//             ierr = VecDestroy(&b); CHKERRQ(ierr);
//             ierr = VecDestroy(&u); CHKERRQ(ierr);
//             ierr = VecDestroy(&All); CHKERRQ(ierr);
//             ierr = MatDestroy(&A); CHKERRQ(ierr);
            
//             if(val <= tolerance){
//                 break;            
//             };          

//             //Updates SUPG Parameter
//             // for (int i = 0; i < fMeshVector[1]->NElements(); i++){
//             //     (*fMeshVector[1]->ElementVec()[i]) -> getParameterSUPG();
//             // };   
//             // for (int i = 0; i < fMeshVector[0]->NElements(); i++){
//             //     (*fMeshVector[0]->ElementVec()[i]) -> getParameterSUPG();
//             // };

//         };

//         //Compute real velocity
        
//         ShapeFunction shapeQuad(DIM,DEG);
//         int nElNodes = fMeshVector[1]->NElNodes(); 
//         VecDouble phi_(nElNodes);
        
//         for (int i = 0; i<fMeshVector[1]->NNodes(); i++){
//             fMeshVector[1]->NodeVec()[i] -> setVelocityArlequin(0,fMeshVector[1]->NodeVec()[i] -> getVelocity(0));
//             fMeshVector[1]->NodeVec()[i] -> setVelocityArlequin(1,fMeshVector[1]->NodeVec()[i] -> getVelocity(1));
//             fMeshVector[1]->NodeVec()[i] -> setPressureArlequin(fMeshVector[1]->NodeVec()[i] ->getPressure());
//         };
        
//         for (int i = 0; i<numNodesGlueZoneFine; i++){
//             double u_coarse[nElNodes], v_coarse[nElNodes], p_coarse[nElNodes];
            
//             double u = 0.;
//             double v = 0.;
//             double p = 0.;
            
//             int elCoarse = fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> getNodalElemCorrespondence();
//             VecDouble xsi = fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> getNodalXsiCorrespondence();
            
//             VecInt connecCoarse = fMeshVector[0]->ElementVec()[elCoarse] -> getConnectivity();
            
//             for (int j = 0; j < nElNodes; j++){
//                 u_coarse[j] = fMeshVector[0]->NodeVec()[connecCoarse[j]] -> getVelocity(0);
//                 v_coarse[j] = fMeshVector[0]->NodeVec()[connecCoarse[j]] -> getVelocity(1);
//                 p_coarse[j] = fMeshVector[0]->NodeVec()[connecCoarse[j]] -> getPressure();
//             };
            
//             shapeQuad.evaluate(xsi,phi_);
            
//             for (int j = 0; j < nElNodes; j++){
//                 u += u_coarse[j] * phi_[j];
//                 v += v_coarse[j] * phi_[j];
//                 p += p_coarse[j] * phi_[j];
//             };
            
//             double wFunc = fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> 
//                 getWeightFunction();
            
//             double u_int = fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> getVelocity(0) * wFunc + u * (1. - wFunc);
//             double v_int = fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> getVelocity(1) * wFunc + v * (1. - wFunc);
//             double p_int = fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> getPressure() * wFunc + p * (1. - wFunc);
            
//             fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> setVelocityArlequin(0,u_int);
//             fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> setVelocityArlequin(1,v_int);
//             fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> setPressureArlequin(p_int);
            
//         };
        
//         for (int i=0; i<fMeshVector[0]->NNodes(); i++){
//             fMeshVector[0]->NodeVec()[i] -> setVelocityArlequin(0,fMeshVector[0]->NodeVec()[i] -> getVelocity(0));
//             fMeshVector[0]->NodeVec()[i] -> setVelocityArlequin(1,fMeshVector[0]->NodeVec()[i] -> getVelocity(1));
//             fMeshVector[0]->NodeVec()[i] -> setPressureArlequin(fMeshVector[0]->NodeVec()[i] -> getPressure());
//         };
        
//         // // Compute and print drag and lift coefficients
//         // if (fMeshVector[1]->getComputeDragAndLift()){
//         //     dragAndLiftCoefficients(dragLift);
//         // };

//         if (rank == 0) {
                       
//             //Printing results
//             printResultsCoarse(iTimeStep);
//             printResultsFine(iTimeStep);
//         };
//     };
        
//     return 0;

// };



// //------------------------------------------------------------------------------
// //----------------COMPUTE ARLEQUIN COUPLED NAVIER-STOKES PROBLEM----------------
// //------------------------------------------------------------------------------

// int Arlequin::solveFSIArlequin(int iterNumber, double tolerance,
//                                         int problem_type, int iTimeStep){


//     std::ofstream dragLift;
//     int DIM = fMeshVector[0]->Dimension();
//     int DEG = fMeshVector[0]->GetDefaultOrder();
//     dragLift.open("dragLift.dat", std::ofstream::out | std::ofstream::app);

//     if ((problem_type < 1) || (problem_type > 2)){
//         std::cout << "WRONG PROBLEM TYPE." << std::endl;
//         return 0;
//     };

//     // Computes the system size
//     int sysSize = 3 * fMeshVector[0]->NNodes() + 3 * fMeshVector[1]->NNodes() + 2 * numNodesGlueZoneFine;
    
//     setNodalCorrespondenceFine();
//     setSignaledDistance();
//     setWeightFunction(1.);

//     double &alpha_f = fMeshVector[1]->getProblemParameters().getAlphaF();
//     double &alpha_m = fMeshVector[1]->getProblemParameters().getAlphaM();
//     double &gamma = fMeshVector[1]->getProblemParameters().getGamma();
    
//     //STARTS NEWTON-RAPHSON
//     for (int inewton = 0; inewton < iterNumber; inewton++){
        
//         std::clock_t t1 = std::clock();
        
//         // Preallocates the matrix
//         ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
//                             sysSize, sysSize, 400, NULL, 600, NULL, &A); 
        
//         CHKERRQ(ierr);
        
//         // Divides the matrix between the processes
//         ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
        
//         //Create PETSc vectors
//         ierr = VecCreate(PETSC_COMM_WORLD, &b); CHKERRQ(ierr);
//         ierr = VecSetSizes(b, PETSC_DECIDE, sysSize); CHKERRQ(ierr);
//         ierr = VecSetFromOptions(b); CHKERRQ(ierr); 
//         ierr = VecDuplicate(b, &u); CHKERRQ(ierr);
                            
//         assembleArlequinSystem();
        
//         //Assemble matrices and vectors
//         ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//         ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        
//         ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
//         ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
        
//         //Create KSP context to solve the linear system
//         ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
        
//         ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
        
//         // // ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE,0, NULL, &nullsp);
//         // // ierr = MatSetNullSpace(A, nullsp);
//         // // ierr = MatNullSpaceDestroy(&nullsp);


//         // // if (iTimeStep > 5){

//         //     ierr = KSPSetTolerances(ksp,1.e-7,1.e-10,PETSC_DEFAULT,
//         //                             15);CHKERRQ(ierr);
            
//         //     //ierr = KSPGMRESSetRestart(ksp, 10); CHKERRQ(ierr);
            
//         //     ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
            
//         //     ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
            
//         //     //ierr = KSPSetPCSide(ksp, PC_RIGHT);
//         //     ierr = KSPSetType(ksp,KSPGMRES); CHKERRQ(ierr);
            
//         //     ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
//         //     // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
//         // // }else{
       
// #if defined(PETSC_HAVE_MUMPS)
//         ierr = KSPSetType(ksp,KSPPREONLY);
//         ierr = KSPGetPC(ksp,&pc);
//         ierr = PCSetType(pc, PCLU);
        
//         // ierr = PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
//         // PCFactorSetUpMatSolverType(pc);
//         // PCFactorGetMatrix(pc,&F);
        
//         // PetscInt ival,icntl;
//         // icntl = 14; ival = 80;
//         // MatMumpsSetIcntl(F,icntl,ival);
//         // icntl = 28; ival = 2;
//         // MatMumpsSetIcntl(F,icntl,ival);
//         // icntl = 29; ival = 2;
//         // MatMumpsSetIcntl(F,icntl,ival);
//         // icntl = 16; ival = 0;
//         // MatMumpsSetIcntl(F,icntl,ival);
//         // icntl = 4; ival = 3;
//         // MatMumpsSetIcntl(F,icntl,ival);
//         // icntl = 11; ival = 1;
//         // MatMumpsSetIcntl(F,11,1);

//         //MatMumpsSetIcntl(F,21,0);

        
// #endif
//         ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
//         ierr = KSPSetUp(ksp);
        
// // #if defined(PETSC_HAVE_MUMPS)
// //         PetscInt  info1,info2,icntl14;
        
// //         MatMumpsGetInfo(F,1,&info1);
// //         MatMumpsGetInfo(F,2,&info2);

//         // PetscReal info5,info6,info7,info8,info9,info10,info11;
//         // MatMumpsGetRinfog(F,5,&info5);
//         // MatMumpsGetRinfog(F,6,&info6);
//         // MatMumpsGetRinfog(F,7,&info7);
//         // MatMumpsGetRinfog(F,8,&info8);
//         // MatMumpsGetRinfog(F,9,&info9);
//         // MatMumpsGetRinfog(F,10,&info10);
//         // MatMumpsGetRinfog(F,11,&info11);

//         // PetscInt info21,info32;
//         // MatMumpsGetIcntl(F,32,&info32);
//         // MatMumpsGetIcntl(F,21,&info21);

        
//         // if(rank==0) std::cout << "ICNTL = " << info21 << " " << info32 << std::endl;
      
//         // // if(rank==0) std::cout << "INFOG = " << info5 << " " << info6 << " " << info7 << " " << info8 << " " << info9 << " " << info10 << " " << info11 << std::endl;
//         // MatMumpsGetIcntl(F,14,&icntl14);    
//         // if((rank == 0) && (info1 != 0)) std::cout << " INFO(1) = " << info1
//         //                                           << " " << info2 << " " 
//         //                                           << icntl14 << std::endl;
// // #endif
//         // }
//         // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
        
//         ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
        
//         ierr = KSPGetTotalIterations(ksp, &iterations);
        
//         //if (rank == 0)std::cout << "GMRES Iterations = " << iterations << std::endl;
    
//         //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        
//         //Gathers the solution vector to the master process
//         ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
        
//         ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
//         CHKERRQ(ierr);
        
//         ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
//         CHKERRQ(ierr);
        
//         ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
        
//         //Updates nodal values
//         double u_[DIM];
//         double normU = 0.;
//         double normP = 0.;
//         double normL = 0.;
//         double normT = 0.;
//         double p_;
//         Ione = 1;

//         for (int i = 0; i < fMeshVector[0]->NNodes(); ++i){
//             for (int k = 0; k < DIM; k++){
//                 Ii = (DIM+1) * i + k;
//                 ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                 // if (nodesCoarse_[i] -> getDistFunction() > -1.2) val *= 1000.e0;
//                 u_[k] = val;
//                 normU += val*val;
//                 fMeshVector[0]->NodeVec()[i] -> incrementAcceleration(k,u_[k]);
//                 fMeshVector[0]->NodeVec()[i] -> incrementVelocity(k,u_[k]*gamma*fMeshVector[1]->getProblemParameters().GetTimeStep());
//             }
//             Ii = (DIM+1) * i + DIM;
//             ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
//             p_ = val;
//             normP += val*val;
//             fMeshVector[0]->NodeVec()[i] -> incrementPressure(p_);
//         };
//         for (int i = 0; i < fMeshVector[1]->NNodes(); ++i){
//             for (int k = 0; k < DIM; k++){
//                 Ii = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * i + k;
//                 ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                 u_[k] = val;
//                 normU += val*val;
//                 fMeshVector[1]->NodeVec()[i] -> incrementAcceleration(k,u_[k]);
//                 fMeshVector[1]->NodeVec()[i] -> incrementVelocity(k,u_[k]*gamma*fMeshVector[1]->getProblemParameters().GetTimeStep());
//             }        
//             Ii = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * i + DIM;
//             ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
//             p_ = val;
//             normP += val*val;
//             fMeshVector[1]->NodeVec()[i] -> incrementPressure(p_);
//         };
//         for (int i = 0; i < numNodesGlueZoneFine; ++i){
//             for (int k = 0; k < DIM; k++){
//                 Ii = (DIM+1) * fMeshVector[0]->NNodes() + (DIM+1) * fMeshVector[1]->NNodes() + DIM * i + k;
//                 ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
//                 u_[k] = val;
//                 fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> incrementLagrangeMultiplier(k,u_[k]);
//                 normL += val*val;
//             }
//         };
        
//         //Computes the solution vector norm
//         ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
        
//         std::clock_t t2 = std::clock();
        
//         if(rank == 0){
//             std::cout<<"Iteration = " << inewton << " (" << iterations <<  
//                 ")  Du Norm = " << std::scientific << sqrt(normU) 
//                      << " " << sqrt(normP) 
//                      << " " << sqrt(normL) 
//                      << " " << val << 
//                 "  Time (s) = " << std::fixed << 
//                 1000.*(t2-t1)/CLOCKS_PER_SEC/1000. << std::endl;
//         };
        
//         ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
//         ierr = VecDestroy(&b); CHKERRQ(ierr);
//         ierr = VecDestroy(&u); CHKERRQ(ierr);
//         ierr = VecDestroy(&All); CHKERRQ(ierr);
//         ierr = MatDestroy(&A); CHKERRQ(ierr);
//         // ierr = MatDestroy(&F); CHKERRQ(ierr);
        
        
//         if(val <= tolerance){
//             break;            
//         };          
        
//     };

//     if (rank == 0){
//         double normUUprev, normU;
//         normUUprev = 0.;
//         normU = 0.;
//         for (int i = 0; i < fMeshVector[0]->NNodes(); ++i){
//             for (int k = 0; k < DIM; k++){
//                 double u, uPr, weight;
//                 weight = fMeshVector[0]->NodeVec()[i] -> getWeightFunction();
//                 u = fMeshVector[0]->NodeVec()[i] -> getVelocity(k) * weight;
//                 uPr = fMeshVector[0]->NodeVec()[i] -> getPreviousVelocity(k) * weight;
//                 normUUprev += (u-uPr) * (u-uPr);
//                 normU += u*u;
//             }
//         }
//         for (int i = 0; i < fMeshVector[1]->NNodes(); ++i){
//             for (int k = 0; k < DIM; k++){
//                 double u, uPr, weight;
//                 weight = fMeshVector[1]->NodeVec()[i] -> getWeightFunction();
//                 u = fMeshVector[1]->NodeVec()[i] -> getVelocity(k) * weight;
//                 uPr = fMeshVector[1]->NodeVec()[i] -> getPreviousVelocity(k) * weight;
//                 normUUprev += (u-uPr) * (u-uPr);
//                 normU += u*u;
//             }
//         }
//         std::cout << "NORM U  " << std::scientific <<  sqrt(normUUprev / normU) << std::endl;
//     }


//     // std::cout << "AQUI1 " << rank << std::endl;


//     //Compute real velocity
//     ShapeFunction shapeQuad(DIM,DEG);
//     int nElNodes = fMeshVector[1]->NElNodes(); 
//     VecDouble phi_(nElNodes);
    
//     for (int i = 0; i<fMeshVector[1]->NNodes(); i++){
//         fMeshVector[1]->NodeVec()[i] -> setVelocityArlequin(0,fMeshVector[1]->NodeVec()[i] -> getVelocity(0));
//         fMeshVector[1]->NodeVec()[i] -> setVelocityArlequin(1,fMeshVector[1]->NodeVec()[i] -> getVelocity(1));
//         fMeshVector[1]->NodeVec()[i] -> setPressureArlequin(fMeshVector[1]->NodeVec()[i] ->getPressure());
//     };
    
//     for (int i = 0; i<numNodesGlueZoneFine; i++){
//         double u_coarse[nElNodes], v_coarse[nElNodes], p_coarse[nElNodes];
        
//         double u = 0.;
//         double v = 0.;
//         double p = 0.;
        
//         int elCoarse = fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> getNodalElemCorrespondence();
//         VecDouble xsi = fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> getNodalXsiCorrespondence();
        
//         VecInt connecCoarse = fMeshVector[0]->ElementVec()[elCoarse] -> getConnectivity();
        
//         for (int j=0; j<nElNodes; j++){
//             u_coarse[j] = fMeshVector[0]->NodeVec()[connecCoarse[j]] -> getVelocity(0);
//             v_coarse[j] = fMeshVector[0]->NodeVec()[connecCoarse[j]] -> getVelocity(1);
//             p_coarse[j] = fMeshVector[0]->NodeVec()[connecCoarse[j]] -> getPressure();
//         };
        
//         shapeQuad.evaluate(xsi,phi_);
        
//         for (int j=0; j<nElNodes; j++){
//             u += u_coarse[j] * phi_[j];
//             v += v_coarse[j] * phi_[j];
//             p += p_coarse[j] * phi_[j];
//         };
        
//         double wFunc = fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> 
//             getWeightFunction();
        
//         double u_int = fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> getVelocity(0) * wFunc + u * (1. - wFunc);
//         double v_int = fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> getVelocity(1) * wFunc + v * (1. - wFunc);
//         double p_int = fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> getPressure() * wFunc + p * (1. - wFunc);
        
//         fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> setVelocityArlequin(0,u_int);
//         fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> setVelocityArlequin(1,v_int);
//         fMeshVector[1]->NodeVec()[nodesGlueZoneFine_[i]] -> setPressureArlequin(p_int);
        
//     };
//     // std::cout << "AQUI2 " << rank << std::endl;

//     for (int i=0; i<fMeshVector[0]->NNodes(); i++){
//         fMeshVector[0]->NodeVec()[i] -> setVelocityArlequin(0,fMeshVector[0]->NodeVec()[i] -> getVelocity(0));
//         fMeshVector[0]->NodeVec()[i] -> setVelocityArlequin(1,fMeshVector[0]->NodeVec()[i] -> getVelocity(1));
//         fMeshVector[0]->NodeVec()[i] -> setPressureArlequin(fMeshVector[0]->NodeVec()[i] -> getPressure());
//     };
//     // std::cout << "AQUI3 " << rank << std::endl;
    
//     return 0;

};


void Arlequin::stabilizeArlequin(MatrixDouble &A0, MatrixDouble &A1, 
                                          MatrixDouble &C0, MatrixDouble &C1,
                                          MatrixDouble &E, VecDouble &b0, 
                                          VecDouble &b1, double &tArlq0, double &tArlq1){
    
    //There are in general three main options for taking the norm of a matrix: the L2, L2 and Linfty norms.
    //In eigen they can be simply obtained by:
    // normL2 = mat.norm(); 
    // normL1 = mat.lpNorm<1>(); 
    // normInfty = mat.lpNorm<Infinity>();

    switch (fArlequinStab)
    {
    case ArlequinStabType::ENoStab:
        {
            tArlq0 = 0.;
            tArlq1 = 0.;
            break;
        }
    
    case ArlequinStabType::EOption1:
        {
            double normC0 = C0.norm();
            double normC1 = C1.norm();
            double normA0 = A0.norm();
            double normA1 = A1.norm();
            double normB0 = b0.norm();
            double normB1 = b1.norm();
            double normE = E.norm();
            tArlq0 = std::min({normC0/normA0, normC0/normE, normC0/normB0});
            tArlq1 = std::min({normC1/normA1, normC1/normE, normC1/normB1});
            break;
        }

    case ArlequinStabType::EOption2:
        {
            double normC0 = C0.norm();
            double normC1 = C1.norm();
            double normA0 = A0.norm();
            double normA1 = A1.norm();
            double normB0 = b0.norm();
            double normB1 = b1.norm();
            double normE = E.norm();
            double aux1 = std::min({normC0/normA0, normC0/normE, normC0/normB0, normC1/normA1, normC1/normE, normC1/normB1});
            tArlq0 = aux1;
            tArlq1 = aux1;
            break;
        }
    
    case ArlequinStabType::EOption3:
        {
            double normC0 = C0.norm();
            double normC1 = C1.norm();
            double normA0 = A0.norm();
            double normA1 = A1.norm();
            double normE = E.norm();
            double normB0 = b0.norm();
            double normB1 = b1.norm();
            
            tArlq0 = std::min({normC0/normA0, normC0/normE, normC0/normB0});
            tArlq1 = std::min({normC1/normA1, normC1/normE, normC1/normB1});

            double normC = std::min(normC0,normC1);
            tArlq0 *= normC/normC0;
            tArlq1 *= normC/normC1;
            break;
        }

    case ArlequinStabType::EOption4:
        {
            double normC1 = C1.norm();
            double normA1 = A1.norm();
            double normE = E.norm();
            double normB1 = b1.norm();
            tArlq0 = 0.;
            tArlq1 = std::min({normC1/normA1, normC1/normE, normC1/normB1});
            break;
        }

    case ArlequinStabType::EOption5:
        {
            double normC0 = C0.norm();
            double normC1 = C1.norm();
            double normA1 = A1.norm();
            double normE = E.norm();
            double normB1 = b1.norm();
            double normC = std::min(normC0,normC1);
            tArlq0 = 0.;
            tArlq1 = std::min({normC/normA1, normC/normE, normC/normB1});
            break;
        }
    default:
        PanicButton();
        break;
    }

}


void Arlequin::stabilizeArlequin(std::vector<MatrixDouble> &Stiffness, int64_t &element){
    
    //There are in general three main options for taking the norm of a matrix: the L2, L2 and Linfty norms.
    //In eigen they can be simply obtained by:
    // normL2 = mat.norm(); 
    // normL1 = mat.lpNorm<1>(); 
    // normInfty = mat.lpNorm<Infinity>();
    if (fMeshVector[0]->getProblemParameters().ProbType() != EPoisson){
        PanicButton();
    }

    switch (fArlequinStab)
    {
    case ArlequinStabType::ENoStab:
        {
            break;
        }
    
    case ArlequinStabType::EOption1:
        {
            double normC = Stiffness[0].norm();
            double normE = Stiffness[1].norm();
            double normA = Stiffness[2].norm();
            double tArlq = std::min({normC/normA, normC/normE});
            Stiffness[1] *= tArlq;
            Stiffness[2] *= tArlq;
            // tArlq1 = std::min({normC1/normA1, normC1/normE, normC1/normB1});
            break;
        }

    case ArlequinStabType::EOption2:
        {
            // double normC0 = C0.norm();
            // double normC1 = C1.norm();
            // double normA0 = A0.norm();
            // double normA1 = A1.norm();
            // double normB0 = b0.norm();
            // double normB1 = b1.norm();
            // double normE = E.norm();
            // double aux1 = std::min({normC0/normA0, normC0/normE, normC0/normB0, normC1/normA1, normC1/normE, normC1/normB1});
            // double tArlq0 = aux1;
            // double tArlq1 = aux1;
            break;
        }
    
    case ArlequinStabType::EOption3:
        {
            // double normC0 = C0.norm();
            // double normC1 = C1.norm();
            // double normA0 = A0.norm();
            // double normA1 = A1.norm();
            // double normE = E.norm();
            // double normB0 = b0.norm();
            // double normB1 = b1.norm();
            
            // tArlq0 = std::min({normC0/normA0, normC0/normE, normC0/normB0});
            // tArlq1 = std::min({normC1/normA1, normC1/normE, normC1/normB1});

            // double normC = std::min(normC0,normC1);
            // tArlq0 *= normC/normC0;
            // tArlq1 *= normC/normC1;
            break;
        }

    case ArlequinStabType::EOption4:
        {
            // double normC1 = C1.norm();
            // double normA1 = A1.norm();
            // double normE = E.norm();
            // double normB1 = b1.norm();
            // tArlq0 = 0.;
            // tArlq1 = std::min({normC1/normA1, normC1/normE, normC1/normB1});
            break;
        }

    case ArlequinStabType::EOption5:
        {
            // double normC0 = C0.norm();
            // double normC1 = C1.norm();
            // double normA1 = A1.norm();
            // double normE = E.norm();
            // double normB1 = b1.norm();
            // double normC = std::min(normC0,normC1);
            // tArlq0 = 0.;
            // tArlq1 = std::min({normC/normA1, normC/normE, normC/normB1});
            break;
        }
    default:
        PanicButton();
        break;
    }

}