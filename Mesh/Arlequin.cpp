#include "Arlequin.h"
// #include "hdf5.h"
#include "CouplingLocal.h"
#include "CouplingGlobal.h"
#include <set>
#include "ShapeHexahedron.h"
#include "ShapeOneDLin.h"
#include "ShapeOneDQua.h"
#include "ShapeOneDCub.h"
#include "ShapeQuadrilateralLin.h"
#include "ShapePoint.h"
#include "ShapeTetrahedronLin.h"
#include "ShapeTetrahedronQua.h"
#include "ShapeTetrahedronCub.h"
#include "ShapeTriangleLin.h"
#include "ShapeTriangleQua.h"
#include "ShapeTriangleCub.h"
#include "VTUGenerator.h"

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

Arlequin::Arlequin(std::vector<CompMesh *> &meshvec, double k0, double k1, ArlequinStabType stab){
    fMeshVector = meshvec;
    fMeshVector.resize(3);
    fMeshVector[2] = new CompMesh();
    fMeshVector[2]->SetNStateVariables(1);
    fArlequinStab = stab;
    fK0 = k0;
    fK1 = k1;
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
                if (fMeshVector[0]->ElementVec()[ielcoarse]->Dimension() != fMeshVector[0]->Dimension()) continue;
                CouplingGlobal *cglobal = new CouplingGlobal(fMeshVector[0]->Dimension(),ielcoarse,fMeshVector[0],fK0,fK1);
                cglobal->SetGlobalXsi(fLocalIntPointToGlobalXsi[iel]);
                cglobal->SetGlobalElemCorresp(fLocalIntPointToGlobalElement[iel]);
                fGlobalElToLocalEl[ielcoarse].insert(index-1);
                auto *el = fMeshVector[0]->ElementVec()[ielcoarse]->Clone();
                el->SetMesh(fMeshVector[2]);
                el->SetWeakForm(cglobal);
                el->setConnectivity(fMeshVector[2]->ElementVec()[iel]->getConnectivity());
                fMeshVector[2]->ElementVec().push_back(el);
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
    VecDouble xk(3), Xk(3);

    //Compute element boxes for coarse model
    //Only function for straight elements
    for (int jel = 0; jel < fMeshVector[0]->NElements(); jel++){
        connec = fMeshVector[0]->ElementVec()[jel] -> getConnectivity();
        int ncorner = fMeshVector[0]->ElementVec()[jel]->NCornerNodes();
        VecDouble xi(ncorner), yi(ncorner), zi(ncorner);
        for (int i = 0; i < ncorner; i++){
            xi[i] = fMeshVector[0]->NodeVec()[connec[i]] -> getCoordinateValue(0);        
            yi[i] = fMeshVector[0]->NodeVec()[connec[i]] -> getCoordinateValue(1);        
            zi[i] = fMeshVector[0]->NodeVec()[connec[i]] -> getCoordinateValue(2);        
        }
        xk[0] = xi.minCoeff();
        Xk[0] = xi.maxCoeff();
        xk[1] = yi.minCoeff();
        Xk[1] = yi.maxCoeff();
        xk[2] = zi.minCoeff();
        Xk[2] = zi.maxCoeff();      
        
        fMeshVector[0]->ElementVec()[jel] -> setIntersectionParameters(xk, Xk);
    };

    //Compute element boxes for fine model
    //Only function for straight elements
    for (int jel = 0; jel < fMeshVector[1]->NElements(); jel++){
        connec = fMeshVector[1]->ElementVec()[jel] -> getConnectivity();
        int ncorner = fMeshVector[1]->ElementVec()[jel]->NCornerNodes();
        VecDouble xi(ncorner), yi(ncorner), zi(ncorner);
        for (int i = 0; i < ncorner; i++){
            xi[i] = fMeshVector[1]->NodeVec()[connec[i]] -> getCoordinateValue(0);        
            yi[i] = fMeshVector[1]->NodeVec()[connec[i]] -> getCoordinateValue(1);        
            zi[i] = fMeshVector[1]->NodeVec()[connec[i]] -> getCoordinateValue(2);        
        }
        xk[0] = xi.minCoeff();
        Xk[0] = xi.maxCoeff();
        xk[1] = yi.minCoeff();
        Xk[1] = yi.maxCoeff();
        xk[2] = zi.minCoeff();
        Xk[2] = zi.maxCoeff();      
        
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
    VecDouble xsiCC(3);
    std::pair<VecDouble,VecDouble> XK;

    elCorr = 150000;
    VecDouble xsi(DIM);
    VecDouble x_(DIM);
    VecDouble deltaX(DIM);
    VecDouble deltaXsi(DIM);
    xsi.setZero(); x_.setZero(); deltaX.setZero(); deltaXsi.setZero();
    bool flg = true;
    
    
    xsiCC.fill(1.e10);
    xsiCorr.fill(1.e50);
    xsi.fill(1./3.);
    x_.fill(0.);
    
    Element *elemsearch = nullptr;
    if (cmesh->ElementVec()[elSearch]->Dimension() != cmesh->Dimension()){
        for (int i=0; i<cmesh->NElements(); i++){
            if(cmesh->ElementVec()[i]->Dimension() != cmesh->Dimension()) continue;
            elemsearch = cmesh->ElementVec()[i];
            break;
        }
    } else {
        elemsearch = cmesh->ElementVec()[elSearch];
    }
    VecInt connec = elemsearch -> getConnectivity();
    
    auto &integdata = elemsearch->IntegrationData();
    integdata.fAdimCoord = xsi;
    elemsearch->ComputeJacobian();
    int nElNodes = integdata.fPhi.size();

    for (int i = 0; i < nElNodes; i++){
        VecDouble xint = cmesh->NodeVec()[connec[i]] -> getCoordinates();
        for (int k = 0; k < DIM; k++){
            x_[k] += xint[k] * integdata.fPhi[i];
        }        
    };

    double error = 1.e6;
    int iterations = 0;

    while ((error > 1.e-8) && (iterations < 4)) {
        
        iterations++;
        
        for (int k = 0; k < DIM; k++) deltaX[k] = x[k] - x_[k];
        deltaXsi.setZero();
        
        elemsearch -> ComputeJacobian();

        // for (int i = 0; i < DIM; i++)
        //     for (int j = 0; j < DIM; j++)
        //         deltaXsi[i] += ainv(j,i) * deltaX[j];
        deltaXsi = elemsearch->IntegrationData().fA0Inv*deltaX;    

        xsi += deltaXsi;
        x_.setZero();
        
        integdata.fAdimCoord = xsi;
        // shapeQuad.Shape(xsi,phi_);
        
        for (int i=0; i<nElNodes; i++){
            VecDouble xint = cmesh->NodeVec()[connec[i]] -> getCoordinates();
            for (int k = 0; k < DIM; k++)x_[k] += xint[k] * integdata.fPhi[i];
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
        elCorr = elemsearch->Index();
        // return;
    } else {

        int nEl;
        if (fFirstSearch){
            nEl = cmesh->NElements();
        } else {
            nEl = elemsearch -> getNumberOfNeighborElements(); 
        }   

        for (int jel = 0; jel < nEl; jel++){
            if (cmesh->ElementVec()[jel]->Dimension()!= cmesh->Dimension()) continue;
            if (fFirstSearch){
                connec = cmesh->ElementVec()[jel] -> getConnectivity();
            } else {
                connec = cmesh->ElementVec()[elemsearch -> getNeighborElement(jel)] -> getConnectivity();
            }   

            //get boxes information        
            XK = cmesh->ElementVec()[jel] -> getXIntersectionParameter();

            //Chech if the node is inside the element box
            if ((x[0] < XK.first[0]) || (x[0] > XK.second[0]) ||
                (x[1] < XK.first[1]) || (x[1] > XK.second[1])) continue;
            
            //Compute nodal correspondence
            for (int i = DIM+1; i--; ) xsiCC[i] = 1.e10;
    
            for (int i = DIM; i--; ){
                xsi[i] = 1. / 3.;
                x_[i] = 0.;
            }

            // shapeQuad.Shape(xsi,phi_);
            cmesh->ElementVec()[jel]->IntegrationData().fAdimCoord = xsi;
            cmesh->ElementVec()[jel]->ComputeJacobian();

            for (int i = 0; i < nElNodes; i++){
                VecDouble xint = cmesh->NodeVec()[connec[i]] -> getCoordinates();
                for (int k = DIM; k--; )
                    x_[k] += xint[k] * cmesh->ElementVec()[jel]->IntegrationData().fPhi[i];
            };
            
            double error = 1.e6;
            int iterations = 0;

            while ((error > 1.e-8) && (iterations < 4)) {
                
                iterations++;

                for (int k = DIM; k--; ){
                    deltaX[k] = x[k] - x_[k];
                    deltaXsi[k] = 0.;
                }
                
                cmesh->ElementVec()[jel] -> ComputeJacobian();
                auto ainv = cmesh->ElementVec()[jel]->IntegrationData().fA0Inv;
                for (int i = 0; i < DIM; i++)
                    for (int j = 0; j < DIM; j++)
                        deltaXsi[i] += ainv(i,j) * deltaX[j];
                
                for (int k = DIM; k--; ){
                    xsi[k] += deltaXsi[k];
                    x_[k] = 0.;
                }
            
                cmesh->ElementVec()[jel]->IntegrationData().fAdimCoord = xsi;
                cmesh->ElementVec()[jel] -> ComputeJacobian();
                for (int i=0; i<nElNodes; i++){
                    VecDouble xint = cmesh->NodeVec()[connec[i]] -> getCoordinates();
                    for (int k = DIM; k--; ) x_[k] += xint[k] * cmesh->ElementVec()[jel]->IntegrationData().fPhi[i];
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
        auto *el = fMeshVector[2]->ElementVec()[ielem];
        VecInt connec = fMeshVector[2]->ElementVec()[ielem] -> getConnectivity();
        int nElNodes = connec.size();
        VecDouble x1(nElNodes), x2(nElNodes), x3(nElNodes);
        
        
        for (int i = 0; i < nElNodes; i++){
            VecDouble x = fMeshVector[2]->NodeVec()[connec[i]] -> getCoordinates();
            // VecDouble xp = fMeshVector[2]->NodeVec()[connec[i]] -> getPreviousCoordinates();
            
            x1[i] = x[0];//alpha_f * x[0] + (1. - alpha_f) * xp[0];
            x2[i] = x[1];//alpha_f * x[1] + (1. - alpha_f) * xp[1];
            if (DIM == 3) x3[i] = x[2];//alpha_f * x[2] + (1. - alpha_f) * xp[2];
        };

        // std::cout << "XX1 " << x2 << " " << x22 << " " << x222 << std::endl;

        int numberIntPoints = fMeshVector[2]->ElementVec()[ielem] -> getNumberOfIntegrationPoints();
        for (int i = 0; i < numberIntPoints; i++){

            VecDouble x(DIM);
            
            x[0] = el->InterpolateVariable(x1,i);
            x[1] = el->InterpolateVariable(x2,i);

            int elCorr = 0;
            VecDouble xsiCorr(DIM);

            if (fLocalIntPointToGlobalElement[ielem].size()==0) {
                fLocalIntPointToGlobalElement[ielem].resize(numberIntPoints);
                fLocalIntPointToGlobalElement[ielem].setZero();
            };
            if (fLocalIntPointToGlobalXsi[ielem].rows()==0) {
                fLocalIntPointToGlobalXsi[ielem].resize(numberIntPoints,DIM);
                fLocalIntPointToGlobalXsi[ielem].setZero();
            }
            searchNodeCorrespondence(x,fMeshVector[0],elCorr,xsiCorr,
                                     fLocalIntPointToGlobalElement[ielem][i]);
                            
            fLocalIntPointToGlobalElement[ielem][i] = elCorr;
            fLocalIntPointToGlobalXsi[ielem](i,0)=xsiCorr[0];
            fLocalIntPointToGlobalXsi[ielem](i,1)=xsiCorr[1];
            if (DIM == 3)fLocalIntPointToGlobalXsi[ielem](i,2)=xsiCorr[2];
        };

    }; 
    fFirstSearch = false;

};

//------------------------------------------------------------------------------
//--------------------SETS THE COUPLING ZONE IN COARSE MODEL--------------------
//------------------------------------------------------------------------------
// Given a line with coordinates 'start' and 'end' and the
// coordinates of a point 'pnt' the proc returns the shortest 
// distance from pnt to the line and the coordinates of the 
// nearest point on the line.
// 1  Convert the line segment to a vector ('line_vec').
// 2  Create a vector connecting start to pnt ('pnt_vec').
// 3  Find the length of the line vector ('line_len').
// 4  Convert line_vec to a unit vector ('line_unitvec').
// 5  Scale pnt_vec by line_len ('pnt_vec_scaled').
// 6  Get the dot product of line_unitvec and pnt_vec_scaled ('t').
// 7  Ensure t is in the range 0 to 1.
// 8  Use t to get the nearest location on the line to the end
//    of vector pnt_vec_scaled ('nearest').
// 9  Calculate the distance from nearest to pnt_vec_scaled.
// 10 Translate nearest back to the start/end line. 
// Malcolm Kesson 16 Dec 2012
double Arlequin::ShortestDistance(VecDouble &point, VecDouble &startSeg, VecDouble &endSeg){

    VecDouble line_vec = endSeg - startSeg;
    VecDouble pnt_vec = point - startSeg;
    double line_len = line_vec.norm();
    VecDouble line_unitvec = line_vec/line_len;
    VecDouble pnt_vec_scaled = pnt_vec/line_len;
    double t = line_unitvec.dot(pnt_vec_scaled);
    if (t < 0.0){
        t = 0.0;
    } else if (t > 1.){
        t = 1.;
    }
        
    VecDouble nearest = line_vec * t;
    double dist = (pnt_vec - nearest).norm();
    nearest += startSeg;
    return dist;
};



void Arlequin::setSignaledDistance(){

    int nSegments = fMeshVector[0]->GetDefaultOrder();

    //Fine mesh nodes
    for (int ino = 0; ino < fMeshVector[1]->NNodes(); ino++){
        VecDouble x = fMeshVector[1]->NodeVec()[ino] -> getCoordinates();
        fLocalSignaledDistance[ino] = 1.e15;

        for (int i = 0; i < fMeshVector[1]->NElements(); i++){
            auto *el = fMeshVector[1]->ElementVec()[i];

            if (fGlueMatID.find(el->GetWeakForm()->Id())==fGlueMatID.end()) continue;
            VecInt bconnec = el -> getConnectivity();
            
            // Loop over the 1D element segments
            switch (nSegments){
            case 1:
                {
                    VecDouble start(3),end(3);
                    start = fMeshVector[1]->NodeVec()[bconnec[0]] -> getCoordinates();
                    end = fMeshVector[1]->NodeVec()[bconnec[1]] -> getCoordinates();
                    double dist = ShortestDistance(x,start,end);
                    if (fabs(dist)<fabs(fLocalSignaledDistance[ino])){
                        fLocalSignaledDistance[ino] = dist;
                    }
                }
                break;
            case 2:
                {
                    VecDouble start(3),end(3);
                    //Segment 1
                    start = fMeshVector[1]->NodeVec()[bconnec[0]] -> getCoordinates();
                    end = fMeshVector[1]->NodeVec()[bconnec[2]] -> getCoordinates();
                    double dist1 = ShortestDistance(x,start,end);
                    //Segment 2
                    start = fMeshVector[1]->NodeVec()[bconnec[2]] -> getCoordinates();
                    end = fMeshVector[1]->NodeVec()[bconnec[1]] -> getCoordinates();
                    double dist2 = ShortestDistance(x,start,end);

                    if (fabs(dist1)<fabs(fLocalSignaledDistance[ino])){
                        fLocalSignaledDistance[ino] = dist1;
                    }
                    if (fabs(dist2)<fabs(fLocalSignaledDistance[ino])){
                        fLocalSignaledDistance[ino] = dist2;
                    }
                }
                break;
            case 3:
                {
                    VecDouble start(3),end(3);
                    //Segment 1
                    start = fMeshVector[1]->NodeVec()[bconnec[0]] -> getCoordinates();
                    end = fMeshVector[1]->NodeVec()[bconnec[2]] -> getCoordinates();
                    double dist1 = ShortestDistance(x,start,end);
                    //Segment 2
                    start = fMeshVector[1]->NodeVec()[bconnec[2]] -> getCoordinates();
                    end = fMeshVector[1]->NodeVec()[bconnec[3]] -> getCoordinates();
                    double dist2 = ShortestDistance(x,start,end);
                    //Segment 3
                    start = fMeshVector[1]->NodeVec()[bconnec[3]] -> getCoordinates();
                    end = fMeshVector[1]->NodeVec()[bconnec[1]] -> getCoordinates();
                    double dist3 = ShortestDistance(x,start,end);

                    if (fabs(dist1)<fabs(fLocalSignaledDistance[ino])){
                        fLocalSignaledDistance[ino] = dist1;
                    }
                    if (fabs(dist2)<fabs(fLocalSignaledDistance[ino])){
                        fLocalSignaledDistance[ino] = dist2;
                    }
                    if (fabs(dist3)<fabs(fLocalSignaledDistance[ino])){
                        fLocalSignaledDistance[ino] = dist3;
                    }
                }
                break;
            
            default:
                PanicButton();
                break;
            }
        }
    }

    //Coarse mesh nodes
    for (int ino = 0; ino < fMeshVector[0]->NNodes(); ino++){
        VecDouble x = fMeshVector[0]->NodeVec()[ino] -> getCoordinates();
        fGlobalSignaledDistance[ino] = 1.e15;

        for (int i = 0; i < fMeshVector[1]->NElements(); i++){
            auto *el = fMeshVector[1]->ElementVec()[i];

            if (fGlueMatID.find(el->GetWeakForm()->Id())==fGlueMatID.end()) continue;
            VecInt bconnec = el -> getConnectivity();
            
            // Loop over the 1D element segments
            switch (nSegments){
            case 1:
                {
                    VecDouble start(3),end(3);
                    start = fMeshVector[1]->NodeVec()[bconnec[0]] -> getCoordinates();
                    end = fMeshVector[1]->NodeVec()[bconnec[1]] -> getCoordinates();
                    double dist = ShortestDistance(x,start,end);
                    if (fabs(dist)<fabs(fGlobalSignaledDistance[ino])){
                        fGlobalSignaledDistance[ino] = dist;
                    }
                }
                break;
            case 2:
                {
                    VecDouble start(3),end(3);
                    //Segment 1
                    start = fMeshVector[1]->NodeVec()[bconnec[0]] -> getCoordinates();
                    end = fMeshVector[1]->NodeVec()[bconnec[2]] -> getCoordinates();
                    double dist1 = ShortestDistance(x,start,end);
                    //Segment 2
                    start = fMeshVector[1]->NodeVec()[bconnec[2]] -> getCoordinates();
                    end = fMeshVector[1]->NodeVec()[bconnec[3]] -> getCoordinates();
                    double dist2 = ShortestDistance(x,start,end);
                    //Segment 3
                    start = fMeshVector[1]->NodeVec()[bconnec[3]] -> getCoordinates();
                    end = fMeshVector[1]->NodeVec()[bconnec[1]] -> getCoordinates();
                    double dist3 = ShortestDistance(x,start,end);
                    
                    if (fabs(dist1)<fabs(fGlobalSignaledDistance[ino])){
                        fGlobalSignaledDistance[ino] = dist1;
                    }
                    if (fabs(dist2)<fabs(fGlobalSignaledDistance[ino])){
                        fGlobalSignaledDistance[ino] = dist2;
                    }
                    if (fabs(dist3)<fabs(fGlobalSignaledDistance[ino])){
                        fGlobalSignaledDistance[ino] = dist3;
                    }
                }
                break;
            
            default:
                PanicButton();
                break;
            }
        }
    }

    for (int jel = 0; jel < fMeshVector[0]->NElements(); jel++){
        if(fMeshVector[0]->ElementVec()[jel]->Dimension()!=fMeshVector[0]->Dimension())continue;
        auto connect = fMeshVector[0]->ElementVec()[jel]->getConnectivity();
        VecDouble distfunction(connect.size());
        for (int i = 0; i < connect.size(); i++){
            distfunction[i] = fGlobalSignaledDistance[connect[i]];
        }
        fMeshVector[0]->ElementVec()[jel] -> ComputeIntPointDistFunction(distfunction);        
    };
    for (int jel = 0; jel < fMeshVector[1]->NElements(); jel++){
        if(fMeshVector[1]->ElementVec()[jel]->Dimension()!=fMeshVector[1]->Dimension())continue;
        auto connect = fMeshVector[1]->ElementVec()[jel]->getConnectivity();
        VecDouble distfunction(connect.size());
        for (int i = 0; i < connect.size(); i++){
            distfunction[i] = fLocalSignaledDistance[connect[i]];
        }
        fMeshVector[1]->ElementVec()[jel] -> ComputeIntPointDistFunction(distfunction);        
    };
    // // for (int jel = 0; jel < fMeshVector[1]->NElements(); jel++){
    // //     fMeshVector[1]->ElementVec()[jel] -> ComputeIntPointDistFunction();        
    // // };
    
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
        if (fMeshVector[1]->ElementVec()[jel]->Dimension() != fMeshVector[1]->Dimension()) continue; 
        VecInt connec = fMeshVector[1]->ElementVec()[jel] -> getConnectivity();
        flag = 0;
        int nElNodes = connec.size();
        for (int ino = 0; ino < nElNodes; ino++){
            VecDouble x = fMeshVector[1]->NodeVec()[connec[ino]] -> getCoordinates();
            double dist = fLocalSignaledDistance[connec[ino]];
            //  std::cout << "DIST " << dist << std::endl;
            if (dist <= fGlueZoneThickness + 0.001) flag++;
            
        };

        if (flag == nElNodes) {
            CouplingLocal *clocal = new CouplingLocal(fMeshVector[0]->Dimension(),jel,fMeshVector[1],fK0,fK1);
            auto *el = fMeshVector[1]->ElementVec()[jel]->Clone();
            el->SetMesh(fMeshVector[2]);
            el->SetWeakForm(clocal);
            fMeshVector[2]->ElementVec().push_back(el);
            fGluingElementIndex.insert(el->Index()); 
        };        
    };
    int nElNodes = 0;
    //Defines which nodes are in the glue zone
    numElemGlueZoneFine = fMeshVector[2]->NElements();
    for (int i = 0; i < numElemGlueZoneFine; i++){
        auto elcoup = fMeshVector[2]->ElementVec()[i];
        CouplingLocal * coupling = dynamic_cast<CouplingLocal*> (elcoup->GetWeakForm());
        auto locindex = coupling->GetLocalIndex();
        VecInt connec = fMeshVector[1]->ElementVec()[locindex] -> getConnectivity();
        nElNodes = connec.size();
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


    std::cout << "GLUE ZONE - Number of Nodes = " 
              << fMeshVector[2]->NNodes() 
              << " - Number of Elements = " 
              << fMeshVector[2]->NElements() << std::endl;
    

    for (int i = 0; i < fMeshVector[2]->NElements(); i++){
        auto *el =fMeshVector[2]->ElementVec()[i];
        auto connec1 = el->getConnectivity();
        
        VecInt connecAux(nElNodes);
        CouplingLocal * coupling = dynamic_cast<CouplingLocal*> (el->GetWeakForm());
        auto locindex = coupling->GetLocalIndex();

        VecInt connec = fMeshVector[1]->ElementVec()[locindex] -> getConnectivity();
        
        // for (int ino = 0; ino < numNodesGlueZoneFine; ino++)
        for (int k = 0; k < nElNodes; k++) connecAux[k] = FineToGluing[connec[k]];
            // if (nodesGlueZoneFine_[ino] == connec[k]) connecAux[k] = ino;

        fMeshVector[2]->ElementVec()[i] -> setConnectivity(connecAux);
    };

     
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
        if(fMeshVector[0]->ElementVec()[jel]->Dimension() != fMeshVector[0]->Dimension()) continue;
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
        if(fMeshVector[1]->ElementVec()[jel]->Dimension() != fMeshVector[1]->Dimension()) continue;
        for (int i = 0; i < fMeshVector[1]->ElementVec()[jel]->getNumberOfIntegrationPoints(); i++){
            double dist = fMeshVector[1]->ElementVec()[jel]->GetIntPointDistFunction(i);
            double wFuncValue = LocalWeightFunction(dist);
            fMeshVector[1]->ElementVec()[jel] -> setIntegPointWeightFunction(i,wFuncValue);    
        }
            
    };

    VTUGenerator::PrintResults(this,"ArlequinMesh");

    return;
};
