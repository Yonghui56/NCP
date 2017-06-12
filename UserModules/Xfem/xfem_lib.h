
#pragma once

#include <vector>
#include <cmath>

#include "logog.hpp"
#include "MathLib/DataType.h"

//namespace xfem
//{

void exactSol_Mode1(
        const double xx, const double yy,
        const double k1, const double kappa, const double mu, const double /*lambda*/,
        double &uu, double &vv)
{
    const static double pi = 4*atan(1.0);

    double rr = sqrt(xx*xx+yy*yy); // Crack tip is at (0, 0).
//    double drrdx = xx/rr;
//    double drrdy = yy/rr;
//
//    double dthdx = -yy/(xx*xx+yy*yy);
//    double dthdy = 0;
//    if (xx == 0.)
//        dthdy = 0;
//    else
//        dthdy = 1/(xx+yy*yy/xx);

    double th = 0;
    if (xx == 0) {
        if (yy == 0) {
            INFO("Theta is undetermined for (x,y) = (0,0)!");
            th = 0;
        } else if (yy > 0) {
            th = 0.5*pi;
        } else {
            th = 1.5*pi;
        }
    } else if (yy == 0) {
        if (xx > 0)
            th = 0;
        else
            th = pi;
    } else if (xx > 0) {
        th = atan(yy/xx);
    } else if (xx < 0) {
        th = pi+atan(yy/xx);
    } else {
        INFO("Internal error!");
    }

    if (th > pi) // This is important due to the multiplication with 1/2 below!!!
        th = th-2*pi;

    // Get exact solution.
    uu = k1/(2*mu)*sqrt(rr/(2*pi))*cos(0.5*th)*(kappa-1+2*sin(0.5*th)*sin(0.5*th));
    vv = k1/(2*mu)*sqrt(rr/(2*pi))*sin(0.5*th)*(kappa+1-2*cos(0.5*th)*cos(0.5*th));

//    duudrr = 1/8*k1*2^(1/2)*cos(1/2*th)*(kappa+1-2*cos(1/2*th)^2)/mu/pi^(1/2)/rr^(1/2);
//    duudth = -1/8*k1*2^(1/2)/pi^(1/2)*rr^(1/2)*sin(1/2*th)*(kappa+1-6*cos(1/2*th)^2)/mu;
//    duudx = duudrr * drrdx + duudth * dthdx;
//    duudy = duudrr * drrdy + duudth * dthdy;
//
//    dvvdrr = 1/8*k1/mu*2^(1/2)/pi^(1/2)/rr^(1/2)*sin(1/2*th)*(kappa+1-2*cos(1/2*th)^2);
//    dvvdth = 1/8*k1*2^(1/2)/pi^(1/2)*rr^(1/2)*cos(1/2*th)*(kappa+5-6*cos(1/2*th)^2)/mu;
//    dvvdx = dvvdrr * drrdx + dvvdth * dthdx;
//    dvvdy = dvvdrr * drrdy + dvvdth * dthdy;
//
//    Sigma11 = k1/sqrt(2*pi*rr)*cos(0.5*th)*(1-sin(0.5*th)*sin(1.5*th));
//    Sigma22 = k1/sqrt(2*pi*rr)*cos(0.5*th)*(1+sin(0.5*th)*sin(1.5*th));
//    Sigma12 = k1/sqrt(2*pi*rr)*cos(0.5*th)*(sin(0.5*th)*cos(1.5*th));
//
//    % Sigma = [Sigma11 Sigma12; Sigma12 Sigma22];
//    % Eps = -0.25*lambda/(mu*(lambda+mu))*trace(Sigma)*[1 0; 0 1] + 1/(2*mu)*Sigma;
//    % Eps11 = Eps(1,1);
//    % Eps12 = Eps(1,2);
//    % Eps22 = Eps(2,2);
//
//    Eps11 = duudx;
//    Eps12 = 0.5 * (duudy + dvvdx);
//    Eps22 = dvvdy;
}

void searchMinMaxNodes(const MeshLib::IMesh &msh, std::vector<size_t> &found_nodes)
{
    double x_min=1e+99, x_max = -1e+99;
    double y_min=1e+99, y_max = -1e+99;
    //double pt[3];
    //search x min/max
    for (size_t i=0; i<msh.getNumberOfNodes(); i++) {
      const GeoLib::Point *pt = msh.getNodeCoordinatesRef(i);
      if ((*pt)[0]<x_min) x_min = (*pt)[0];
      if ((*pt)[0]>x_max) x_max = (*pt)[0];
      if ((*pt)[1]<y_min) y_min = (*pt)[1];
      if ((*pt)[1]>y_max) y_max = (*pt)[1];
    }

    //search nodes on min/max
    for (size_t i=0; i<msh.getNumberOfNodes(); i++) {
        const GeoLib::Point *pt = msh.getNodeCoordinatesRef(i);
      if (abs((*pt)[0]-x_min)<std::numeric_limits<double>::epsilon()) {
          found_nodes.push_back(i);
      } else if (abs((*pt)[0]-x_max)<std::numeric_limits<double>::epsilon()) {
          found_nodes.push_back(i);
      } else if (abs((*pt)[1]-y_min)<std::numeric_limits<double>::epsilon()) {
          found_nodes.push_back(i);
      } else if (abs((*pt)[1]-y_max)<std::numeric_limits<double>::epsilon()) {
          found_nodes.push_back(i);
      }
    }
}

void getEnrichedNodesElems(const MeshLib::IMesh &msh, const std::vector<double> &ff, std::vector<size_t> &ElemsEnriched, std::vector<size_t> &NodesEnriched)
{
    // Get cut elements and correspoding nodes.

    for (size_t i=0; i<msh.getNumberOfElements(); i++) {
        MeshLib::IElement* e = msh.getElement(i);
        bool hasNegativeNode = false;
        bool hasNonNegativeNode = false;
        size_t cnt = 0;
        for (size_t j=0; j<e->getNumberOfNodes(); j++) {
            if (ff[e->getNodeID(j)] < .0) hasNegativeNode = true;
            else hasNonNegativeNode = true;
            if (msh.getNodeCoordinatesRef(e->getNodeID(j))->getData()[0]<.0) {
                cnt++;
            }
        }
        if (hasNegativeNode && hasNonNegativeNode && cnt==e->getNumberOfNodes())
            ElemsEnriched.push_back(i);
    }

    std::set<size_t> set_nodes;
    for (size_t i=0; i<ElemsEnriched.size(); i++) {
        MeshLib::IElement* e = msh.getElement(ElemsEnriched[i]);
        for (size_t j=0; j<e->getNumberOfNodes(); j++) {
            set_nodes.insert(e->getNodeID(j));
        }
    }
    NodesEnriched.assign(set_nodes.begin(), set_nodes.end());
}

void ShapeFctsStrdFEM(
            const MathLib::LocalVector &xxElem, const MathLib::LocalVector &yyElem,
            const MathLib::LocalVector &xxIntRef, const MathLib::LocalVector &yyIntRef, const MathLib::LocalVector &wwIntRef,
            size_t nQ,
            MathLib::LocalMatrix &N, MathLib::LocalMatrix &dNdx, MathLib::LocalMatrix &dNdy,
            MathLib::LocalVector &xxInt, MathLib::LocalVector &yyInt, MathLib::LocalVector &wwInt
            )
{
    // Compute standard bi-linear shape functions at integration points in
    // the real element.
    const size_t n_ele_nodes = 4;

    // Initialization.
    N     = MathLib::LocalMatrix::Zero(n_ele_nodes, nQ);
    dNdx  = MathLib::LocalMatrix::Zero(n_ele_nodes, nQ);
    dNdy  = MathLib::LocalMatrix::Zero(n_ele_nodes, nQ);

    xxInt = MathLib::LocalVector::Zero(nQ);
    yyInt = MathLib::LocalVector::Zero(nQ);
    wwInt = MathLib::LocalVector::Zero(nQ);

    // Define N, dNdr, dNds at the integration points in the REFERENCE element.
    Eigen::ArrayXd N1 = 0.25*(1-xxIntRef.array())*(1-yyIntRef.array());
    Eigen::ArrayXd N2 = 0.25*(1+xxIntRef.array())*(1-yyIntRef.array());
    Eigen::ArrayXd N3 = 0.25*(1+xxIntRef.array())*(1+yyIntRef.array());
    Eigen::ArrayXd N4 = 0.25*(1-xxIntRef.array())*(1+yyIntRef.array());
    Eigen::ArrayXd dN1dr = -0.25*(1-yyIntRef.array());
    Eigen::ArrayXd dN2dr = +0.25*(1-yyIntRef.array());
    Eigen::ArrayXd dN3dr = +0.25*(1+yyIntRef.array());
    Eigen::ArrayXd dN4dr = -0.25*(1+yyIntRef.array());
    Eigen::ArrayXd dN1ds = -0.25*(1-xxIntRef.array());
    Eigen::ArrayXd dN2ds = -0.25*(1+xxIntRef.array());
    Eigen::ArrayXd dN3ds = +0.25*(1+xxIntRef.array());
    Eigen::ArrayXd dN4ds = +0.25*(1-xxIntRef.array());

    Eigen::MatrixXd NN(n_ele_nodes, nQ); //MatrixXd is column-base
    NN  << N1.transpose(), N2.transpose(), N3.transpose(),  N4.transpose();
    Eigen::MatrixXd dNNdr(n_ele_nodes, nQ);
    dNNdr << dN1dr.transpose(), dN2dr.transpose(), dN3dr.transpose(), dN4dr.transpose();
    Eigen::MatrixXd dNNds(n_ele_nodes, nQ);
    dNNds << dN1ds.transpose(), dN2ds.transpose(), dN3ds.transpose(), dN4ds.transpose();
//    std::cout << "N1=" << N1;
//    std::cout << "dN1dr=" << dN1dr;
//    std::cout << "dN1ds=" << dN1ds;
//    std::cout << "NN=" << NN;
//    std::cout << "dNNdr=" << dNNdr;
//    std::cout << "dNNds=" << dNNds;
//    std::flush(std::cout);

    // Define N, dNdx, dNdy at the integration points in the REAL element.
    N = NN;

    for (size_t i=0; i<nQ; i++) {
        double dxdr = dNNdr.col(i).transpose()*xxElem; // J11
        double dydr = dNNdr.col(i).transpose()*yyElem; // J12
        double dxds = dNNds.col(i).transpose()*xxElem; // J21
        double dyds = dNNds.col(i).transpose()*yyElem; // J22

        // det(J) = J11*J22 - J21*J12
        double detJ = dxdr*dyds - dxds*dydr;
        assert(detJ>.0);
        wwInt(i) = wwIntRef(i) * detJ;

        dNdx.col(i) = ( dyds*dNNdr.col(i) - dydr*dNNds.col(i)) / detJ;
        dNdy.col(i) = (-dxds*dNNdr.col(i) + dxdr*dNNds.col(i)) / detJ;
    }

    // Map integration points from reference to real element.
    xxInt = NN.transpose() * xxElem;
    yyInt = NN.transpose() * yyElem;

}

void IntPoints2DRealElemQuad(
            const MathLib::LocalVector &xxElem, const MathLib::LocalVector &yyElem,
            size_t nQ,
            const MathLib::LocalVector &xxIntRef, const MathLib::LocalVector &yyIntRef, const MathLib::LocalVector &wwIntRef,
            /* MathLib::LocalMatrix &N, MathLib::LocalMatrix &dNdx, MathLib::LocalMatrix &dNdy,*/
            MathLib::LocalVector &xxInt, MathLib::LocalVector &yyInt, MathLib::LocalVector &wwInt
            )
{
    MathLib::LocalMatrix N, dNdx, dNdy;
    ShapeFctsStrdFEM(xxElem, yyElem, xxIntRef, yyIntRef, wwIntRef, nQ, N, dNdx, dNdy, xxInt, yyInt, wwInt);
}

void IntPoints2DLevelSet(const MathLib::LocalVector &ff, const std::vector<GeoLib::Point> &vec_int_ref_xx, const std::vector<double> &vec_int_ref_w, const size_t nQnQ, MathLib::LocalVector &xxIntTotal, MathLib::LocalVector &yyIntTotal, MathLib::LocalVector &wwIntTotal)
{
    const size_t n_ele_nodes = ff.rows();
    std::vector<double> SignVect(n_ele_nodes);
    for (size_t j=0; j<n_ele_nodes; j++)
        SignVect[j] = (ff[j] >= .0 ? 1 : -1);

    if (SignVect[0]==SignVect[1] && SignVect[0]==SignVect[2] && SignVect[0]==SignVect[3]) {
        // Treatment of Type-1-elements (no disc., i.e. all level set values have the same sign).
        xxIntTotal.resize(vec_int_ref_xx.size());
        yyIntTotal.resize(vec_int_ref_xx.size());
        wwIntTotal.resize(vec_int_ref_xx.size());
        for (size_t i=0; i<vec_int_ref_xx.size(); i++) {
            xxIntTotal(i) = vec_int_ref_xx[i][0];
            yyIntTotal(i) = vec_int_ref_xx[i][1];
            wwIntTotal(i) = vec_int_ref_w[i];
        }
    } else if (SignVect[0]*SignVect[1]*SignVect[2]*SignVect[3]<.0) {
        // Treatment of Type-2-elements (one level function value is on the other side than the other three).

        // Find the one node on the other side.
        size_t Pos;
        double xA, xB, xC, xD, xE, xF;
        double yA, yB, yC, yD, yE, yF;
        if (SignVect[0]!=SignVect[1] && SignVect[0]!=SignVect[2] && SignVect[0]!=SignVect[3]) {
            Pos = 1;
            xA = -1+2*ff[3]/(ff[3]-ff[0]);
            yB = -1+2*ff[0]/(ff[0]-ff[1]);
            //IntSide = -SignVect[0] * [ones(nQnQ, 1); ones(nQnQ, 1);-ones(nQnQ, 1);-ones(nQnQ, 1);-ones(nQnQ, 1)];
        } else if (SignVect[1]!=SignVect[0] && SignVect[1]!=SignVect[2] && SignVect[1]!=SignVect[3]) {
            Pos = 2;
            xA = -1+2*ff[0]/(ff[0]-ff[1]);
            yB = -1+2*ff[1]/(ff[1]-ff[2]);
            //IntSide = -SignVect[1] * [ones(nQnQ, 1); ones(nQnQ, 1);-ones(nQnQ, 1);-ones(nQnQ, 1);-ones(nQnQ, 1)];
        } else if (SignVect[2]!=SignVect[0] && SignVect[2]!=SignVect[1] && SignVect[2]!=SignVect[3]) {
            Pos = 3;
            xA = -1+2*ff[1]/(ff[1]-ff[2]);
            yB = -1+2*ff[2]/(ff[2]-ff[3]);
            //IntSide = -SignVect[2] * [ones(nQnQ, 1); ones(nQnQ, 1);-ones(nQnQ, 1);-ones(nQnQ, 1);-ones(nQnQ, 1)];
        } else if (SignVect[3]!=SignVect[0] && SignVect[3]!=SignVect[1] && SignVect[3]!=SignVect[2]) {
            Pos = 4;
            xA = -1+2*ff[2]/(ff[2]-ff[3]);
            yB = -1+2*ff[3]/(ff[3]-ff[0]);
            //IntSide = -SignVect[3] * [ones(nQnQ, 1); ones(nQnQ, 1);-ones(nQnQ, 1);-ones(nQnQ, 1);-ones(nQnQ, 1)];
        } else {
            ERR("Internal error.");
        }

        double x1, x2, x3, x4;
        double y1, y2, y3, y4;
        x1 =-1; y1 =-1;
        x2 = 1; y2 =-1;
        x3 = 1; y3 = 1;
        x4 =-1; y4 = 1;

        // Make sub-elements and int. points as if the SECOND node were on the other side (Pos=2).
        yA = -1;
        xB = 1;
        xC = 0.5*(xA+xB);
        yC = 0.5*(yA+yB);
        xD = 1/3+2/3*xC;
        yD = -1/3+2/3*yC;
        xE = 0.5*(1+xA);
        yE = -1;
        xF = 1;
        yF = 0.5*(-1+yB);

        // Project integration points into sub-elements.
//        MeshLib::Quadrirateral q1;
//        FemLib::FemIntegrationGaussQuad int_q;
//        int_q.initialize(q1, 2);
//        int_q.getNumberOfSamplingPoints()
        MathLib::LocalVector xxIntRef(vec_int_ref_xx.size()), yyIntRef(vec_int_ref_xx.size()), wwIntRef(vec_int_ref_xx.size());
        for (size_t i=0; i<vec_int_ref_xx.size(); i++) {
            xxIntRef(i) = vec_int_ref_xx[i].getData()[0];
            yyIntRef(i) = vec_int_ref_xx[i].getData()[1];
            wwIntRef(i) = vec_int_ref_w[i];
        }
        MathLib::LocalVector xxElem1(4), yyElem1(4), xxInt1, yyInt1, wwInt1;
        xxElem1 << x1, xA, xC, x4;
        yyElem1 << y1, yA, yC, y4;
        MathLib::LocalVector xxElem2(4), yyElem2(4), xxInt2, yyInt2, wwInt2;
        xxElem2 << x4, xC, xB, x3;
        yyElem2 << y4, yC, yB, y3;
        MathLib::LocalVector xxElem3(4), yyElem3(4), xxInt3, yyInt3, wwInt3;
        xxElem3 << xA, xE, xD, xC;
        yyElem3 << yA, yE, yD, yC;
        MathLib::LocalVector xxElem4(4), yyElem4(4), xxInt4, yyInt4, wwInt4;
        xxElem4 << xB, xC, xD, xF;
        yyElem4 << yB, yC, yD, yF;
        MathLib::LocalVector xxElem5(4), yyElem5(4), xxInt5, yyInt5, wwInt5;
        xxElem5 << x2, xF, xD, xE;
        yyElem5 << y2, yF, yD, yE;
        IntPoints2DRealElemQuad(xxElem1, yyElem1, nQnQ, xxIntRef, yyIntRef, wwIntRef, xxInt1, yyInt1, wwInt1);
        IntPoints2DRealElemQuad(xxElem2, yyElem2, nQnQ, xxIntRef, yyIntRef, wwIntRef, xxInt2, yyInt2, wwInt2);
        IntPoints2DRealElemQuad(xxElem3, yyElem3, nQnQ, xxIntRef, yyIntRef, wwIntRef, xxInt3, yyInt3, wwInt3);
        IntPoints2DRealElemQuad(xxElem4, yyElem4, nQnQ, xxIntRef, yyIntRef, wwIntRef, xxInt4, yyInt4, wwInt4);
        IntPoints2DRealElemQuad(xxElem5, yyElem5, nQnQ, xxIntRef, yyIntRef, wwIntRef, xxInt5, yyInt5, wwInt5);

        MathLib::LocalVector xxEleTotal(xxInt1.rows()*5);
        MathLib::LocalVector yyEleTotal(yyInt1.rows()*5);
        MathLib::LocalVector wwEleTotal(wwInt1.rows()*5);
        xxEleTotal << xxInt1, xxInt2, xxInt3, xxInt4, xxInt5;
        yyEleTotal << yyInt1, yyInt2, yyInt3, yyInt4, yyInt5;
        wwEleTotal << wwInt1, wwInt2, wwInt3, wwInt4, wwInt5;


//    %     patch([x1 xA xC x4], [y1 yA yC y4], 'b')
//    %     patch([x4 xC xB x3], [y4 yC yB y3], 'b')
//    %     patch([xA xE xD xC], [yA yE yD yC], 'b')
//    %     patch([xB xC xD xF], [yB yC yD yF], 'b')
//    %     patch([x2 xF xD xE], [y2 yF yD yE], 'b')
//    %
//    %     plot(xxInt1, yyInt1, 'k*')
//    %     plot(xxInt2, yyInt2, 'k*')
//    %     plot(xxInt3, yyInt3, 'k*')
//    %     plot(xxInt4, yyInt4, 'k*')
//    %     plot(xxInt5, yyInt5, 'k*')

        // Rotate integration points according to which of the nodes is really on the other side.
        MathLib::LocalMatrix T(2,2);
        if (Pos == 1)
            T << 0, 1, -1, 0; // Rotation -90 degree.
        else if (Pos == 2)
            T << 1, 0, 0, 1; // Nothing to do.
        else if (Pos == 3)
            T << 0, -1, 1, 0; // Rotation 90 degree.
        else if (Pos == 4)
            T << -1, 0, 0, -1; // Rotation 180 degree.

        MathLib::LocalMatrix xxyyIntTotal(2, xxEleTotal.rows());
        xxyyIntTotal << xxEleTotal.transpose(), yyEleTotal.transpose();
        MathLib::LocalMatrix HelpMat = T * xxyyIntTotal;
        xxIntTotal = HelpMat.row(0).transpose();
        yyIntTotal = HelpMat.row(1).transpose();
        wwIntTotal = wwEleTotal;
    } else {
        // Treatment of Type-3-elements (two level function values are on the other side than the other two).
        size_t Pos;
        double xA, xB;
        double yA, yB;
        if (SignVect[0]==SignVect[1] && SignVect[2]==SignVect[3]) {
            Pos = 1;
            xA = -1+2*ff[3]/(ff[3]-ff[0]);
            xB = -1+2*ff[2]/(ff[2]-ff[1]);
            //IntSide = SignVect[0] * [-ones(nQnQ, 1); ones(nQnQ, 1)];
        } else if (SignVect[0]==SignVect[3] && SignVect[1]==SignVect[2]) {
            Pos = 2;
            xA = -1+2*ff[0]/(ff[0]-ff[1]);
            xB = -1+2*ff[3]/(ff[3]-ff[2]);
            //IntSide = SignVect[0] * [ones(nQnQ, 1); -ones(nQnQ, 1)];
        } else if (SignVect[0]==SignVect[2] && SignVect[1]==SignVect[3]) {
            ERR("Special situation, disc. can not be uniquely determined from level set.");
        } else {
            ERR("Internal error.");
        }

        double x1, x2, x3, x4;
        double y1, y2, y3, y4;
        x1 =-1; y1 =-1;
        x2 = 1; y2 =-1;
        x3 = 1; y3 = 1;
        x4 =-1; y4 = 1;

        // Make sub-elements and int. points as if the SECOND node were on the other side (Pos=2).
        yA =-1; yB = 1;

        // Project integration points into sub-elements.
        MathLib::LocalVector xxIntRef(vec_int_ref_xx.size()), yyIntRef(vec_int_ref_xx.size()), wwIntRef(vec_int_ref_xx.size());
        for (size_t i=0; i<vec_int_ref_xx.size(); i++) {
            xxIntRef(i) = vec_int_ref_xx[i].getData()[0];
            yyIntRef(i) = vec_int_ref_xx[i].getData()[1];
            wwIntRef(i) = vec_int_ref_w[i];
        }
        MathLib::LocalVector xxElem1(4), yyElem1(4), xxInt1, yyInt1, wwInt1;
        xxElem1 << x1, xA, xB, x4;
        yyElem1 << y1, yA, yB, y4;
        MathLib::LocalVector xxElem2(4), yyElem2(4), xxInt2, yyInt2, wwInt2;
        xxElem2 << xB, xA, x2, x3;
        yyElem2 << yB, yA, y2, y3;
        IntPoints2DRealElemQuad(xxElem1, yyElem1, nQnQ, xxIntRef, yyIntRef, wwIntRef, xxInt1, yyInt1, wwInt1);
        IntPoints2DRealElemQuad(xxElem2, yyElem2, nQnQ, xxIntRef, yyIntRef, wwIntRef, xxInt2, yyInt2, wwInt2);

        MathLib::LocalVector xxEleTotal(xxInt1.rows()*2);
        MathLib::LocalVector yyEleTotal(yyInt1.rows()*2);
        MathLib::LocalVector wwEleTotal(wwInt1.rows()*2);
        xxEleTotal << xxInt1, xxInt2;
        yyEleTotal << yyInt1, yyInt2;
        wwEleTotal << wwInt1, wwInt2;

        // Rotate integration points according to which of the nodes is really on the other side.
        MathLib::LocalMatrix T(2,2);
        if (Pos == 1)
            T << 0, 1, -1, 0; // Rotation -90 degree.
        else if (Pos == 2)
            T << 1, 0, 0, 1; // Nothing to do.

        MathLib::LocalMatrix xxyyIntTotal(2, xxEleTotal.rows());
        xxyyIntTotal << xxEleTotal.transpose(), yyEleTotal.transpose();
        MathLib::LocalMatrix HelpMat = T * xxyyIntTotal;
        xxIntTotal = HelpMat.row(0).transpose();
        yyIntTotal = HelpMat.row(1).transpose();
        wwIntTotal = wwEleTotal;
    }
}

void ShapeFctsXFEMSign(
        const MathLib::LocalVector &xxElem, const MathLib::LocalVector &yyElem, const MathLib::LocalVector &ffElem,
        const MathLib::LocalVector &NodesAct,
        const MathLib::LocalVector &xxIntRef, const MathLib::LocalVector &yyIntRef, const MathLib::LocalVector &wwIntRef,
        const size_t nQ,
        MathLib::LocalMatrix &N, MathLib::LocalMatrix &dNdx, MathLib::LocalMatrix &dNdy,
        MathLib::LocalMatrix &M, MathLib::LocalMatrix &dMdx, MathLib::LocalMatrix &dMdy,
        MathLib::LocalVector &xxInt, MathLib::LocalVector &yyInt, MathLib::LocalVector &wwInt, MathLib::LocalVector &ffInt)
{
    // Compute shape functions and enrichment functions at integration
    // points in the real element.

    // Get standard FE shape functions and level-set values at int.points.
    ShapeFctsStrdFEM(xxElem, yyElem, xxIntRef, yyIntRef, wwIntRef, nQ, N, dNdx, dNdy, xxInt, yyInt, wwInt);
    ffInt = N.transpose() * ffElem;

    // Define the enrichment functions M, dMdx, dMdy for enr. element nodes.
    M     = MathLib::LocalMatrix::Zero(4, nQ);
    dMdx  = MathLib::LocalMatrix::Zero(4, nQ);
    dMdy  = MathLib::LocalMatrix::Zero(4, nQ);
    if (NodesAct.array().sum() == 0)
        return; // Go back if enr. is not active in this element.

    for (size_t i=0; i<nQ; i++) {
        // Strd. FEM shape fcts. at int. point.
        MathLib::LocalVector NInt  = N.col(i);
        MathLib::LocalVector NxInt = dNdx.col(i);
        MathLib::LocalVector NyInt = dNdy.col(i);
        // Sign(Levelset) at int. point. and nodes.
        MathLib::LocalVector PsiElem = ffElem.array() / ffElem.array().abs(); //sign(ffElem)
        double PsiInt = ffInt(i) / ffInt.array().abs()(i); //sign(ffInt(i));
        double dPsidxInt = 0;
        double dPsidyInt = 0;
        // Enrichment function at int. point.
        M.col(i)    = NInt.array()*(PsiInt-PsiElem.array()).array();
        dMdx.col(i) = NxInt.array()*(PsiInt-PsiElem.array()).array() + NInt.array()*dPsidxInt;
        dMdy.col(i) = NyInt.array()*(PsiInt-PsiElem.array()).array() + NInt.array()*dPsidyInt;
    }

    MathLib::LocalVector Pos;
    //pos.array() = NodesAct.array() == 0;  //find(NodesAct == 0);
    for (int i=0; i<NodesAct.rows(); i++) {
        if (NodesAct(i)!=0) continue;
        M.row(i) *= 0.;
        dMdx.row(i) *= 0.;
        dMdy.row(i) *= 0.;
    }

}

void BuildMatRhs_Hooke(
        const MathLib::LocalMatrix &NMat, const MathLib::LocalMatrix &dNdxMat, const MathLib::LocalMatrix &dNdyMat,
        const MathLib::LocalMatrix &MMat, const MathLib::LocalMatrix &dMdxMat, const MathLib::LocalMatrix &dMdyMat,
        const MathLib::LocalVector &/*xxInt*/, const MathLib::LocalVector &/*yyInt*/, const MathLib::LocalVector &wwInt,
        const MathLib::LocalVector &ffInt, const Eigen::VectorXi &Nodes,
        double lambda1, double lambda2, double mu1, double mu2,
        double fx, double fy, size_t nQ, size_t NodeNum,
        MathLib::ILinearEquation &leqs)
{
    // Initialization.
    MathLib::LocalMatrix ElemMAT11 = MathLib::LocalMatrix::Zero(8, 8);
    MathLib::LocalMatrix ElemMAT12 = MathLib::LocalMatrix::Zero(8, 8);
    MathLib::LocalMatrix ElemMAT21 = MathLib::LocalMatrix::Zero(8, 8);
    MathLib::LocalMatrix ElemMAT22 = MathLib::LocalMatrix::Zero(8, 8);
    MathLib::LocalVector ElemRHS1  = MathLib::LocalVector::Zero(8);
    MathLib::LocalVector ElemRHS2  = MathLib::LocalVector::Zero(8);

    // Loop over integration points.
    for (size_t i=0; i<nQ; i++) {

        double lambda =.0 , mu = .0;
        if (ffInt(i) > 0) {
            lambda = lambda1;
            mu = mu1;
        } else if (ffInt(i) < 0) {
            lambda = lambda2;
            mu = mu2;
        }

        MathLib::LocalVector N  = NMat.col(i);
        MathLib::LocalVector Nx = dNdxMat.col(i);
        MathLib::LocalVector Ny = dNdyMat.col(i);

        MathLib::LocalVector M  = MMat.col(i);
        MathLib::LocalVector Mx = dMdxMat.col(i);
        MathLib::LocalVector My = dMdyMat.col(i);

        MathLib::LocalVector NxMx(2*Nx.rows());
        NxMx << Nx, Mx;
        MathLib::LocalVector NyMy(2*Ny.rows());
        NxMx << Ny, My;
        MathLib::LocalMatrix NxNxT = NxMx * NxMx.transpose();
        MathLib::LocalMatrix NxNyT = NxMx * NyMy.transpose();
        MathLib::LocalMatrix NyNxT = NyMy * NxMx.transpose();
        MathLib::LocalMatrix NyNyT = NyMy * NyMy.transpose();

        // Compute element matrices.
        ElemMAT11 += wwInt(i) * ( (lambda+2*mu) * NxNxT + mu * NyNyT );
        ElemMAT12 += wwInt(i) * ( lambda * NxNyT + mu  * NyNxT );
        ElemMAT21 += wwInt(i) * ( lambda * NyNxT + mu  * NxNyT );
        ElemMAT22 += wwInt(i) * ( (lambda+2*mu) * NyNyT + mu * NxNxT );

        // Compute right hand side.
        MathLib::LocalVector NM(2*N.rows());
        NM << N, M;
        ElemRHS1 += wwInt(i) * ( NM * fx );
        ElemRHS2 += wwInt(i) * ( NM * fy );

    }

    // Add element contribution to global matrix.
//    uuNodes = [Nodes          Nodes+2*NodeNum];
//    vvNodes = [Nodes+NodeNum  Nodes+3*NodeNum];
    std::vector<size_t> uuNodes(2*4), vvNodes(2*4);
    for (size_t i=0; i<4; i++) {
        uuNodes[i] = Nodes[i];
        uuNodes[i+4] = Nodes[i] + 2*NodeNum;
        vvNodes[i] = Nodes[i] + NodeNum;
        vvNodes[i+4] = Nodes[i] + 3*NodeNum;
    }
    leqs.addAsub(uuNodes, uuNodes, ElemMAT11);
    leqs.addAsub(uuNodes, vvNodes, ElemMAT12);
    leqs.addAsub(vvNodes, uuNodes, ElemMAT21);
    leqs.addAsub(vvNodes, vvNodes, ElemMAT22);
    leqs.addRHSsub(uuNodes, ElemRHS1);
    leqs.addRHSsub(vvNodes, ElemRHS2);
}

//}
