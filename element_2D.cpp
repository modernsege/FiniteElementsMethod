#include <iostream>
#include <vector>
#include <math.h>  
#include "element_2D.h"
#include "grid.h"

using namespace std;

Element4_2D::Element4_2D() {
    for (int i = 0; i < 4; i++) {                                       //4pc   3pc
        if (i == 2 || i == 3) {//i == punkt calkowania                  //1pc   2pc
            eta = 1 / sqrt(3);
        }
        else {
            eta = -1 / sqrt(3);
        }

        if (i == 1 || i == 2) {
            ksi = 1 / sqrt(3);
        }
        else {
            ksi = -1 / sqrt(3);
        }

        dNdKSI[i][0] = -0.25 * (1.0 - eta);
        dNdKSI[i][1] = 0.25 * (1.0 - eta);
        dNdKSI[i][2] = 0.25 * (1.0 + eta);
        dNdKSI[i][3] = -0.25 * (1.0 + eta);

        dNdETA[i][0] = -0.25 * (1.0 - ksi);
        dNdETA[i][1] = -0.25 * (1.0 + ksi);
        dNdETA[i][2] = 0.25 * (1.0 + ksi);
        dNdETA[i][3] = 0.25 * (1.0 - ksi);
    }
}

void Element4_2D::printFor4IntegrationPoints() {
    cout << "\ndN/dKsi oraz dN/dEta dla 4 punktow calkowania:" << endl;
    cout << "Funkcje ksztaltu dN/dKsi" << endl;
    for (int i = 0; i < 4; i++) {
        cout << i + 1 << " ";
        for (int j = 0; j < 4; j++) {
            cout << dNdKSI[i][j] << " ";
        }
        cout << endl;
    }

    cout << "\nFunkcje ksztaltu dN/dEta" << endl;
    for (int i = 0; i < 4; i++) {
        cout << i + 1 << " ";
        for (int j = 0; j < 4; j++) {
            cout << dNdETA[i][j] << " ";
        }
        cout << endl;
    }
}


Jacobian Element4_2D::jacobian(int i, int j, Element4_2D element4_2D, std::vector<Node> nodes, std::vector<Element> elements, Jacobian struct_jacobian) {
    // FOR TEST 
    //double xHardcode[4] = { 0, 0.025, 0.025, 0 };
    //double yHardcode[4] = { 0, 0, 0.025, 0.025 };
    double dXdKsi = 0;
    double dYdEta = 0;
    double dYdKsi = 0;
    double dXdEta = 0;


    /*for (int k = 0; k < 4; k++) {
        dXdKsi += element4_2D.dNdKSI[j][k] * xHardcode[k];
        dYdEta += element4_2D.dNdETA[j][k] * yHardcode[k];
        dYdKsi += element4_2D.dNdKSI[j][k] * yHardcode[k];
        dXdEta += element4_2D.dNdETA[j][k] * xHardcode[k];
    }*/

    //cout << endl << "!!!   "; elements[i].printIdOfElement(); cout << "   !!!" << endl;
    for (int k = 0; k < 4; k++) {
        dXdKsi += element4_2D.dNdKSI[j][k] * nodes[elements[i].id[k] - 1].x;
        dYdEta += element4_2D.dNdETA[j][k] * nodes[elements[i].id[k] - 1].y;
        dYdKsi += element4_2D.dNdKSI[j][k] * nodes[elements[i].id[k] - 1].y;
        dXdEta += element4_2D.dNdETA[j][k] * nodes[elements[i].id[k] - 1].x;
    }

    struct_jacobian.jacobian[0][0] = dXdKsi;
    struct_jacobian.jacobian[0][1] = dYdKsi;
    struct_jacobian.jacobian[1][0] = dXdEta;
    struct_jacobian.jacobian[1][1] = dYdEta;

    struct_jacobian.jacobian_inv[0][0] = dYdEta;
    struct_jacobian.jacobian_inv[0][1] = -dYdKsi;
    struct_jacobian.jacobian_inv[1][0] = -dXdEta;
    struct_jacobian.jacobian_inv[1][1] = dXdKsi;

    return struct_jacobian;
}
  

Element9_2D:: Element9_2D() {
    for (int i = 0; i < 9; i++) {
        if (i == 0 || i == 3 || i == 6) {              //7pc    8pc    9pc 
           ksi = -sqrt(3.0 / 5.0);                     //4pc    5pc    6pc
        }                                              //1pc    2pc    3pc
        else if (i == 1 || i == 4 || i == 7) {
            ksi = 0;
        }
        else {
            ksi = sqrt(3.0 / 5.0);
        }

        if (i == 0 || i == 1 || i == 2) {
            eta = -sqrt(3.0 / 5.0);
        }
        else if (i == 3 || i == 4 || i == 5) {
            eta = 0;
        }
        else {
            eta = sqrt(3.0 / 5.0);
        }

        dNdKSI[i][0] = -0.25 * (1.0 - eta);
        dNdKSI[i][1] = 0.25 * (1.0 - eta);
        dNdKSI[i][2] = 0.25 * (1.0 + eta);
        dNdKSI[i][3] = -0.25 * (1.0 + eta);

        dNdETA[i][0] = -0.25 * (1.0 - ksi);
        dNdETA[i][1] = -0.25 * (1.0 + ksi);
        dNdETA[i][2] = 0.25 * (1.0 + ksi);
        dNdETA[i][3] = 0.25 * (1.0 - ksi);
    }
}

void Element9_2D::printFor9IntegrationPoints() {
    cout << "\ndN/dKsi oraz dN/dEta dla 9 punktow calkowania:" << endl;
    cout << "Funkcje ksztaltu dN/dKsi" << endl;
    for (int i = 0; i < 9; i++) {
        cout << i + 1 << " ";
        for (int j = 0; j < 4; j++) {
            cout << dNdKSI[i][j] << " ";
        }
        cout << endl;
    }

    cout << "\nFunkcje ksztaltu dN/dEta" << endl;
    for (int i = 0; i < 9; i++) {
        cout << i + 1 << " ";
        for (int j = 0; j < 4; j++) {
            cout << dNdETA[i][j] << " ";
        }
        cout << endl;
     }
}

Jacobian Element9_2D::jacobian(int i, int j, Element9_2D element9_2D, std::vector<Node> nodes, std::vector<Element> elements, Jacobian struct_jacobian) {
   // FOR TEST 
   // double xHardcode[4] = { 0, 0.025, 0.025, 0 };
   // double yHardcode[4] = { 0, 0, 0.025, 0.025 };
    double dXdKsi = 0;
    double dYdEta = 0;
    double dYdKsi = 0;
    double dXdEta = 0;


    /*for (int k = 0; k < 4; k++) {
        dXdKsi += element4_2D.dNdKSI[j][k] * xHardcode[k];
        dYdEta += element4_2D.dNdETA[j][k] * yHardcode[k];
        dYdKsi += element4_2D.dNdKSI[j][k] * yHardcode[k];
        dXdEta += element4_2D.dNdETA[j][k] * xHardcode[k];
    }*/

    //cout << endl << "!!!   "; elements[i].printIdOfElement(); cout << "   !!!" << endl;
    for (int k = 0; k < 4; k++) {
        dXdKsi += element9_2D.dNdKSI[j][k] * nodes[elements[i].id[k] - 1].x;
        dYdEta += element9_2D.dNdETA[j][k] * nodes[elements[i].id[k] - 1].y;
        dYdKsi += element9_2D.dNdKSI[j][k] * nodes[elements[i].id[k] - 1].y;
        dXdEta += element9_2D.dNdETA[j][k] * nodes[elements[i].id[k] - 1].x;
    }

    struct_jacobian.jacobian[0][0] = dXdKsi;
    struct_jacobian.jacobian[0][1] = dYdKsi;
    struct_jacobian.jacobian[1][0] = dXdEta;
    struct_jacobian.jacobian[1][1] = dYdEta;

    struct_jacobian.jacobian_inv[0][0] = dYdEta;
    struct_jacobian.jacobian_inv[0][1] = -dYdKsi;
    struct_jacobian.jacobian_inv[1][0] = -dXdEta;
    struct_jacobian.jacobian_inv[1][1] = dXdKsi;

    return struct_jacobian;
}



