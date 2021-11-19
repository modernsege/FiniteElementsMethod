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


    this->sides[0].integrationPointsBC[0].x = -1.0 / sqrt(3.0);
    this->sides[0].integrationPointsBC[0].y = -1.0;
    this->sides[0].integrationPointsBC[1].x = 1.0 / sqrt(3.0);
    this->sides[0].integrationPointsBC[1].y = -1.0;
                                    
    this->sides[1].integrationPointsBC[0].x = 1.0;
    this->sides[1].integrationPointsBC[0].y = -1.0 / sqrt(3.0);
    this->sides[1].integrationPointsBC[1].x = 1.0;
    this->sides[1].integrationPointsBC[1].y = 1.0 / sqrt(3.0);
                                    
    this->sides[2].integrationPointsBC[0].x = 1.0 / sqrt(3.0);
    this->sides[2].integrationPointsBC[0].y = 1.0;
    this->sides[2].integrationPointsBC[1].x = -1.0 / sqrt(3.0);
    this->sides[2].integrationPointsBC[1].y = 1.0; 
                                    
    this->sides[3].integrationPointsBC[0].x = -1.0;
    this->sides[3].integrationPointsBC[0].y = 1.0/sqrt(3.0);
    this->sides[3].integrationPointsBC[1].x = -1.0;
    this->sides[3].integrationPointsBC[1].y = -1.0 / sqrt(3.0);
                  


    for (int i = 0; i < 4; i++) {
        this->sides[i].N[0][0] = 0.25 * (1.0 - this->sides[i].integrationPointsBC[0].x) * (1.0 - this->sides[i].integrationPointsBC[0].y);
        this->sides[i].N[0][1] = 0.25 * (1.0 + this->sides[i].integrationPointsBC[0].x) * (1.0 - this->sides[i].integrationPointsBC[0].y);
        this->sides[i].N[0][2] = 0.25 * (1.0 + this->sides[i].integrationPointsBC[0].x) * (1.0 + this->sides[i].integrationPointsBC[0].y);
        this->sides[i].N[0][3] = 0.25 * (1.0 - this->sides[i].integrationPointsBC[0].x) * (1.0 + this->sides[i].integrationPointsBC[0].y);

        this->sides[i].N[1][0] = 0.25 * (1.0 - this->sides[i].integrationPointsBC[1].x) * (1.0 - this->sides[i].integrationPointsBC[1].y);
        this->sides[i].N[1][1] = 0.25 * (1.0 + this->sides[i].integrationPointsBC[1].x) * (1.0 - this->sides[i].integrationPointsBC[1].y);
        this->sides[i].N[1][2] = 0.25 * (1.0 + this->sides[i].integrationPointsBC[1].x) * (1.0 + this->sides[i].integrationPointsBC[1].y);
        this->sides[i].N[1][3] = 0.25 * (1.0 - this->sides[i].integrationPointsBC[1].x) * (1.0 + this->sides[i].integrationPointsBC[1].y);
        
        this->sides[i].wages[0] = 1.0;
        this->sides[i].wages[1] = 1.0;
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

    double dXdKsi = 0;
    double dYdEta = 0;
    double dYdKsi = 0;
    double dXdEta = 0;


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


Element9_2D::Element9_2D() {
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

    this->sides[0].integrationPointsBC[0].x = -sqrt(3.0/5.0);
    this->sides[0].integrationPointsBC[0].y = -1;
    this->sides[0].integrationPointsBC[1].x = 0.0;
    this->sides[0].integrationPointsBC[1].y = -1;
    this->sides[0].integrationPointsBC[2].x = sqrt(3.0 / 5.0);
    this->sides[0].integrationPointsBC[2].y = -1;
                                    
    this->sides[1].integrationPointsBC[0].x = 1;
    this->sides[1].integrationPointsBC[0].y = -sqrt(3.0 / 5.0);
    this->sides[1].integrationPointsBC[1].x = 1;
    this->sides[1].integrationPointsBC[1].y = 0.0;
    this->sides[1].integrationPointsBC[2].x = 1;
    this->sides[1].integrationPointsBC[2].y = sqrt(3.0 / 5.0);
                                    
    this->sides[2].integrationPointsBC[0].x = sqrt(3.0 / 5.0);
    this->sides[2].integrationPointsBC[0].y = 1;
    this->sides[2].integrationPointsBC[1].x = 0.0;
    this->sides[2].integrationPointsBC[1].y = 1;
    this->sides[2].integrationPointsBC[2].x = -sqrt(3.0 / 5.0);
    this->sides[2].integrationPointsBC[2].y = 1;
                                    
    this->sides[3].integrationPointsBC[0].x = -1;
    this->sides[3].integrationPointsBC[0].y = sqrt(3.0 / 5.0);
    this->sides[3].integrationPointsBC[1].x = -1;
    this->sides[3].integrationPointsBC[1].y = 0.0;
    this->sides[3].integrationPointsBC[2].x = -1;
    this->sides[3].integrationPointsBC[2].y = -sqrt(3.0 / 5.0);

    for (int i = 0; i < 4; i++) {

        this->sides[i].N[0][0] = 0.25 * (1.0 - this->sides[i].integrationPointsBC[0].x) * (1.0 - this->sides[i].integrationPointsBC[0].y);
        this->sides[i].N[0][1] = 0.25 * (1.0 + this->sides[i].integrationPointsBC[0].x) * (1.0 - this->sides[i].integrationPointsBC[0].y);
        this->sides[i].N[0][2] = 0.25 * (1.0 + this->sides[i].integrationPointsBC[0].x) * (1.0 + this->sides[i].integrationPointsBC[0].y);
        this->sides[i].N[0][3] = 0.25 * (1.0 - this->sides[i].integrationPointsBC[0].x) * (1.0 + this->sides[i].integrationPointsBC[0].y);

        this->sides[i].N[1][0] = 0.25 * (1.0 - this->sides[i].integrationPointsBC[1].x) * (1.0 - this->sides[i].integrationPointsBC[1].y);
        this->sides[i].N[1][1] = 0.25 * (1.0 + this->sides[i].integrationPointsBC[1].x) * (1.0 - this->sides[i].integrationPointsBC[1].y);
        this->sides[i].N[1][2] = 0.25 * (1.0 + this->sides[i].integrationPointsBC[1].x) * (1.0 + this->sides[i].integrationPointsBC[1].y);
        this->sides[i].N[1][3] = 0.25 * (1.0 - this->sides[i].integrationPointsBC[1].x) * (1.0 + this->sides[i].integrationPointsBC[1].y);

        this->sides[i].N[2][0] = 0.25 * (1.0 - this->sides[i].integrationPointsBC[2].x) * (1.0 - this->sides[i].integrationPointsBC[2].y);
        this->sides[i].N[2][1] = 0.25 * (1.0 + this->sides[i].integrationPointsBC[2].x) * (1.0 - this->sides[i].integrationPointsBC[2].y);
        this->sides[i].N[2][2] = 0.25 * (1.0 + this->sides[i].integrationPointsBC[2].x) * (1.0 + this->sides[i].integrationPointsBC[2].y);
        this->sides[i].N[2][3] = 0.25 * (1.0 - this->sides[i].integrationPointsBC[2].x) * (1.0 + this->sides[i].integrationPointsBC[2].y);
        
        this->sides[i].wages[0] = 5.0 / 9.0;
        this->sides[i].wages[1] = 8.0 / 9.0;
        this->sides[i].wages[2] = 5.0 / 9.0;
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
    double dXdKsi = 0;
    double dYdEta = 0;
    double dYdKsi = 0;
    double dXdEta = 0;


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





