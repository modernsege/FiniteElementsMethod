#include <iostream>
#include <vector>
#include <math.h>  
#include"grid.h"
#include"element_2D.h"

using namespace std;



Node::Node(float x, float y, short int bcFlag) {
    this->x = x; //coordinate x of element's node
    this->y = y; //coordinate y of element's node
    this->bcFlag = bcFlag;
}

void Node::printCoordinates() {
    cout << "x = " << x << "\ny = " << y;
}

void Node::printBcFlag() {
    cout << "Boundry Condition flag: " << bcFlag;
}


Element::Element(int id1, int id2, int id3, int id4) {
    id[0] = id1;
    id[1] = id2;
    id[2] = id3;
    id[3] = id4;

    for (int iH = 0; iH < 4; iH++) {
        for (int jH = 0; jH < 4; jH++) {
            H[iH][jH] = 0;
            hbc[iH][jH] = 0;
        }
    }
}

void Element::printIdOfElement()
{
    cout << "\nId1: " << id[0] << "\nId2: " << id[1] << "\nId3: " << id[2] << "\nId4: " << id[3];
}




Grid::Grid(double H, double B, int nH, int nB) {
    this->H = H;
    this->B = B;
    this->nH = nH;
    this->nB = nB;

    this->nN = nH * nB;
    this->nE = (nH - 1) * (nB - 1);

    deltaH = H / (nH - 1);
    deltaB = B / (nB - 1);

    short int bcFlag = 0;

    for (int i = 0; i < nB; i++) {
        for (int j = 0; j < nH; j++) {
            bcFlag = 0;
            if (i == 0 || i == nB - 1 || j == 0 || j == nH - 1) {
                bcFlag = 1;
            }
            nodes.push_back(Node(i * deltaB, j * deltaH, bcFlag));
        }
    }

    int ID1, ID2, ID3, ID4;
    int elementinColumn = 1; //ordinal number of Element in each column
    ID1 = 1;
    for (int i = 0; i < nE; i++) {
        if (elementinColumn < nH) {
            ID2 = ID1 + nH;
            ID3 = ID2 + 1;
            ID4 = ID1 + 1;
            elements.push_back(Element(ID1, ID2, ID3, ID4));
            elementinColumn++;
        }
        else if (elementinColumn == nH) {
            ID1++;
            ID2 = ID1 + nH;
            ID3 = ID2 + 1;
            ID4 = ID1 + 1;
            elements.push_back(Element(ID1, ID2, ID3, ID4));
            elementinColumn = 2;
        }
        ID1++;
    }
}

void Grid::printElements() {
    for (int i = 0; i < elements.size(); i++) {
        cout << "\n\nElement: " << i + 1;
        elements[i].printIdOfElement();
    }
}

void Grid::printElementsAndNodesCoordinates() {
    for (int i = 0; i < elements.size(); i++) {
        cout << "\n_______________________";
        cout << "\n\nElement: " << i + 1;
        elements[i].printIdOfElement();
        cout << "\n\nNodes have coordinates: ";
        for (int j = 0; j < 4; j++) {
            cout << "\nID" << j + 1 << ": " << endl;
            nodes[elements[i].id[j] - 1].printCoordinates();
        }
    }
}

void Grid::printNodes() {
    for (int i = 0; i < nodes.size(); i++) {
        cout << "\n\nNode: " << i + 1;
        nodes[i].printCoordinates();
    }
    for (int i = 0; i < nodes.size(); i++) {
        cout << "\n\nNode: " << i + 1;
        nodes[i].printBcFlag();
    }
}






void Grid::calculateJacobianForGrid(int nPc) {
    Jacobian newJacobian = Jacobian();
    double detJ;
    double detJinv;
    double multipliedMatrix[2][2];
    Jacobian temp;

    if (nPc == 4) {
        Element4_2D newElement = Element4_2D();
        for (int i = 0; i < this->nE; i++) { //i -> number of element 
            for (int j = 0; j < nPc; j++) {//j -> number of integration point in one element

                Jacobian temp;
                Jacobian newJacobian = Jacobian();
                temp = newElement.jacobian(i, j, newElement, nodes, elements, newJacobian);
                cout << "Element: " << i + 1 << " Integration Point: " << j + 1 << endl;
                cout << "\nJacobian: " << endl;
                for (int k = 0; k < 2; k++) {
                    for (int l = 0; l < 2; l++) {
                        cout << temp.jacobian[k][l] << "\t";
                    }
                    cout << "\n";
                }
                cout << "\n";
                detJ = (temp.jacobian[0][0] * temp.jacobian[1][1]) - (temp.jacobian[1][0] * temp.jacobian[0][1]);
                detJinv = 1.0 / detJ;
                cout << "detJ = " << detJ << endl;
                cout << "1/detJ = " << detJinv << endl << endl;

                cout << "Multiplied Matrix: " << endl;
                for (int k = 0; k < 2; k++) {
                    for (int l = 0; l < 2; l++) {
                        multipliedMatrix[k][l] = temp.jacobian_inv[k][l] * detJinv;
                        cout << multipliedMatrix[k][l] << "\t";
                    }
                    cout << "\n";
                }




                //================================ calculate dNi/dNx, dNi/dy ===============================
                double dNdx[4][4]; //store value of dNi/dNx and its transposition
                double dNdy[4][4]; //store value of dNi/dy and its transposition


               // Element4_2D elem = Element4_2D();
                for (int iH = 0; iH < 4; iH++) {
                    for (int jH = 0; jH < 4; jH++) {
                        dNdx[iH][jH] = (multipliedMatrix[0][0] * newElement.dNdKSI[j][iH] + multipliedMatrix[0][1] * newElement.dNdETA[j][iH]) * (multipliedMatrix[0][0] * newElement.dNdKSI[j][jH] + multipliedMatrix[0][1] * newElement.dNdETA[j][jH]);
                        dNdy[iH][jH] = (multipliedMatrix[1][0] * newElement.dNdKSI[j][iH] + multipliedMatrix[1][1] * newElement.dNdETA[j][iH]) * (multipliedMatrix[1][0] * newElement.dNdKSI[j][jH] + multipliedMatrix[1][1] * newElement.dNdETA[j][jH]);
                    }
                }
                /* for (int iH = 0; iH < 4; iH++) {
                     for (int jH = 0; jH < 4; jH++) {
                         cout << dNdy[iH][jH] << "\t";
                     }
                     cout << endl;
                 }*/

                double Hmatrix[4][4];
                double k = elements[i].k = 30.0;


                cout << "\nH matrix in " << j + 1 << " integration point\n";
                for (int iH = 0; iH < 4; iH++) {
                    for (int jH = 0; jH < 4; jH++) {
                        Hmatrix[iH][jH] = ((dNdx[iH][jH] + dNdy[iH][jH]) * detJ * k);
                        cout << Hmatrix[iH][jH] << "  ";
                        elements[i].H[iH][jH] += Hmatrix[iH][jH];
                    }
                    cout << "\n";
                }
                cout << "===========================" << endl;

            }
            cout << "===========================\n";
            cout << "\nH matrix for element " << i + 1 << "\n";
            for (int iH = 0; iH < 4; iH++) {
                for (int jH = 0; jH < 4; jH++) {
                    cout << elements[i].H[iH][jH] << "  ";
                }
                cout << "\n";
            }
            cout << "\n===========================\n===========================\n";


        }

        //HBC


        for (int i = 0; i < this->nE; i++) {

            double alpha = elements[i].alpha = 25.0;

            for (int i1 = 0; i1 < 4; i1++) {
                if (nodes[elements[i].id[i1] - 1].bcFlag == 1 && nodes[elements[i].id[(i1 + 1) % 4] - 1].bcFlag == 1) {
                    //double Pc1[4] = { 0,0,0,0 };
                    //double Pc2[4] = { 0,0,0,0 };
                    double Pc1[4] = { 0,0,0,0 };
                    double Pc2[4] = { 0,0,0,0 };


                    for (int l = 0; l < 4; l++) {
                        Pc1[l] = newElement.sides[i1].N[0][l];
                        Pc2[l] = newElement.sides[i1].N[1][l];
                    }




                    //**** STARA WERSJA - bez wartoœci wpisanych w strukture elementu 2D
                    //Pc1[i1] = (1.0 - 1.0 / sqrt(3)) / 2.0; //1 to wartosc f. ksztaltu w swoim wezle odjac 1/sqrt(3) czyli przesuniecie wzgledem 100% wartoœci f. //ksztaltu,dzielone .na 2 bo jestesmy w ukladzie (-1,1) a chce sie odniesc do ukladu (0, 1)
                    //Pc1[(i1 + 1) % 4] = 1 - Pc1[i1];
                    //
                    //
                    //Pc2[(i1 + 1) % 4] = Pc1[i1];
                    //Pc2[i1] = Pc1[(i1 + 1) % 4];
                    //****************

                    double detJ = sqrt(pow(nodes[elements[i].id[i1] - 1].x - nodes[elements[i].id[(i1 + 1) % 4] - 1].x, 2) + pow(nodes[elements[i].id[i1] - 1].y - nodes[elements[i].id[(i1 + 1) % 4] - 1].y, 2)) / 2.0; //Jacobian = Dlugosc bokku/2, dlugosc boku liczona z wekotra, potem dlugosc wektora

                    for (int hbcRow = 0; hbcRow < 4; hbcRow++) {
                        for (int hbcColumn = 0; hbcColumn < 4; hbcColumn++) {
                            elements[i].hbc[hbcRow][hbcColumn] += newElement.sides[i1].wages[0] * Pc1[hbcRow] * Pc1[hbcColumn] * detJ * elements[i].alpha;
                            elements[i].hbc[hbcRow][hbcColumn] += newElement.sides[i1].wages[1] * Pc2[hbcRow] * Pc2[hbcColumn] * detJ * elements[i].alpha;
                        }
                    }

                }

            }
            cout << "===========================\n";
            cout << "\nHBC matrix for element " << i + 1 << endl;
            for (int hbcRow = 0; hbcRow < 4; hbcRow++) {
                for (int hbcColumn = 0; hbcColumn < 4; hbcColumn++) {
                    cout << elements[i].hbc[hbcRow][hbcColumn] << "\t";
                }
                cout << "\n";
            }
            cout << "\n===========================\n";

        }


    }
    else if (nPc == 9) {
        Element9_2D newElement = Element9_2D();
        for (int i = 0; i < this->nE; i++) { //i -> number of element 
            for (int j = 0; j < nPc; j++) { //j -> number of integration point in one element
                temp = newElement.jacobian(i, j, newElement, nodes, elements, newJacobian);
                cout << "Element: " << i + 1 << " Integration Point: " << j + 1 << endl;
                cout << "\nJacobian: " << endl;
                for (int k = 0; k < 2; k++) {
                    for (int l = 0; l < 2; l++) {
                        cout << temp.jacobian[k][l] << "\t";
                    }
                    cout << "\n";
                }
                cout << "\n";
                detJ = (temp.jacobian[0][0] * temp.jacobian[1][1]) - (temp.jacobian[1][0] * temp.jacobian[0][1]);
                detJinv = 1.0 / detJ;
                cout << "detJ = " << detJ << endl;
                cout << "1/detJ = " << detJinv << endl << endl;

                cout << "Multiplied Matrix: " << endl;
                for (int k = 0; k < 2; k++) {
                    for (int l = 0; l < 2; l++) {
                        multipliedMatrix[k][l] = temp.jacobian_inv[k][l] * detJinv;
                        cout << multipliedMatrix[k][l] << "\t";
                    }
                    cout << "\n";
                }
                cout << "===========================" << endl;
                //================================ calculate dNi/dNx, dNi/dy ===============================
                double dNdx[4][4]; //store value of dNi/dNx and its transposition
                double dNdy[4][4]; //store value of dNi/dy and its transposition


                double weight = 1;
                if (nPc == 9) {
                    if (j == 0 || j == 2 || j == 6 || j == 8) {
                        weight = 25.0 / 81.0;
                    }
                    else if (j == 4) {
                        weight = 64.0 / 81.0;
                    }
                    else
                    {
                        weight = 40.0 / 81.0;
                    }
                }

                // Element4_2D elem = Element4_2D();
                for (int iH = 0; iH < 4; iH++) {
                    for (int jH = 0; jH < 4; jH++) {
                        dNdx[iH][jH] = (multipliedMatrix[0][0] * newElement.dNdKSI[j][iH] + multipliedMatrix[0][1] * newElement.dNdETA[j][iH]) * (multipliedMatrix[0][0] * newElement.dNdKSI[j][jH] + multipliedMatrix[0][1] * newElement.dNdETA[j][jH]);
                        dNdy[iH][jH] = (multipliedMatrix[1][0] * newElement.dNdKSI[j][iH] + multipliedMatrix[1][1] * newElement.dNdETA[j][iH]) * (multipliedMatrix[1][0] * newElement.dNdKSI[j][jH] + multipliedMatrix[1][1] * newElement.dNdETA[j][jH]);
                    }
                }
                double Hmatrix[4][4];
                double k = elements[i].k = 30.0;


                cout << "\nH matrix in " << j + 1 << " integration point\n";
                for (int iH = 0; iH < 4; iH++) {
                    for (int jH = 0; jH < 4; jH++) {
                        Hmatrix[iH][jH] = weight * ((dNdx[iH][jH] + dNdy[iH][jH]) * detJ * k);
                        cout << Hmatrix[iH][jH] << "  ";
                        elements[i].H[iH][jH] += Hmatrix[iH][jH];
                    }
                    cout << "\n";
                }





            }
            cout << "===========================\n";
            cout << "\nH matrix for element " << i + 1 << "\n";
            for (int iH = 0; iH < 4; iH++) {
                for (int jH = 0; jH < 4; jH++) {
                    cout << elements[i].H[iH][jH] << "  ";
                }
                cout << "\n";
            }
            cout << "\n===========================\n===========================\n";
        }


        for (int i = 0; i < this->nE; i++) {

            double alpha = elements[i].alpha = 25.0;

            for (int i1 = 0; i1 < 4; i1++) {
                if (nodes[elements[i].id[i1] - 1].bcFlag == 1 && nodes[elements[i].id[(i1 + 1) % 4] - 1].bcFlag == 1) {
                    double Pc1[4] = { 0,0,0,0 };
                    double Pc2[4] = { 0,0,0,0 };
                    double Pc3[4] = { 0,0,0,0 };


                    for (int l = 0; l < 4; l++) {
                        Pc1[l] = newElement.sides[i1].N[0][l];
                        Pc2[l] = newElement.sides[i1].N[1][l];
                        Pc3[l] = newElement.sides[i1].N[2][l];
                    }



                   //Pc1[i1] = (1.0 - sqrt(3.0 / 5.0)) / 2.0;
                   //Pc1[(i1 + 1) % 4] = 1 - Pc1[i1];
                   //
                   //
                   //Pc2[(i1 + 1) % 4] = (1.0) / 2.0;
                   //Pc2[i1] = (1.0) / 2.0;
                   //
                   //Pc3[(i1 + 1) % 4] = Pc1[i1];
                   //Pc3[i1] = Pc1[(i1 + 1) % 4];

                    double detJ = sqrt(pow(nodes[elements[i].id[i1] - 1].x - nodes[elements[i].id[(i1 + 1) % 4] - 1].x, 2) + pow(nodes[elements[i].id[i1] - 1].y - nodes[elements[i].id[(i1 + 1) % 4] - 1].y, 2)) / 2.0; //Jacobian = Dlugosc bokku/2, dlugosc boku liczona z wekotra, potem dlugosc wektora

                    for (int hbcRow = 0; hbcRow < 4; hbcRow++) {
                        for (int hbcColumn = 0; hbcColumn < 4; hbcColumn++) {
                            elements[i].hbc[hbcRow][hbcColumn] += newElement.sides->wages[0] * Pc1[hbcRow] * Pc1[hbcColumn] * detJ * elements[i].alpha;
                            elements[i].hbc[hbcRow][hbcColumn] += newElement.sides->wages[1] * Pc2[hbcRow] * Pc2[hbcColumn] * detJ * elements[i].alpha;
                            elements[i].hbc[hbcRow][hbcColumn] += newElement.sides->wages[2] * Pc3[hbcRow] * Pc3[hbcColumn] * detJ * elements[i].alpha;
                        }
                    }

                }

            }

            cout << "HBC matrix for element " << i + 1 << endl;
            for (int hbcRow = 0; hbcRow < 4; hbcRow++) {
                for (int hbcColumn = 0; hbcColumn < 4; hbcColumn++) {
                    cout << elements[i].hbc[hbcRow][hbcColumn] << "\t";
                }
                cout << "\n";
            }

        }
    }

}
