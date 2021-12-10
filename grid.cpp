#include <iostream>
#include <vector>
#include <math.h>  
#include"grid.h"
#include"element_2D.h"

using namespace std;

///// NODE

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



//ELEMENT
Element::Element(int id1, int id2, int id3, int id4, double alpha, double k, double temprature, double specificHeat, double density) {
    id[0] = id1;
    id[1] = id2;
    id[2] = id3;
    id[3] = id4;

    for (int iH = 0; iH < 4; iH++) {
        for (int jH = 0; jH < 4; jH++) {
            H[iH][jH] = 0;
            hbc[iH][jH] = 0;
            C[iH][jH] = 0;
        }
    }

    this->alpha = alpha;
    this->k = k;
    this->temperature = temprature;
    this->c = specificHeat;
    this->ro = density;
}

void Element::printIdOfElement()
{
    cout << "\nId1: " << id[0] << "\nId2: " << id[1] << "\nId3: " << id[2] << "\nId4: " << id[3];
}









//GRID

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



void Grid::agregation() {
    //HGlobalne
    for (int elementNum = 0; elementNum < this->nE; elementNum++) {
        for (int iH = 0; iH < 4; iH++) {
            for (int jH = 0; jH < 4; jH++) {
                HGlobalne[elements[elementNum].id[iH] - 1][elements[elementNum].id[jH] - 1] += elements[elementNum].H[iH][jH] + elements[elementNum].hbc[iH][jH]; //H Global
                CGlobal[elements[elementNum].id[iH] - 1][elements[elementNum].id[jH] - 1] += elements[elementNum].C[iH][jH]; //C Global
            }
        }
    }
}

void Grid::pVectorAgregation() {
    for (int elementNum = 0; elementNum < this->nE; elementNum++) {      
        for (int pRows = 0; pRows < 4; pRows++) {
            PGlobal[elements[elementNum].id[pRows]-1] += elements[elementNum].p_vector[pRows];
        }
    }
}


void sulution_t();


void Grid::print_H_matrix() {
    for (int i = 0; i < this->nE; i++) {
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

}

void Grid::print_C_matrix() {
    for (int i = 0; i < this->nE; i++) {
        cout << "===========================\n";
        cout << "\nC matrix for element " << i + 1 << "\n";
        for (int iH = 0; iH < 4; iH++) {
            for (int jH = 0; jH < 4; jH++) {
                cout << elements[i].C[iH][jH] << "  ";
            }
            cout << "\n";
        }
        cout << "\n===========================\n===========================\n";
    }
}

void Grid::print_HBC_matrix() {
    for (int i = 0; i < this->nE; i++) {
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


void Grid::print_P_vector() {
    for (int i = 0; i < this->nE; i++) {
        cout << "===========================\n";
        cout << "\nP vector for element " << i + 1 << endl;
        for (int pRow = 0; pRow < 4; pRow++) {

            cout << elements[i].p_vector[pRow] << "\t";
        }
        cout << "\n===========================\n";
    }
}

void Grid::print_H_with_HBC_Global() {
    cout << "\nH GLOBAL: \n";
    for (int iH = 0; iH < HGlobalne.size(); iH++) {
        for (int jH = 0; jH < HGlobalne[iH].size(); jH++) {
            cout << HGlobalne[iH][jH] << "\t";
        }
        cout << endl;
    }
}

void Grid::print_C_Global() {
    cout << "\nC GLOBAL: \n";
    for (int iC = 0; iC < CGlobal.size(); iC++) {
        for (int jC = 0; jC < CGlobal[iC].size(); jC++) {
            cout << CGlobal[iC][jC] << "\t";
        }
        cout << endl;
    }
}


void Grid::print_P_Global_vector() {
    cout << "\nP GLOBAL: \n";
    for (int iH = 0; iH < PGlobal.size(); iH++) {
        cout << PGlobal[iH] << "\t";
        cout << endl;
    }
}


Grid::Grid(double H, double B, int nH, int nB, int nPc, double alpha, double k, double temprature,double specificHeat, double density) {

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
            elements.push_back(Element(ID1, ID2, ID3, ID4, alpha, k, temprature, specificHeat, density));
            elementinColumn++;
        }
        else if (elementinColumn == nH) {
            ID1++;
            ID2 = ID1 + nH;
            ID3 = ID2 + 1;
            ID4 = ID1 + 1;
            elements.push_back(Element(ID1, ID2, ID3, ID4, alpha, k, temprature, specificHeat, density));
            elementinColumn = 2;
        }
        ID1++;
    }


    HGlobalne.resize(nN, vector<double>(nN, 0)); //macierz H globalne dla grida ma miec wymiary nN na nN - wstepnie wypelniona zerami
    PGlobal.resize(nN, 0);
    CGlobal.resize(nN, vector<double>(nN, 0)); //macierz C globalne dla grida ma miec wymiary nN na nN - wstepnie wypelniona zerami


    Jacobian newJacobian = Jacobian();
    double detJ;
    double detJinv;
    double multipliedMatrix[2][2];
    Jacobian jacobianStructure;

    if (nPc == 4) {
        //////////calculate Jacobian for each element
        Element4_2D newElement = Element4_2D();
        for (int i = 0; i < this->nE; i++) { //i -> number of element 
            for (int j = 0; j < nPc; j++) {//j -> number of integration point in one element

                Jacobian newJacobian = Jacobian();
                jacobianStructure = newElement.jacobian(i, j, newElement, nodes, elements, newJacobian);
                detJ = (jacobianStructure.jacobian[0][0] * jacobianStructure.jacobian[1][1]) - (jacobianStructure.jacobian[1][0] * jacobianStructure.jacobian[0][1]);
                detJinv = 1.0 / detJ;
                for (int k = 0; k < 2; k++) {
                    for (int l = 0; l < 2; l++) {
                        multipliedMatrix[k][l] = jacobianStructure.jacobian_inv[k][l] * detJinv;
                    }
                }

                //Print jacobian in each element and integration point
              /* cout << "Element: " << i + 1 << " Integration Point: " << j + 1 << endl;
               cout << "\nJacobian: " << endl;
               for (int k = 0; k < 2; k++) {
                   for (int l = 0; l < 2; l++) {
                       cout << jacobianStructure.jacobian[k][l] << "\t";
                   }
                   cout << "\n";
               }
               cout << "\n";

               cout << "detJ = " << detJ << endl;
               cout << "1/detJ = " << detJinv << endl << endl;
              
               cout << "Multiplied Matrix: " << endl;
                for (int k = 0; k < 2; k++) {
                    for (int l = 0; l < 2; l++) {
                       cout << multipliedMatrix[k][l] << "\t";
                    }
                    cout << "\n";
                }
                */



                //================================ calculate dNi/dNx, dNi/dy ===============================
                double dNdx[4][4]; //store value of dNi/dNx and its transposition
                double dNdy[4][4]; //store value of dNi/dy and its transposition
                double c_matrix_temp[4][4]; //store value of N_valuses_in_integration_pionts and its transposition


                for (int iH = 0; iH < 4; iH++) {
                    for (int jH = 0; jH < 4; jH++) {
                        dNdx[iH][jH] = (multipliedMatrix[0][0] * newElement.dNdKSI[j][iH] + multipliedMatrix[0][1] * newElement.dNdETA[j][iH]) * (multipliedMatrix[0][0] * newElement.dNdKSI[j][jH] + multipliedMatrix[0][1] * newElement.dNdETA[j][jH]);
                        dNdy[iH][jH] = (multipliedMatrix[1][0] * newElement.dNdKSI[j][iH] + multipliedMatrix[1][1] * newElement.dNdETA[j][iH]) * (multipliedMatrix[1][0] * newElement.dNdKSI[j][jH] + multipliedMatrix[1][1] * newElement.dNdETA[j][jH]);
                        c_matrix_temp[iH][jH] = newElement.N_valuses_in_integration_pionts[j][iH] * newElement.N_valuses_in_integration_pionts[j][jH];
                    }
                }

                double Hmatrix[4][4];
                double Cmatrix[4][4];
                double k = elements[i].k;


                for (int iH = 0; iH < 4; iH++) {
                    for (int jH = 0; jH < 4; jH++) {
                        Hmatrix[iH][jH] = ((dNdx[iH][jH] + dNdy[iH][jH]) * detJ * k);
                        elements[i].H[iH][jH] += Hmatrix[iH][jH];

                        Cmatrix[iH][jH] = (newElement.sides->weight[0] * elements[i].c * elements[i].ro * c_matrix_temp[iH][jH] * detJ);
                        elements[i].C[iH][jH] += Cmatrix[iH][jH];
                        //cout << elements[i].H[iH][jH]<<" "<<elements[i].c << " " << elements[i].ro << " " << c_matrix_temp[iH][jH] << " "<< detJ<< endl;
                    }
                }


                /*cout << "\nH matrix in " << j + 1 << " integration point\n";
                for (int iH = 0; iH < 4; iH++) {
                    for (int jH = 0; jH < 4; jH++) {
                        cout << Hmatrix[iH][jH] << "  ";
                    }
                     cout << "\n";
                }
                 cout << "===========================" << endl;
                 */
           

            }
        }

        //******* HBC *******
        for (int i = 0; i < this->nE; i++) {

            double alpha = elements[i].alpha;

            for (int vector_iterator = 0; vector_iterator < 4; vector_iterator++) {
                elements[i].p_vector[vector_iterator] = 0;
            }


            for (int i1 = 0; i1 < 4; i1++) {
                if (nodes[elements[i].id[i1] - 1].bcFlag == 1 && nodes[elements[i].id[(i1 + 1) % 4] - 1].bcFlag == 1) {
                    double Pc1[4] = { 0,0,0,0 };
                    double Pc2[4] = { 0,0,0,0 };


                    for (int l = 0; l < 4; l++) {
                        Pc1[l] = newElement.sides[i1].N[0][l]; //wartorsci f. ksztaltu w piierwszym pkt calkowania
                        Pc2[l] = newElement.sides[i1].N[1][l]; //wartorsci f. ksztaltu w drugim pkt calkowania
                    }

                    double detJ = sqrt(pow(nodes[elements[i].id[i1] - 1].x - nodes[elements[i].id[(i1 + 1) % 4] - 1].x, 2) + pow(nodes[elements[i].id[i1] - 1].y - nodes[elements[i].id[(i1 + 1) % 4] - 1].y, 2)) / 2.0; //Jacobian = Dlugosc bokku/2, dlugosc boku liczona z wekotra, potem dlugosc wektora

                    for (int hbcRow = 0; hbcRow < 4; hbcRow++) {
                        elements[i].p_vector[hbcRow] += newElement.sides[i1].weight[0] * Pc1[hbcRow] * elements[i].temperature * elements[i].alpha * detJ;
                        elements[i].p_vector[hbcRow] += newElement.sides[i1].weight[1] * Pc2[hbcRow] * elements[i].temperature * elements[i].alpha * detJ;
                        for (int hbcColumn = 0; hbcColumn < 4; hbcColumn++) {
                            elements[i].hbc[hbcRow][hbcColumn] += newElement.sides[i1].weight[0] * Pc1[hbcRow] * Pc1[hbcColumn] * detJ * elements[i].alpha;
                            elements[i].hbc[hbcRow][hbcColumn] += newElement.sides[i1].weight[1] * Pc2[hbcRow] * Pc2[hbcColumn] * detJ * elements[i].alpha;
                        }
                    }

                }

            }

        }


    }
    else if (nPc == 9) {
        Element9_2D newElement = Element9_2D();
        for (int i = 0; i < this->nE; i++) { //i -> number of element 
            for (int j = 0; j < nPc; j++) { //j -> number of integration point in one element
                jacobianStructure = newElement.jacobian(i, j, newElement, nodes, elements, newJacobian);
                detJ = (jacobianStructure.jacobian[0][0] * jacobianStructure.jacobian[1][1]) - (jacobianStructure.jacobian[1][0] * jacobianStructure.jacobian[0][1]);
                detJinv = 1.0 / detJ;
                for (int k = 0; k < 2; k++) {
                    for (int l = 0; l < 2; l++) {
                        multipliedMatrix[k][l] = jacobianStructure.jacobian_inv[k][l] * detJinv;

                    }
                }


                //Print jacobian in each element and integration point
                /*cout << "Element: " << i + 1 << " Integration Point: " << j + 1 << endl;
                cout << "\nJacobian: " << endl;
                for (int k = 0; k < 2; k++) {
                    for (int l = 0; l < 2; l++) {
                        cout << jacobianStructure.jacobian[k][l] << "\t";
                    }
                    cout << "\n";
                }
                cout << "\n";

                cout << "detJ = " << detJ << endl;
                cout << "1/detJ = " << detJinv << endl << endl;

                cout << "Multiplied Matrix: " << endl;
                for (int k = 0; k < 2; k++) {
                    for (int l = 0; l < 2; l++) {
                               cout << multipliedMatrix[k][l] << "\t";
                    }
                        cout << "\n";
                }
                cout << "===========================" << endl;*/



                //================================ calculate dNi/dNx, dNi/dy ===============================
                double dNdx[4][4]; //store value of dNi/dNx and its transposition
                double dNdy[4][4]; //store value of dNi/dy and its transposition
                double c_matrix_temp[4][4]; //store value of N_valuses_in_integration_pionts and its transposition


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
                        c_matrix_temp[iH][jH] = newElement.N_valuses_in_integration_pionts[j][iH] * newElement.N_valuses_in_integration_pionts[j][jH];
                    }
                }
                double Hmatrix[4][4];
                double Cmatrix[4][4];
                double k = elements[i].k;


                for (int iH = 0; iH < 4; iH++) {
                    for (int jH = 0; jH < 4; jH++) {
                        Hmatrix[iH][jH] = weight * ((dNdx[iH][jH] + dNdy[iH][jH]) * detJ * k);
                        elements[i].H[iH][jH] += Hmatrix[iH][jH];

                        Cmatrix[iH][jH] = (weight * elements[i].c * elements[i].ro * c_matrix_temp[iH][jH] * detJ);
                        elements[i].C[iH][jH] += Cmatrix[iH][jH];
                    }
                }

                // PRINT H MATRIX IN EACH INTEGRATION POINT
                /*cout << "\nH matrix in " << j + 1 << " integration point\n";
                for (int iH = 0; iH < 4; iH++) {
                    for (int jH = 0; jH < 4; jH++) {
                        Hmatrix[iH][jH] = weight * ((dNdx[iH][jH] + dNdy[iH][jH]) * detJ * k);
                        cout << Hmatrix[iH][jH] << "  ";
                        elements[i].H[iH][jH] += Hmatrix[iH][jH];
                    }
                    cout << "\n";
                }*/

            }
        }

        // *** HBC ***
        for (int i = 0; i < this->nE; i++) {

            double alpha = elements[i].alpha;
            for (int vextor_iterator = 0; vextor_iterator < 4; vextor_iterator++) {
                elements[i].p_vector[vextor_iterator] = 0;
            }

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

                    double detJ = sqrt(pow(nodes[elements[i].id[i1] - 1].x - nodes[elements[i].id[(i1 + 1) % 4] - 1].x, 2) + pow(nodes[elements[i].id[i1] - 1].y - nodes[elements[i].id[(i1 + 1) % 4] - 1].y, 2)) / 2.0; //Jacobian = Dlugosc bokku/2, dlugosc boku liczona z wekotra, potem dlugosc wektora

                    for (int hbcRow = 0; hbcRow < 4; hbcRow++) {
                        elements[i].p_vector[hbcRow] += newElement.sides[i1].weight[0] * Pc1[hbcRow] * elements[i].temperature * elements[i].alpha * detJ;
                        elements[i].p_vector[hbcRow] += newElement.sides[i1].weight[1] * Pc2[hbcRow] * elements[i].temperature * elements[i].alpha * detJ;                      
                        elements[i].p_vector[hbcRow] += newElement.sides[i1].weight[2] * Pc3[hbcRow] * elements[i].temperature * elements[i].alpha * detJ;
                        for (int hbcColumn = 0; hbcColumn < 4; hbcColumn++) {
                            elements[i].hbc[hbcRow][hbcColumn] += newElement.sides->weight[0] * Pc1[hbcRow] * Pc1[hbcColumn] * detJ * elements[i].alpha;
                            elements[i].hbc[hbcRow][hbcColumn] += newElement.sides->weight[1] * Pc2[hbcRow] * Pc2[hbcColumn] * detJ * elements[i].alpha;
                            elements[i].hbc[hbcRow][hbcColumn] += newElement.sides->weight[2] * Pc3[hbcRow] * Pc3[hbcColumn] * detJ * elements[i].alpha;
                        }
                    }

                }

            }

        }
    }
    agregation();
    pVectorAgregation();

}





std::vector<double> Grid::solution_t(double tau, std::vector<double> init_temp){
    int i, j, k;
    double m, s;

    std::vector<double> min_max_result;

    std::vector<std::vector<double>> array = HGlobalne;

    //array.push_back(PGlobal);
    int n = HGlobalne.size();
    T_vector_of_solutions.resize(n, 0);

    std::vector<double> initialTempVector;
    initialTempVector.resize(n, 0);

    for (int i = 0; i<HGlobalne.size(); i++) {
        for (int j = 0; j<HGlobalne.size(); j++) {
            array[i][j] += CGlobal[i][j] / tau;
            initialTempVector[i] += CGlobal[i][j] / tau * init_temp[j];
            
            //cout << array[i][j] << "\t";
        }
       // cout << initialTempVector[i] << "\n";
    }


    for (int i = 0; i < nN; i++) {
        array[i].push_back(PGlobal[i]+initialTempVector[i]);
        //cout << PGlobal[i] + initialTempVector[i] << endl;
    }

    // eliminacja wspó³czynników
    for (i = 0; i < n - 1; i++){
        for (j = i + 1; j < n; j++){
            if (array[i][i] == 0.0) {
                return min_max_result;
            }
            m = -array[j][i] / array[i][i];
            for (k = i + 1; k <= n; k++)
                array[j][k] += m * array[i][k];
        }
    }

    // wyliczanie niewiadomych

    for (i = n - 1; i >= 0; i--){
        s = array[i][n];
        for (j = n - 1; j >= i + 1; j--)
            s -= array[i][j] * T_vector_of_solutions[j];
        if (array[i][i] == 0.0) { 
            cout << "Operacja rozwiazywania ukladu rownan nie powiodla sie";
            return min_max_result;
        }
        T_vector_of_solutions[i] = s / array[i][i];
    }




    //print_T_vector_of_solutions();
    return T_vector_of_solutions;
}


void Grid::print_T_vector_of_solutions() {
    cout << "\nT vector of solutions: " << endl;
    for (int i = 0; i < nN ; i++) {
        cout << T_vector_of_solutions[i]<<endl;
    }
}

void Grid::calculate_T_vector_in_iterations(double tau, std::vector<double> iteration_init_temp, int number_of_iterations) {
    this->T_vector_of_solutions = iteration_init_temp;
    for (int i = 0; i < number_of_iterations; i++) {
        std::vector<double> resultin_iteration;
        resultin_iteration = solution_t(tau, T_vector_of_solutions);
        cout << "\nIteration " << i + 1;
        print_T_vector_of_solutions();
        double min = resultin_iteration[0];
        double max = resultin_iteration[0];
        for (int i = 1; i < T_vector_of_solutions.size(); i++) {
            if (resultin_iteration[i] < min) {
                min = resultin_iteration[i];
            }
            if (resultin_iteration[i] > max) {
                max = resultin_iteration[i];
            }
        }
        cout << "MIN = " << min << "\tMAX = " << max <<endl;

    }
}