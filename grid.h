#include <iostream>
#include <vector>
#include <math.h>  


struct Node {
    float x, y;
    short int bcFlag; //deteminate if Node shoud have boundy condition point
    Node(float x, float y, short int bcFlag);
    void printCoordinates();
    void printBcFlag();
};

struct Element {
    int id[4]; //array of node IDs
    double H[4][4]; //H matrix for each element
    double k; //stala dla danego materialu (przy obliczaniu macierzy H)
    double alpha; //wspolczynnik wyminay ciepla
    double hbc[4][4]; //hbc matrix for each element
    double detJ;
    double detJinv;
    double multipliedMatrix[2][2];
    double p_vector[4];
    double temperature;

    Element(int id1, int id2, int id3, int id4, double alpha, double k, double temprature);
    void printIdOfElement();
};



class Grid {
public:
    double H; //height
    double B; //width
    int nH; //number of elements on height
    int nB; //number of elements on width
    int nN; //number of nodes
    int nE; //number of elements
    double deltaH; //value of the height (y coord) change between nodes
    double deltaB; //value of the width (x coord) change between nodes
    std::vector <Node> nodes; //will contain all Nodes of Grid
    std::vector <Element> elements; //will contain all Elements instances of Grid
    std::vector<std::vector<double>> HGlobalne;
    std::vector<double> PGlobal;
    std::vector<double> T_vector_of_solutions;
    Grid(double H, double B, int nH, int nB, int nPc, double alpha, double k, double temprature);
    void printElements();
    void printElementsAndNodesCoordinates();
    void printNodes();
    void print_H_matrix();
    void print_HBC_matrix();
    void print_P_vector();
    void print_H_with_HBC_Global();
    void print_P_Global_vector();
    void agregation();
    void pVectorAgregation();
    void solution_t();
    void print_T_vector_of_solutions();
};