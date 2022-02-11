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
    double c; //cieplo wlascie
    double ro; //gestosc
    double C[4][4]; //C-matrix

    Element(int id1, int id2, int id3, int id4, double alpha, double k, double temprature, double specificHeat, double density);
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
    int nPc; //liczba punktow calkowania
    std::vector <Node> nodes; //will contain all Nodes of Grid
    std::vector <Element> elements; //will contain all Elements instances of Grid
    std::vector<std::vector<double>> HGlobalne;
    std::vector<double> PGlobal;
    std::vector<double> T_vector_of_solutions;
    std::vector<std::vector<double>> CGlobal;
    Grid(double H, double B, int nH, int nB, int nPc, double alpha, double k, double temprature, double specificHeat, double density);
    void printElements();
    void printElementsAndNodesCoordinates();
    void printNodes();
    void calculate_H_matrix(int);
    void calculate_C_matrix(int);
    void calculate_P_vector(int);
    void calculate_HBC_matrix(int);
    void reset_structures();//fill with 0

    void print_H_matrix();
    void print_C_matrix();
    void print_HBC_matrix();
    void print_P_vector();

    void print_H_with_HBC_Global();
    void print_C_Global();
    void print_P_Global_vector();
    void agregation();
    void pVectorAgregation();
    std::vector<double> solution_t(double tau, std::vector<double> init_temp);
    void print_T_vector_of_solutions();
    void calculate_T_vector_in_iterations(double tau, std::vector<double> iteration_init_temp, int nuber_of_iterations);
};