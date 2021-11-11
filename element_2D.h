#include <iostream>
#include <vector>
#include <math.h> 


struct Jacobian {
    double jacobian[2][2];
    double jacobian_inv[2][2];

};

class Grid;
struct Element;
struct Node;


class Element4_2D {
public:
    double dNdKSI[4][4], dNdETA[4][4];
    double eta, ksi;
    Element4_2D();
    void printFor4IntegrationPoints();
    Jacobian jacobian(int,int, Element4_2D, std::vector<Node>, std::vector<Element>, Jacobian);
    
};


class Element9_2D {
public:
    double dNdKSI[9][4], dNdETA[9][4];
    double eta, ksi;
    Element9_2D();
    void printFor9IntegrationPoints();
    Jacobian jacobian(int, int, Element9_2D, std::vector<Node>, std::vector<Element>, Jacobian);


};