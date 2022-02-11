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

struct coordinates {
    double x;
    double y;
};

struct sideOfElement4_2d {
    coordinates integrationPointsBC[2]; //pkt ca³kowania
    double weight[2]; // wagi
    double N[2][4]; //wart. f. kszta³tu w pkt calkowania (pierwzy wymiar 2, bo na scianie 2 pkt calkowania, czyli 4 wartosci f ksztaltu w pierwszym pkt calkowania i 4 w drugim)
};


class Element4_2D {
public:
    double dNdKSI[4][4], dNdETA[4][4], N_valuses_in_integration_pionts[4][4];
    double eta, ksi;
    Element4_2D();
    void printFor4IntegrationPoints();
    Jacobian jacobian(int, int, Element4_2D, std::vector<Node>, std::vector<Element>, Jacobian);
    sideOfElement4_2d sides[4];
};


struct sideOfElement9_2d {
    coordinates integrationPointsBC[3]; //pkt ca³kowania
    double weight[3]; // wagi
    double N[3][4]; //wart. f. kszta³tu w pkt calkowania
};

class Element9_2D {
public:
    double dNdKSI[9][4], dNdETA[9][4], N_valuses_in_integration_pionts[9][4];
    double eta, ksi;
    Element9_2D();
    void printFor9IntegrationPoints();
    Jacobian jacobian(int, int, Element9_2D, std::vector<Node>, std::vector<Element>, Jacobian);
    sideOfElement9_2d sides[4];

};