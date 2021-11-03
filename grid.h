#include <iostream>
#include <vector>
#include <math.h>  

struct Node {
    float x, y;
    Node(float x, float y);
    void printCoordinates();
};

struct Element {
    int id[4]; //array of node IDs

    Element(int id1, int id2, int id3, int id4);
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
    Grid(double H, double B, int nH, int nB);
    void printElements();
    void printElementsAndNodesCoordinates();
    void printNodes();
    void calculateJacobianForGrid(int);
};