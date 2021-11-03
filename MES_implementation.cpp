#include <iostream>
#include <vector>
#include <math.h>
#include "element_2D.h"
#include "gauss_Quadrature.h"
#include"grid.h"

using namespace std;




int main()
{
    Grid newGrid = Grid(0.025, 0.025, 2, 2);
    newGrid.calculateJacobianForGrid(4);
    //newGrid.printElements();
    //newGrid.printNodes();
    //newGrid.printElementsAndNodesCoordinates();
    //cout << "\n" << gaussQuadrature1D(1);
    //cout << "\n" << gaussQuadrature1D(2);    
    //cout << "\n" << gaussQuadrature2D(1);
    //cout << "\n" << gaussQuadrature2D(2);

    Element4_2D newElement4_2D = Element4_2D();
    newElement4_2D.printFor4IntegrationPoints();

    Element9_2D newElement9_2D = Element9_2D();
    //newElement9_2D.printFor9IntegrationPoints();

    return 0;
}

