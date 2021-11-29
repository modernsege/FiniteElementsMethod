#include <iostream>
#include <vector>
#include <math.h>
#include "element_2D.h"
#include "gauss_Quadrature.h"
#include"grid.h"

using namespace std;



int main()
{
    Grid newGrid = Grid(0.1, 0.1, 4, 4, 4, 300.0, 25.0, 1200.0); //double H, double B, int nH, int nB, int nPc, double alpha, double k, double temprature 
    //0.2, 0.1, 5, 4      0.025, 0.025, 2, 2
    //newGrid.print_H_matrix();
    //newGrid.print_HBC_matrix();
    //newGrid.print_P_vector();
    newGrid.print_H_with_HBC_Global();
    newGrid.print_P_Global_vector();
    newGrid.print_T_vector_of_solutions();

    //newGrid.printElements();
    //newGrid.printNodes();
    //newGrid.printElementsAndNodesCoordinates();

    return 0;
}

