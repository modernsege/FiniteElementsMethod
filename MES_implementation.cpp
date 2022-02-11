#include <iostream>
#include <vector>
#include <math.h>
#include "element_2D.h"
#include "gauss_Quadrature.h"
#include"grid.h"
#include <fstream>
#include <string>

using namespace std;



int main()
{

    fstream file;
    file.open("data.txt", ios::in);
    string line;
    vector<string> tab;
    while (getline(file, line))
    {
        line.erase(std::find(line.begin(), line.end(), '#'), line.end());
        tab.push_back(line);
    }

    file.close();
    double H = atof(tab[0].c_str());
    double B = atof(tab[1].c_str());
    int nH = stoi(tab[2]);
    int nB = stoi(tab[3]);
    int nPc = stoi(tab[4]);
    double alpha = atof(tab[5].c_str());
    double k = atof(tab[6].c_str());
    double temprature = atof(tab[7].c_str());
    double specificHeat = atof(tab[8].c_str());
    double density = atof(tab[9].c_str());
    double initial_temperature = atof(tab[10].c_str());
    double simulation_time = atof(tab[11].c_str());
    double simulation_step_time = atof(tab[12].c_str());
    int number_of_iterations = floor(simulation_time / simulation_step_time);



    Grid newGrid = Grid(H, B, nH, nB, nPc, alpha, k, temprature, specificHeat, density); //double H, double B, int nH, int nB, int nPc, double alpha, doublek,doubletemprature double specificHeat, double density
      //0.1, 0.1, 4, 4, 4, 300.0, 25.0, 1200.0, 700.0, 7800.0
      //0.2, 0.1, 5, 4      0.025, 0.025, 2, 2
      //newGrid.print_H_matrix();
      //newGrid.print_C_matrix();
      //newGrid.print_HBC_matrix();
      //newGrid.print_P_vector();
      //newGrid.print_H_with_HBC_Global();
       //newGrid.print_P_Global_vector();
       //newGrid.print_C_Global();
      // newGrid.print_T_vector_of_solutions();
    std::vector<double> init_temp;
    init_temp.resize(newGrid.nN, initial_temperature);
    newGrid.calculate_T_vector_in_iterations(simulation_step_time, init_temp, number_of_iterations);
    


    return 0;
}

