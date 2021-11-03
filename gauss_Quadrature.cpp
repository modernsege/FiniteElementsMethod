#include <iostream>
#include <vector>
#include <math.h>  
#include "gauss_Quadrature.h"

double fun(double x) {
    return(5 * x * x + 3 * x + 6);
}
double fun2var(double x, double y) {
    return(5 * pow(x, 2) * pow(y, 2) + 3 * x * y + 6);
}


N1::N1() {
    nodeXk[0] = -1 / sqrt(3);
    nodeXk[1] = 1 / sqrt(3);
    wages[0] = 1;
    wages[1] = 1;
}



N2::N2() {
    nodeXk[0] = -(sqrt(3.0 / 5.0));
    nodeXk[1] = 0;
    nodeXk[2] = sqrt(3.0 / 5.0);
    wages[0] = 5.0 / 9.0;
    wages[1] = 8.0 / 9.0;
    wages[2] = 5.0 / 9.0;
}


double gaussQuadrature1D(int N) {
    double result = 0;

    if (N == 1) {
        double nodeXk[2] = { -1 / sqrt(3), 1 / sqrt(3) };
        double wages[2] = { 1, 1 };
        for (int i = 0; i < 2; i++) {
            result += wages[i] * fun(nodeXk[i]);
        }
    }
    else if (N == 2) {
        double nodeXk[3] = { -(sqrt(3.0 / 5.0)), 0, sqrt(3.0 / 5.0) };
        double wages[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
        for (int i = 0; i < 3; i++) {
            result += wages[i] * fun(nodeXk[i]);
        }
    }


    return result;
}



double gaussQuadrature2D(int N) {
    double result = 0;

    if (N == 1) {
        double nodeXk[2] = { -1 / sqrt(3), 1 / sqrt(3) };
        double wages[2] = { 1, 1 };
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                result += wages[i] * wages[j] * fun2var(nodeXk[j], nodeXk[i]);
            }
        }
    }
    else if (N == 2) {
        double nodeXk[3] = { -(sqrt(3.0 / 5.0)), 0, sqrt(3.0 / 5.0) };
        double wages[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                result += wages[i] * wages[j] * fun2var(nodeXk[j], nodeXk[i]);
            }
        }
    }


    return result;
}