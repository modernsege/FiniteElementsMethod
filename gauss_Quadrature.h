#include <iostream>
#include <vector>
#include <math.h>  


double fun(double x);

double fun2var(double x, double y);

struct N1 {
    double nodeXk[2];
    double wages[2];
    N1();
};

struct N2 {
    double nodeXk[3];
    double wages[3];
    N2();
};

double gaussQuadrature1D(int N);

double gaussQuadrature2D(int N);
