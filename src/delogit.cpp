//delogit.cpp
#include "delogit.h"
#include <math.h> 

double delogit(double x){
    double temp = exp(x);
    temp /= (1+temp);
    return temp;
}