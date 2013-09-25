//delogit.cpp
#include "delogit.h"
#include <math.h> 

double delogit(double x){
    x = exp(x);
    x /= (1+x);
    return x;
}
