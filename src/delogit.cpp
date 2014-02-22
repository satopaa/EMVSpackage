//delogit.cpp
#include "delogit.h"
#include <math.h> 

vec delogit(vec &x){
    x = exp(x);
    x /= (1+x);
    return x;
}
