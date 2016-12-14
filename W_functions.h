//
// Created by nikita on 13.12.16.
//

#ifndef SPHSM6_W_FUNCTIONS_H
#define SPHSM6_W_FUNCTIONS_H
//
#include <cmath>
#include <iostream>
#include <vector>
#include<map>
#include <assert.h>
#include "algorithm"

class W_functions {
public:
    W_functions() {

    }

    static double dW_Bykov(const double &r, const double &h) {
        double q = fabs(r) / h;
        if (q > 2.0)
            return 0.0;
        double a = -2.0 / (3.0 * h * h);
        double result = 0;

        if (q < 0.66666)
            result = 1;
        else if (q >= 0.66666 && q < 1.0)
            result = 3.0 * q * (4.0 - 3.0 * q) / 4.0;
        else if (q >= 1.0 && q <= 2.0)
            result = 3.0 * (2.0 - q) * (2.0 - q) / 4.0;

        return result * a;
    }
    //
    static double W_Bykov(const double &r, const double &h){
        double q = fabs(r) / h;
        if(q > 2.0)
            return 0.0;

        double a = 2.0 / (3.0 * h);
        double result = 0;

        if(q >= 0 && q <= 1.0)
            result = 0.25 * (4.0 - 6 * q * q + 3 * q * q * q);
        else if(q > 1.0 && q <= 2.0)
            result = 0.25 * (2.0 - q) * (2.0 - q) * (2.0 - q);
        return result * a;
    }

    //
    static double W_cubic_smoothing_function(const double &r, const double &h){
        double q = fabs(r);
        if(q >= 2*h) return 0;
        double b1 = -2.0 / (3.0 * h);
        double b2 = 1.0 / (6.0 * h);
        double a = 1/pow(h,3);
        if( q <= h){
            return  b1 * pow((1 - q/h),3) + b2 * pow((2 - q/h),3);
        }
        if(q > h){
            return  b2 * pow((2 - q/h),3);
        }
    }

    //
    static double dW_cubic_smoothing_function(const double &r, const double &h) {
        if( fabs(h - 0.0075 < 0.00001) ) throw 228;

        double q = fabs(r);
        if(q >= 2*h) return 0;
        double b1 = -2.0 / (3.0 * h);
        double b2 = 1.0 / (6.0 * h);
        double a = 1/pow(h,3);
        if( q <= h){
            return ( pow((h - q),2) * b1 + pow((2*h - q),2) * b2) * a * -3;
        }
        if(q > h){
            return  -3 * pow((2*h - q),2) * b2 * a;
        }
    }
};


#endif //SPHSM6_W_FUNCTIONS_H
