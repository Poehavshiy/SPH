//
// Created by nikita on 13.05.16.
//

#ifndef SPHSM6_PARTICILE_H
#define SPHSM6_PARTICILE_H
//
#include <cmath>
#include <iostream>
#include <vector>
#include<map>

using namespace std;

//
struct Point {
    Point() {
        x=0;
        y=0;
    }
    Point(double X, double Y) {
        x = X;
        y = Y;
    }

    double x;
    double y;
};

class Particle {
    bool isBoundary;
    //1 if its boundary cell
    //data
    double mass;
    //position
    Point pos;
    //
    double density;
    double pressure;
    double energy;
    //
    double vx;
    double vy;
    //
public:
    Particle();

    //
    Particle(bool status, double p = 0,
             double P = 0, double E = 0, double Vx = 0, double Vy = 0, double M = 0);

    //get functions
    double X() {
        return pos.x;
    }

    //
    double Y() {
        return pos.y;
    }

    //
    double p() {
        return density;
    }

    //
    double P() {
        return pressure;
    }

    double E() {
        return energy;
    }

    //
    double Vx() {
        return vx;
    }

    //
    double Vy() {
        return vy;
    }

    //
    double M() {
        return mass;
    }

    //
    void set_pos(double X, double Y) {
        pos.x = X;
        pos.y = Y;
    }
    //
    void set_pos(Point& income) {
        pos.x = income.x;
        pos.y = income.y;
    }
    //
    void set_pressure(double& data) {
        pressure=data;
    }
    //
    friend std::istream& operator<<(ostream& is, Particle& income) {
        is<<income.X()<<' '<<income.Y()<<' '<<income.P()<<endl;
    }

};


#endif //SPHSM6_PARTICILE_H
