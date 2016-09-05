//
// Created by nikita on 13.05.16.
//

#include "Particile.h"

Particle::Particle():pos(0, 0) {

    isBoundary = 0;
    //
    density = 0;
    pressure = 0;
    energy = 0;
    //
    vx = 0;
    vy = 0;
    //
    mass = 0;
};

//
Particle::Particle(bool status, double p,
                   double P, double E, double Vx, double Vy, double M ) {

    isBoundary = status;
    density = p;
    pressure = P;
    energy = E;
    //
    vx = Vx;
    vy = Vy;
    //
    mass = M;

    p0 = p;
};
//
