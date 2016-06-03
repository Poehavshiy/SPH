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

    friend std::istream& operator<<(ostream& is, Point& income) {
        is<<income.x<<' '<<income.y;
    }

    friend bool operator==(const Point& left, const Point& right) {
        if((left.x==right.x) && (left.y==right.y)) return true;
        else return false;
    }
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
    //
   /* Particle(const Particle& left) {
        bool isBoundary=left.isBoundary;
        double mass=left.mass;
        Point pos=left.pos;
        double density=left.density;
        double pressure=left.pressure;
        double energy=left.energy;
        double vx=left.vx;
        double vy=left.vy;
    }*/

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
    Point position() {
        return pos;
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
    void set_vx(double income) {
        vx=income;
    }
    //
    void set_vy(double income) {
        vy=income;
    }
    //
    void set_pressure(double& data) {
        pressure=data;
    }
    //
    void set_density(double& data) {
        density=data;
    }
    //
    void set_energy(double& data) {
        energy=data;
    }
    //
    void set_mass(double& data) {
        mass=data;
    }
    //
    void set_from(Particle& data) {
        mass=data.mass;
        //position
        pos=data.pos;
        //
        density=data.density;
        pressure=data.pressure;
        energy=data.energy;
        //
        vx=data.vx;
        vy=data.vy;
    }
    //
    bool boundary_status() {
        return isBoundary;
    }
    //
    friend std::istream& operator<<(ostream& is, Particle& income) {
        is<<income.X()<<' '<<income.Y()<<' '<<income.P()<<endl;
    }

};

typedef vector<vector<Particle*>*> PartPointers;
typedef vector<vector<Particle>*> PartPointers_add;


#endif //SPHSM6_PARTICILE_H
