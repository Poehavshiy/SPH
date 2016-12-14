//
// Created by nikita on 13.05.16.
//

#ifndef SPHSM6_PARTICILE_H
#define SPHSM6_PARTICILE_H

#include "W_functions.h"
//
//QTшные либы
#include <QMainWindow>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QGraphicsEllipseItem>
#include <QMouseEvent>
#include <QPointF>
#include <QPushButton>
#include <QVBoxLayout>
#include <QMenu>
#include <QtGui>
#include <QMenuBar>
#include <QToolButton>
#include <QtWidgets>
//
using namespace std;

struct Point {
    Point() {
        x = 0;
        y = 0;
    }

    Point(double X, double Y) {
        x = X;
        y = Y;
    }

    double x;
    double y;

    friend std::istream &operator<<(ostream &is, Point &income) {
        is << income.x << ' ' << income.y;
    }

    friend bool operator==(const Point &left, const Point &right) {
        if ((left.x == right.x) && (left.y == right.y)) return true;
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
    double molar_mass = 0.029;
    //
    double density;
    double next_density;
    double pressure;
    double energy;
    //
    double vx;
    double vy;
    //
    double p0;
    //производные
    double density_dir;
    double energy_dir;
    double vx_dir;
    double vy_dir;

    static constexpr double R_costil = 8.31;

    static constexpr double k_costil = 1.4;

    static constexpr double molar = 0.029;

    void delta_coor(const Particle &left, double &deltaX, double &deltaY);

    void delta_velocity(const Particle &left, double &delta_Vx, double &delta_Vy);

//искусскственная вязкость
    double two_part_art_visc(const Particle &a, const Particle &b);

public:

    Particle();

    //
    Particle(bool status, double p = 0,
             double P = 0, double E = 0, double Vx = 0, double Vy = 0, double M = 0);

    void calculate_derivatives(const Particle &left, double &h) ;

    //get functions
    double X() const {
        return pos.x;
    }

    //
    double Y() const {
        return pos.y;
    }

    //
    double p() const {
        return density;
    }

    //
    double P() const {
        return pressure;
    }

    double E() const {
        return energy;
    }

    double V() const {
        return hypot(vx, vy);
    }

    //
    double Vx() const {
        return vx;
    }

    //
    double Vy() const {
        return vy;
    }

    //
    double M() const {
        return mass;
    }

    //
    double T() const {
        return pressure * molar_mass / (density * R_costil);
    }

    //
    double C() const {
        double C = sqrt(k_costil * (k_costil - 1) * energy);
    }

    //
    double h(double &h0) const {
        return 1.3 * pow(mass / density, 0.3333333);
    }

    //
    Point position() const {
        return pos;
    }

    void set_pos(double X, double Y) {
        pos.x = X;
        pos.y = Y;
    }

    //
    void set_pos(Point &income) {
        pos.x = income.x;
        pos.y = income.y;
    }

    //
    void set_vx(double income) {
        vx = income;
    }

    //
    void set_vy(double income) {
        vy = income;
    }

    //
    void set_p_dir(double inc) {
        density_dir = inc;
    }

    //
    void set_vx_dir(double inc) {
        vx_dir = inc;
    }
    //

    void set_vy_dir(double inc) {
        vy_dir = inc;
    }
    //
    void set_e_dir(double inc) {
        energy_dir = inc;
    }

    //
    void set_pressure(double &data) {
        pressure = data;
    }

    //
    void set_density(double &data) {
        density = data;
    }

    //
    void set_energy(double &data) {
        energy = data;
    }

    //
    void set_mass(double &data) {
        mass = data;
    }

    //пересчет значений
    void ronge_kutt(double dt);

    //
    void set_from(Particle &data) {
        mass = data.mass;
        //position
        pos = data.pos;
        //
        density = data.density;
        pressure = data.pressure;
        energy = data.energy;
        //
        vx = data.vx;
        vy = data.vy;
    }

    //
    bool boundary_status() {
        return isBoundary;
    }
    //
};


#endif //SPHSM6_PARTICILE_H
