//
// Created by nikita on 13.05.16.
//

#include "Particile.h"

Particle::Particle() : pos(0, 0) {

    isBoundary = 0;
    //
    density = 0;
    next_density = 0;
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
                   double P, double E, double Vx, double Vy, double M) {

    isBoundary = status;
    density = p;
    next_density = 0;
    pressure = P;
    energy = E;
    full_energy = energy + 0.5 * hypot(vx, vy) * hypot(vx, vy);
    //
    vx = Vx;
    vy = Vy;
    //
    mass = M;

    p0 = p;
};


//У китайцев ВСЕГДА из iй вычитают  jю
void Particle::delta_coor(const Particle &left, double &deltaX, double &deltaY) {
    deltaX = X() - left.X();
    deltaY = Y() - left.Y();
}

void Particle::delta_velocity(const Particle &left, double &delta_Vx, double &delta_Vy) {
    delta_Vx = Vx() - left.Vx();
    delta_Vy = Vy() - left.Vy();
}

double Particle::two_part_art_visc(const Particle &a, const Particle &b) {
    /* double delta_vx, delta_vy;
     double deltaX, deltaY;
     delta_coor(a, b, deltaX, deltaY);
     delta_velocity(b, a, delta_vx, delta_vy);
     double scalar_product = delta_vx * deltaX + delta_vy * deltaY;
     if (scalar_product > 0) return 0;
     double r = hypot(deltaX, deltaY);
     double h_ab = 0.5 * (a.h(h0) + b.h(h0));
     //h_ab = h;//опционально
     double phi = h_ab * scalar_product / (pow(r, 2) + pow(artvis_phi * h_ab, 2));
     double c_ab = 0.5 * (a.C() + b.C());
     double p_ab = 0.5 * (a.p() + b.p());
     double result = (-artvis_alpha * c_ab * phi + artvis_beta * pow(phi, 2)) / p_ab;
 */
    return 1;
}

void Particle::calculate_derivatives(const Particle &left, double &h) {

    double mj = left.mass;
    double delta_x = 0;
    double delta_y = 0;
    double delta_Vx = 0;
    double delta_Vy = 0;
    delta_coor(left, delta_x, delta_y);
    delta_velocity(left, delta_Vx, delta_Vy);

    double H = 0;
    double r = hypot(left.X() - X(), left.Y() - Y());


    if (left.isBoundary) {
        double r0 = 1;
        double ratio = r0 / r;
        if (ratio > 1) return;
        else {
            double n1 = 12;
            double n2 = 4;
            double D = this->V() > left.V() ? V() * V() : left.V() * left.V();
            double forse = D * (pow(ratio, n1) - pow(ratio, n2));
            double x_forse = forse * delta_x / r;
            double y_forse = forse * delta_y / r;
            vx_dir += x_forse / mass;
            vy_dir += y_forse / mass;
        }
        return;
    }

    double dW = W_functions::dW_Bykov(r, h);
    //double dW = W_functions::dW_cubic_smoothing_function(r, h);

    double dWx = (delta_x / r) * dW;
    double dWy = (delta_y / r) * dW;
    if (this == &left) {
        dWx = 0;
        dWy = 0;
    }
    double scalar_dW_V = delta_Vx * dWx + delta_Vy * dWy;

    double brackets = mj *
                      (pressure / (density * density) +
                       left.pressure / (left.density * left.density)
                       + H);

    vx_dir += -brackets * dWx;
    vy_dir += -brackets * dWy;

    energy_dir += 0.5 * brackets * scalar_dW_V;

    //
    next_density += left.mass * W_functions::W_Bykov(r, h);
    //next_density += left.mass * W_functions::W_cubic_smoothing_function(r, h);
    //density_dir += mj * scalar_dW_V;
}

void Particle::calculate_derivatives_dis(const Particle &left, double opt_h) {
    if (left.isBoundary) {
        double r0 = 0.1;
        double r = hypot(left.X() - X(), left.Y() - Y());
        if (r0 / r <= 1) {
            double n1 = 12;
            double n2 = 4;
            double D = V() * V();
            D = D ;
            double brackets = -D * ( pow(r0 / r, n1) - pow(r0 / r, n2) ) / (r * r);
            double PDijx = brackets * (left.X() - X());
            double PDijy = brackets * (left.Y() - Y());
            //
            vx_dir += PDijx / mass;
            vy_dir += PDijy / mass;
            if(fabs(D) > 0.00000000000001){
                int a = 1;
            }
            return;
        }
        return;
    }
    double alpha = 1.2;
    double h = alpha * (D() + left.D()) * 0.5;
    if(opt_h > 0){
        h = opt_h;
    }
    double r = hypot(left.X() - X(), left.Y() - Y());
    if (r == 0) return;
    double dW = W_functions::dW_disser(r, h);
    double m_j = left.mass;
    double Ro_j = left.density;
    double P_j = left.pressure;
    double Cl_j = left.C(); //скорость звука
    double Cl_i = C();
    //
    double Rij_norm_x = (left.pos.x - pos.x) / r;
    double Rij_norm_y = (left.pos.y - pos.y) / r;
    double U_Ri = Rij_norm_x * vx + Rij_norm_y * vy;
    double U_Rj = Rij_norm_x * left.vx + Rij_norm_y * left.vy;
    //
    double U_starRij = (U_Rj * Ro_j * Cl_j + U_Ri * density * Cl_i - P_j + pressure) /
                       (Ro_j * Cl_j + density * Cl_i);// 1.20
    double P_starij =
            (P_j * density * Cl_i + pressure * density * Cl_j + Ro_j * Cl_j * density * Cl_i * (U_Rj - U_Ri)) /
            (Ro_j * Cl_j + density * Cl_i); //1.21
    //
    double mn4VandE = 2 * ((m_j * P_starij) / (density * Ro_j)) * dW;

    eps_disser += -2 * m_j * dW * (U_Ri - U_starRij) / density; //1.83
    vx_dir += mn4VandE * Rij_norm_x;
    vy_dir += mn4VandE * Rij_norm_y;
    full_energy_dir += -mn4VandE * U_starRij * 0.5;
}

void Particle::ronge_kutt(double dt) {
    energy += dt * energy_dir;
    //
    pos.x += dt * vx;
    pos.y += dt * vy;
    pressure = energy * (k_costil - 1) * density;
//    //smooth way
//
//    density = next_density;
//    next_density = 0;

    //traditional way
    density += dt * density_dir;

    vx += dt * vx_dir;
    vy += dt * vy_dir;
    //

}

double Particle::disser_integr(double dt) {
    pos.x += dt * vx;
    pos.y += dt * vy;
    pressure = energy * (k_costil - 1) * density;
    double kin_energy = 0.5 * V() * V();
    vx += dt * vx_dir;
    vy += dt * vy_dir;
    double new_kin_energy = 0.5 * V() * V();
    double delta_full_energy = full_energy_dir * dt;
    //
    energy += delta_full_energy - (new_kin_energy - kin_energy);
    full_energy += dt * full_energy_dir;
    //
    density *= (2 - eps_disser * dt) / (2 + eps_disser * dt);
    //
    density_dir = 0;
    energy_dir = 0;
    vx_dir = 0;
    vy_dir = 0;
    full_energy_dir = 0;
    eps_disser = 0;

    //
    double d = D();
    double c = C();
    return 0.7 * d / sqrt((pow(c, 2) + pow(4 * d * eps_disser, 2)));


}