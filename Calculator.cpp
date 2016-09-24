//
// Created by nikita on 30.05.16.
//
//размер ячейки 50

#include "Flow.h"
#include "Calculator.h"

//для трубы делай r0 = h/5 (максимум 6) D = 2maxV*maxV
QGraphicsScene *scene_debug;

double Calculator::current_time = 0;
double Calculator::h;//выбирается из соображения того, что рядом должна быть 21 частица
// размер ячейки равен 2h
double Calculator::R = 8.31;
double Calculator::k = 1.4;
double Calculator::M = 0.029;
double Calculator::viscous = 0.00023;
double Calculator::alpha = 0.9;
double Calculator::beta = 0.9;
//для подсчета силы от гарничных частиц
double Calculator::r0;
double Calculator::D = 0;//равен квадрату наибольшей скорости
double Calculator::n1 = 12;
double Calculator::n2 = 6;
int Calculator::iteration = 0;

//typedef vector<vector<Particle*>*> PartPointers;
//typedef vector<vector<Particle>*> PartPointers_add;

void Calculator::draw_debug(vector<Cell *> &for_sum_calculating) {
    double rad = 2;
    const vector<Particle *> *real;
    const vector<Particle> *shadow;
    const vector<Particle *> *bound;
    for (int i = 0; i < for_sum_calculating.size(); ++i) {
        real = for_sum_calculating[i]->get_real();
        for (int j = 0; j < real->size(); ++j) {
            const Particle *b = real->operator[](j);
            scene_debug->addEllipse(b->X() - rad, b->Y() - rad, rad * 2.0, rad * 2.0,
                                    QPen(Qt::green), QBrush(Qt::SolidPattern));
        }
        shadow = for_sum_calculating[i]->get_shadow();
        for (int j = 0; j < shadow->size(); ++j) {
            const Particle *b = &shadow->operator[](j);
            scene_debug->addEllipse(b->X() - rad, b->Y() - rad, rad * 2.0, rad * 2.0,
                                    QPen(Qt::blue), QBrush(Qt::SolidPattern));
        }
        //теперь, если считается скорость, посчитаем силу отталкивания
        bound = for_sum_calculating[i]->get_bound();
        for (int j = 0; j < bound->size(); ++j) {
            const Particle *b = bound->operator[](j);
            scene_debug->addEllipse(b->X() - rad, b->Y() - rad, rad * 2.0, rad * 2.0,
                                    QPen(Qt::black), QBrush(Qt::SolidPattern));
        }
    }
}

bool Calculator::to_boundary(Particle &a, Cell &target) {
    Particle *A = target.get_bound()->operator[](0);
    Particle *B = target.get_bound()->back();
    double m = (B->Y() - A->Y()) / (B->X() - A->X());
    double k = m * A->X() + A->Y();
    double x_N = (a.X() + m * a.Y() - m * k) / (pow(m, 2) + 1);
    double y_N = m * (a.X() + m * a.Y() - m * k) / (pow(m, 2) + 1) + k;
    double normal_x = x_N - a.X();
    double normal_y = y_N - a.Y();
    double sc_product = normal_x * a.Vx() + normal_y * a.Vy();
    if (sc_product > 0)
        return true;
    else
        return false;

}
//
bool Calculator::to_boundary(Particle &a, Point &A, Point &B) {
    double m = (B.y - A.y) / (B.x - A.x);
    double k = m * A.x + A.y;
    double x_N = (a.X() + m * a.Y() - m * k) / (pow(m, 2) + 1);
    double y_N = m * (a.X() + m * a.Y() - m * k) / (pow(m, 2) + 1) + k;
    double normal_x = x_N - a.X();
    double normal_y = y_N - a.Y();
    double sc_product = normal_x * a.Vx() + normal_y * a.Vy();
    if (sc_product > 0)
        return true;
    else
        return false;
}

void Calculator::delta_coor(const Particle &a, const Particle &b, double &deltaX, double &deltaY) {
    deltaX = a.X() - b.X();
    deltaY = a.Y() - b.Y();
}

void Calculator::delta_velocity(const Particle &a, const Particle &b, double &delta_Vx, double &delta_Vy) {
    delta_Vx = a.Vx() - b.Vx();
    delta_Vy = a.Vy() - b.Vy();
}

//Производная W по координате direction(0-x 1-y)
//для вязкости
double Calculator::two_part_E(const Particle &a, const Particle &b, bool direct) {
    //direct 0 -> вычисляем Exy, если direct 1 -> вычисляем Eyx
    double mj_pj = b.M() / b.p();

    //произведение Yскорости на Xпроизводную
    double Vj_y_Wi_x;
    //произведение Xскорости на Yпроизводную
    double Vj_x_Wi_y;
    double deltaX;
    double deltaY;
    delta_coor(a, b, deltaX, deltaY);

    if (direct == 0) {
        Vj_y_Wi_x = b.Vy() * grad_w_test(deltaX, deltaY, 0, a.h());
        Vj_x_Wi_y = b.Vx() * grad_w_test(deltaX, deltaY, 1, a.h());
    } else {
        Vj_y_Wi_x = b.Vx() * grad_w_test(deltaX, deltaY, 1, a.h());
        Vj_x_Wi_y = b.Vy() * grad_w_test(deltaX, deltaY, 0, a.h());
    }
    double result = mj_pj * (Vj_y_Wi_x + Vj_x_Wi_y);
    return result;
}

//искусскственная вязкость
double Calculator::two_part_art_visc(const Particle &a, const Particle &b) {
    double delta_vx, delta_vy;
    double deltaX, deltaY;
    delta_coor(a, b, deltaX, deltaY);
    delta_velocity(a, b, delta_vx, delta_vy);
    if( (delta_vx*deltaX + delta_vy*deltaY) < 0) return 0;
    double r = hypot(deltaX, deltaY);
    double V = hypot(deltaX, deltaY);


    double phi = h * r * V / (pow(r, 2) + pow(0.1 * h, 2));
    double result;
    result = 2 * (-alpha * 0.5 * (a.C() + b.C()) * phi + beta * pow(phi, 2)) / (a.p() + b.p());
    return result;
}

//искусственная теплота
double Calculator::two_part_art_heat(const Particle &a, const Particle &b, bool direct = 0) {
    /* double deltaX = b.X() - a.X();
     double deltaY = b.Y() - a.Y();
     double divVa = b.M() *(b.Vx()*grad_w_test(deltaX, deltaY, 0)+b.Vy()*grad_w_test(deltaX, deltaY, 1))/b.p();
     double qi = alpha*h*a.p()*a.C()*abs(divVa)+beta*h*a.p()*pow(divVa,2);

     double divVb = a.M() *(b.Vx()*grad_w_test(deltaX, deltaY, 0)+b.Vy()*grad_w_test(deltaX, deltaY, 1))/b.p();
     double qj = alpha*h*b.p()*b.C()*abs(divV)+beta*h*b.p()*pow(divV,2);*/
    return 0;

}

double Calculator::grad_w_test(double x, double y, bool direction, double h_inc) {
    //double a = 15 / (7*3.14 * pow(h, 2));
    double a = 1 / (3.14 * pow(h, 4));
    double r = hypot(x, y) / h;
    double minR = 0.001;
    if (r < 0 || r > 2 || r < minR) return 0;
    double result;
    double factor;
    if (direction == 0) {
        factor = x / hypot(x, y);
    } else {
        factor = y / hypot(x, y);
    }
    //  a = 15/(7*3.14*h*h);
    // result = -9*r/4 + 19*pow(r,2)/8 - 5*pow(r,3)/8;
    if (r > 0 && r <= 1) result = 9 * pow(r, 2) / 4 - 3 * r;
    else if (r > 1 && r <= 2) result = -(3 / 4) * pow(2 - r, 2);
   // result = a*(-18*r / 8 + 57*pow(r,2)/24 - 20*pow(r,3)/32);

    return result * a * factor;//factor instead of 1
}

//
double Calculator::two_part_p(const Particle &a, const Particle &b, bool direct) {
    double Mj = b.M();
    double delta_vx, delta_vy;
    double deltaX, deltaY;
    delta_coor(a, b, deltaX, deltaY);
    delta_velocity(a, b, delta_vx, delta_vy);
    //скалярное произведение скоростей и градиента
    double VijWij = delta_vx * grad_w_test(deltaX, deltaY, 0, h) + delta_vy * grad_w_test(deltaX, deltaY, 1, h);
    double res = Mj * VijWij;
    return res;
}

//
double Calculator::two_part_v(const Particle &a, const Particle &b, bool direct) {
    double Mj = b.M();
    double deltaX;
    double deltaY;
    delta_coor(a, b, deltaX, deltaY);
    double Wij = grad_w_test(deltaX, deltaY, direct, h);
    double H = two_part_art_visc(a, b);
    double brackets1 = -(a.P() / pow(a.p(), 2) +
                         b.P() / pow(b.p(), 2) + H);
    //
    double res = Mj * brackets1 * Wij;//1е слагаемоеж
    //считаем 2е
    if (EULER == false) {
        double brackets2 = viscous * two_part_E(a, b, direct) / pow(a.p(), 2) +
                           //второе слагаемое  Ej a b поменяли
                           viscous * two_part_E(b, a, direct) / pow(b.p(), 2);
        res += Mj * brackets2 * Wij;
        //добавим иск вязкость
        double H = two_part_art_visc(a, b);

        res += H * Wij;
    }

    return res;
}

//
double Calculator::two_part_energy(const Particle &a, const Particle &b, bool direct) {
    double delta_vx, delta_vy;
    double deltaX, deltaY;
    delta_coor(a, b, deltaX, deltaY);
    delta_velocity(a, b, delta_vx, delta_vy);
    //скалярное произведение Vij и Wij
    double VijWij = delta_vx * grad_w_test(deltaX, deltaY, 0, h) +
                    delta_vy * grad_w_test(deltaX, deltaY, 1, h);
    //
    double Mj = b.M();
    double H = two_part_art_visc(a, b);
    double brackets1 = (a.P() / pow(a.p(), 2) +
                        b.P() / pow(b.p(), 2) + H);
    //

    double res = 0.5 * VijWij * Mj * brackets1;
    if (EULER == false) {
        ///?????????????????????????????????
        double brackets2 = (viscous / (2 * a.p())) * two_part_E(a, b, 0) * two_part_E(a, b, 1);
        res += brackets2;
    }
    //??????????????????????????????????????????
    return res;

}

//высчитывает часть производной скорости от действия граничной частицы
double Calculator::two_part_bforse(const Particle &a, const Particle &b, bool direct = 0) {
    double deltaX, deltaY;
    delta_coor(a, b, deltaX, deltaY);
    double r = hypot(deltaX, deltaY);
    double ratio = r0 / r;
    if (ratio > 1) return 0;
    else {
        double x = deltaX;
        if (direct == 1) x = deltaY;
        double res = D * (pow(ratio, n1) - pow(ratio, n2)) * x / r;
        res = res / a.M();
        return -res;
    }
}

/*
 считаем взаимодействие частицы a cо всеми частицами в HomeCell и  for_sum_calculating
 */
double Calculator::calc_particle_dir(function<double(const Particle &, const Particle &, bool)> &derivative_specific,
                                     Particle &a,
                                     vector<Cell *> &for_sum_calculating, bool direction, DERIVATIVES &what_dir) {

    //debugin
    if (debug_mode == true && &a == debug_draw_part) {
        double rad = 5;
        draw_debug(for_sum_calculating);
        scene_debug->addEllipse(debug_draw_part->X() - rad, debug_draw_part->Y() - rad, rad * 2.0, rad * 2.0,
                                QPen(Qt::green), QBrush(Qt::SolidPattern));
    }

    double res = 0;
    ////пройдем по реальным
    const vector<Particle *> *real;
    const vector<Particle> *shadow;
    for (int i = 0; i < for_sum_calculating.size(); ++i) {
        real = for_sum_calculating[i]->get_real();
        for (int j = 0; j < real->size(); ++j) {
            const Particle *b = real->operator[](j);
            res += derivative_specific(a, *b, direction);
        }
        //теперь от виртуальных частиц
        // if (what_dir == VX || what_dir == VY) {
        shadow = for_sum_calculating[i]->get_shadow();
        for (int j = 0; j < shadow->size(); ++j) {
            const Particle *b = &shadow->operator[](j);
            res += derivative_specific(a, *b, direction);
        }
        //}
    }
    return res;
}

//вот этот метод дергается в методе calculate, который во Flow
void Calculator::calculate_derivatives() {
    /*
     * для каждой i,j области пространства создаются вектора всех реальных, граничных и теневых
     * частиц с которыми потенциально может провзаимодействовать любая частица в этой i,j области
     */
    for (int i = 0; i < parsing->cells_per_y; ++i) {
        for (int j = 0; j < parsing->cells_per_x; ++j) {
            if (parsing->part_groups[i][j].real_group.size() != 0) {//если ячейка содержит реальные частицы

                vector<Cell *> for_sum_calculating;
                //соберем всех соседей i,j Cell
                for (int k = 0; k < 8; ++k) {
                    Cell *current = parsing->get_cell_clockwise(i, j, k);
                    if (current != nullptr)
                        for_sum_calculating.push_back(current);
                }
                for_sum_calculating.push_back(&parsing->part_groups[i][j]);
                //для частиц каждой каждой группы высчитываются производные и запихиваются в соответствующие
                for (int it = DENSITY; it <= ENERGY; ++it) {
                    DERIVATIVES dir = static_cast<DERIVATIVES >(it);
                    bool direction = 0;
                    if (dir == VY) direction = 1;
                    calc_target_derivatives(derivative_functions[it],
                                            parsing->part_groups[i][j], for_sum_calculating, dir,
                                            direction, i, j);
                }
            }
        }
    }
}

//вот эти 4 функции считают прозводные ДЛЯ ВСЕХ ЧАСТИЦ В target
void Calculator::calc_target_derivatives(function<double(const Particle &, const Particle &, bool)> derivative_specific,
                                         Cell &target, vector<Cell *> &for_sum_calculating,
                                         DERIVATIVES &what_dir, bool direction,
                                         int &row, int &column) {
//посчитаем
    int size = target.real_group.size();
    double res = 0;
    for (int i = 0; i < size; ++i) {
        Particle &target_p = *target.real_group[i];
        res = calc_particle_dir(derivative_specific, target_p, for_sum_calculating, direction, what_dir);

        if ((what_dir == VX || what_dir == VY) && target.is_boundary() == true) {
            double b_forse = boundary_forse(target_p, for_sum_calculating, direction);
            res += b_forse;

        }
        target.real_group[i]->set_dir(res, what_dir);
    }
}


void Calculator::recalculate_parameters() {
    int data_size = parsing->data_ptr->size();
    double cur_maxV = 0;
    double cur_maxC = 0;
    double cur_maxP = 0;
    double cur_minP = 0;
    for (int i = 0; i < data_size; ++i) {
        Particle *particle = &parsing->data_ptr->operator[](i);
        particle->ronge_kutt(dt);

        cur_maxV = particle->V();
        cur_maxC = particle->C();
        cur_maxP = particle->P();
        cur_minP = particle->P();
        if (cur_maxV > max_V) max_V = cur_maxV;
        if (cur_maxC > max_C) max_C = cur_maxC;
        if (cur_maxP > max_P) max_P = cur_maxP;
        if (cur_minP < min_P) min_P = cur_minP;
    }
}

double Calculator::boundary_forse(Particle &a, vector<Cell *> &for_sum_calculating, bool direction) {
    double res = 0;
    int size = for_sum_calculating.size();
    for (int i = 0; i < size; ++i) {
        const vector<Particle *> &bound = *for_sum_calculating[i]->get_bound();
        int boundsize = bound.size();
        for (int j = 0; j < boundsize; ++j) {
            res += two_part_bforse(a, *bound[j], direction);
        }
    }
    return res;
}

void Calculator::calculate() {

    max_V = 0;
    max_C = 0;
    //всякие информативные штуки
    max_P = 0;
    min_P = 100500;

    parsing->clear_symetric_groups();
    parsing->create_symetric_groups();
    calculate_derivatives();
    recalculate_parameters();//пересчитали конечные параметры
    parsing->replace();

    D = max_V * max_V * 3;
    dt = h / (6 * max_C);
    current_time += dt;
    iteration++;
    //
}

//
Calculator::Calculator(SpaceParsing *target) {
    //
    derivative_functions[0] = &(Calculator::two_part_p);
    derivative_functions[1] = &(Calculator::two_part_v);
    derivative_functions[2] = &(Calculator::two_part_v);
    derivative_functions[3] = &(Calculator::two_part_energy);
    //
    parsing = target;
    debug_draw_part = &target->data_ptr->operator[](50);
    V_theoretical = sqrt(4 * k * parsing->data_ptr->operator[](0).E() / (k - 1));
    max_V = 0;
    max_C = 0;
    dt = 0.1;
    //всякие информативные штуки
    max_P = 0;
    min_P = 100500;
    //
    debug_mode = false;
    h = parsing->x_size/2;
    r0 =parsing->x_size/5;
}