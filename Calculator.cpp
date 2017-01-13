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
int Calculator::iteration = 0;
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
                                    QPen(Qt::black), QBrush(Qt::SolidPattern));
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

void Calculator::calc_particle_dir(Particle &a, vector<Cell *> &for_sum_calculating) {

    if (debug_mode && &a == debug_draw_part) {
        double rad = 5;
        draw_debug(for_sum_calculating);
        scene_debug->addEllipse(debug_draw_part->X() - rad, debug_draw_part->Y() - rad, rad * 2.0, rad * 2.0,
                                QPen(Qt::red), QBrush(Qt::SolidPattern));
    }
    ////пройдем по реальным
    const vector<Particle *> *real;
    const vector<Particle> *shadow;
    const vector<Particle*> * bound;
    for (int i = 0; i < for_sum_calculating.size(); ++i) {
        real = for_sum_calculating[i] -> get_real();
        for (int j = 0; j < real -> size(); ++j) {
            const Particle &b = *real->operator[](j);
            a.calculate_derivatives_dis(b);
        }
        //пройдем по волшебным
        if(for_sum_calculating[i] -> is_boundary()) {
            shadow = for_sum_calculating[i]->get_shadow();
            for (int j = 0; j < shadow->size(); ++j) {
                const Particle &b = shadow->operator[](j);
                a.calculate_derivatives_dis(b);
            }
            //
            bound = for_sum_calculating[i]->get_bound();
            for (int j = 0; j < bound->size(); ++j) {
                const Particle &b = *bound->operator[](j);
               // a.calculate_derivatives_dis(b);
            }
        }

    }
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
                //жобавили саму себя
                for_sum_calculating.push_back(&parsing->part_groups[i][j]);
                //для частиц каждой каждой группы высчитываются производные
                calc_target_derivatives(parsing->part_groups[i][j], for_sum_calculating);
            }
        }
    }
}


//посчитаем все производные для частиц в target
void Calculator::calc_target_derivatives(Cell &target, vector<Cell *> &for_sum_calculating) {
//посчитаем
    unsigned long size = target.real_group.size();
    double res = 0;
    for (int i = 0; i < size; ++i) {
        Particle &target_p = *target.real_group[i];
        calc_particle_dir(target_p, for_sum_calculating);
    }
}

void Calculator::recalculate_parameters() {
    unsigned long data_size = parsing->data_ptr->size();
    double cur_maxV = 0;
    double cur_maxC = 0;
    double cur_maxP = 0;
    double cur_minP = 0;
    double new_dt = 100;
    for (unsigned long i = 0; i < data_size; ++i) {
        Particle *particle = &parsing->data_ptr->operator[](i);
        double cur_dt = particle->disser_integr(dt);
        if(cur_dt < new_dt) new_dt = cur_dt;
        cur_maxV = particle->V();
        cur_maxC = particle->C();
        cur_maxP = particle->P();
        cur_minP = particle->P();
        if (cur_maxV > max_V) max_V = cur_maxV;
        if (cur_maxC > max_C) max_C = cur_maxC;
        if (cur_maxP > max_P) max_P = cur_maxP;
        if (cur_minP < min_P) min_P = cur_minP;
    }
    //dt = new_dt;
}


void Calculator::calculate() {

    max_V = 0;
    max_C = 0;
    //всякие информативные штуки
    max_P = 0;
    min_P = 2e+12;

    parsing->clear_symetric_groups();
    parsing->create_symetric_groups();
    calculate_derivatives();
    recalculate_parameters();//пересчитали конечные параметры
    parsing->replace();
    //dt = h / (2 * max_C);
    current_time += dt;
    iteration++;
    //
}

//
Calculator::Calculator(SpaceParsing *target, double smooth_length) {
    //

    //
    parsing = target;
    debug_draw_part = &target->data_ptr->operator[](500);
    V_theoretical = sqrt(4 * 1.4 * parsing->data_ptr->operator[](0).E() / (1.4 - 1));
    max_V = 0;
    max_C = 0;
    dt = 1e-3;
    //всякие информативные штуки
    max_P = 0;
    min_P = 2e+14;
    //
    debug_mode = false;
    //h = parsing->x_size / 15;
    //h = fabs((parsing->data_ptr->operator[](0).X() - parsing->data_ptr->operator[](1).X())) * 5;
    h = parsing->h;
}