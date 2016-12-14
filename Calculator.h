//
// Created by nikita on 30.05.16.
//

#ifndef SPHSM6_CALCULATOR_H
#define SPHSM6_CALCULATOR_H

#include "SpaceParsing.h"

const bool EULER = 1;
#define FUNCTION_NUMBER 4;

//
class Calculator {
protected:
    Particle *debug_draw_part;
    bool debug_mode;
    double max_V = 0;
    double max_C = 0;
    //всякие информативные штуки
    double V_theoretical;
    double max_P = 0;
    double min_P = 0;
    SpaceParsing *parsing;

    static int iteration;
    static double h;

    inline void calc_particle_dir(Particle &a, vector<Cell *> &for_sum_calculating);

    //старые функции, что тут всегда были
    inline void calc_target_derivatives(Cell &target, vector<Cell *> &for_sum_calculating);

    inline void calculate_derivatives();

//тупо пересчет новых значений физ величин
    inline void recalculate_parameters();


public:
    static double current_time;
    double dt;

    Calculator(SpaceParsing *target, double smooth_length);
    //основной метод, вызываемый снаружи
/*
 * в нем сначала дергается calculate_derivatives();-метод, высчитывающий производные
 * затем метод calculate_final();-метод, который высчитывает новые значения для частиц на новом шаге по времени
 * calculations::curent_time+=calculations::deltaT; обновляет счетчик времени
*/
    void calculate();

    void draw_debug(vector<Cell *> &for_sum_calculating);

    int get_iteration() {
        return iteration;
    }

    double get_time() {
        return current_time;
    }

    double get_V_theory() {
        return V_theoretical;
    }

    double get_maxV() {
        return max_V;
    }

    double get_maxP() {
        return max_P;
    }

    double get_minP() {
        return min_P;
    }

};


#endif //SPHSM6_CALCULATOR_H
