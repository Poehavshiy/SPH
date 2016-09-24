//
// Created by nikita on 30.05.16.
//

//тут старое вычисление скоростей
/*
    double calc_V_d(Particle &a, Cell &home_cell, vector<Cell *> &for_sum_calculating, bool direct) {
        if (for_debugin == &a) {
            //draw_debug(all_real, all_add, all_bound);
            double rad = 5;
            scene_debug->addEllipse(a.X() - rad, a.Y() - rad, rad * 2.0, rad * 2.0,
                                    QPen(Qt::red), QBrush(Qt::SolidPattern));
        }
        double res = 0;
        //здесь м посчитали часть производной скорости от других реальныз частиц
        for_sum_calculating.push_back(&home_cell);
        for (int i = 0; i < for_sum_calculating.size(); ++i) {
            const vector<Particle *> *real = for_sum_calculating[i]->get_real();
            for (int j = 0; j < real->size(); ++j) {
                res += calculations::two_part_v(a, *real->operator[](j), direct);
            }
        }
        //теперь от виртуальных частиц
        //теперь тут надо посчитать от симетричных частиц за границей
        for (int i = 0; i < for_sum_calculating.size(); ++i) {
            const vector<Particle> *shadow = for_sum_calculating[i]->get_shadow();
            for (int j = 0; j < shadow->size(); ++j) {
                res += calculations::two_part_v(a, shadow->operator[](j), direct);
            }
        }
        //а вот тут надо посчитать производную от действия граничных частиц
        for (int i = 0; i < for_sum_calculating.size(); ++i) {
            const vector<Particle *> *bound = for_sum_calculating[i]->get_bound();
            for (int j = 0; j < bound->size(); ++j) {
                res += calculations::two_part_bforse(a, *bound->operator[](j), direct);
            }
        }
        return res;
    }
    */


#ifndef SPHSM6_CALCULATOR_H
#define SPHSM6_CALCULATOR_H

#include "SpaceParsing.h"

const bool EULER = 1;
#define FUNCTION_NUMBER 4;
//
class Calculator {
protected:
    Particle* debug_draw_part;
    bool debug_mode ;
    double max_V = 0;
    double max_C = 0;
    double dt = 0.001;
    //всякие информативные штуки
    double V_theoretical;
    double max_P = 0;
    double min_P = 0;
    function<double(const Particle &, const Particle &, bool)> derivative_functions[4];

    SpaceParsing *parsing;
    //с нэймспэйса
    static double h ;
    static double R  ;
    static double k  ;
    static double M  ;
    static double viscous ;
    static double alpha ;
    static double beta ;
//для подсчета силы от гарничных частиц
    static double r0;
    static double D;//равен квадрату наибольшей скорости
    static double n1;
    static double n2;
    static  int iteration;
    //бывшие функции нэймспэйса
    static void delta_coor(const Particle &a, const Particle &b, double &deltaX, double &deltaY);

    static void delta_velocity(const Particle &a, const Particle &b, double &delta_Vx, double &delta_Vy);

    //
    static double grad_w_test(double x, double y, bool direction, double h_inc);

    //
    static double two_part_p(const Particle &a, const Particle &b, bool direct );

    static double two_part_v(const Particle &a, const Particle &b, bool direct );

    static double two_part_energy(const Particle &a, const Particle &b, bool direct );

    //
    static double two_part_E(const Particle &a, const Particle &b, bool direct);

    //
    static double two_part_art_visc(const Particle &a, const Particle &b);

    //
    static double two_part_art_heat(const Particle &a, const Particle &b, bool direct);

    //высчитывает часть производной скорости от действия граничной частицы
    inline static double two_part_bforse(const Particle &a, const Particle &b, bool direct);

    inline double calc_particle_dir(function<double(const Particle &, const Particle &, bool)>& derivative_specific,
            Particle &a, vector<Cell *> &for_sum_calculating, bool direction, DERIVATIVES &what_dir);

    //старые функции, что тут всегда были
    inline void calc_target_derivatives(function<double(const Particle &, const Particle &, bool)> derivative_specific,
                                 Cell &target, vector<Cell *> &for_sum_calculating,
                                 DERIVATIVES &what_dir, bool direction,
                                 int &row, int &column);

    inline void calculate_derivatives();

//тупо пересчет новых значений физ величин
    inline void recalculate_parameters();

    double boundary_forse(Particle& a, vector<Cell *> &for_sum_calculating, bool direction);

    bool to_boundary(Particle& a, Cell &target);

    bool to_boundary(Particle& a, Point& A, Point& B);

public:
    static double current_time;

    Calculator(SpaceParsing *target);
    //основной метод, вызываемый снаружи
/*
 * в нем сначала дергается calculate_derivatives();-метод, высчитывающий производные
 * затем метод calculate_final();-метод, который высчитывает новые значения для частиц на новом шаге по времени
 * calculations::curent_time+=calculations::deltaT; обновляет счетчик времени
*/
    void calculate();

    void draw_debug(vector<Cell *> &for_sum_calculating);

    int get_iteration(){
        return iteration;
    }

    double get_time(){
        return current_time;
    }

    double get_V_theory(){
        return V_theoretical;
    }

    double get_maxV(){
        return max_V;
    }

    double get_maxP(){
        return max_P;
    }

    double get_minP(){
        return min_P;
    }

};




#endif //SPHSM6_CALCULATOR_H
