//
// Created by nikita on 30.05.16.
//

#ifndef SPHSM6_CALCULATOR_H
#define SPHSM6_CALCULATOR_H

#include "SpaceParsing.h"

typedef pair<pair<int, int>, int> second_to_first;
//из этого cell нужно переместить такие то частицы туда то туда то
typedef pair<vector<second_to_first>, Cell *> from_second_replace;

namespace calculations {

    extern double current_time;

    extern double h;

    //для подсчета силы от гарничных частиц
    extern double r0;

    extern double D;

    extern double n1, n2;

    extern double deltaT;//требует глубокого переосмысления

    double r_ij(Particle &a, Particle &b);

    //
    double w_test(double r);

    //
    double grad_w_test(double r);

    //
    double two_part_p(Particle &a, Particle &b) ;

    //
    double two_part_v(Particle &a, Particle &b, bool direct) ;

    //
    double two_part_e(Particle &a, Particle &b) ;

    //высчитывает часть производной скорости от действия граничной частицы
    double two_part_bforse(Particle &a, Particle& b, bool direct) ;

    //
    //Вот эти 3 функции считают соответствующую  производную для частицы а векторы, переддаваемые
    //в функции содержат все частицы, которые потенциально могут прореагировать с а

    /*
     * просто без задней мысли идем по всем реальным частицам и находим часть
     * от взаимодействия целевой Particle &a и all_real[i][j]
     *
     */
    double calc_p_d(Particle &a, PartPointers &all_real, PartPointers_add &all_add) ;

    //
    /*
    * просто без задней мысли идем по всем реальным граничным и теневым частицам и находим часть
    * от взаимодействия целевой Particle &a и all_real[i][j]
    * от взаимодействия целевой Particle &a и all_add[i][j]
    * от взаимодействия целевой Particle &a и all_bound[i][j]
    */
    double calc_V_d(Particle &a,
                    PartPointers &all_real,
                    PartPointers &all_bound,
                    PartPointers_add &all_add, bool direct);

    //
    double calc_E_d(Particle &a,
                    PartPointers &all_real,
                    PartPointers_add &all_add);
}

class Calculator {
    vector<from_second_replace> for_replacement;

    SpaceParsing *parsing;

    vector<double> p_derivatives;

    vector<double> vx_derivatives;

    vector<double> vy_derivatives;

    vector<double> e_derivatives;

    //
    PartPointers_add get_sym_part(int &row, int &column);

    //
    PartPointers get_bound_part(int &row, int &column);

    //
    PartPointers get_real_part(int &row, int &column);

    //
    void *get_part_around(int row, int column, int i, int type);

    //
    void calc_p_derivatives(Cell &target,
                            PartPointers &all_real,
                            PartPointers_add &all_add);

    //
    void calc_vx_derivatives(Cell &target,
                             PartPointers &all_real,
                             PartPointers &all_bound,
                             PartPointers_add &all_add);

    //
    void calc_vy_derivatives(Cell &target,
                             PartPointers &all_real,
                             PartPointers &all_bound,
                             PartPointers_add &all_add);

    //
    void calc_e_derivatives(Cell &target,
                            PartPointers &all_real,
                            PartPointers_add &all_add);

    //
    void calc_t_values(Cell &target, int &row, int &colum);

    //
    Particle ronge_cutt(Particle &a, int &index);

    //
    std::pair<int, int> rebaze(Particle &new_part, int &row, int &colum);

    void calculate_derivatives();

    void calculate_final();

    void replace(from_second_replace &replase_inf);

public:

    Calculator(SpaceParsing *target) {
        parsing = target;
    }

    void calculate();

};


#endif //SPHSM6_CALCULATOR_H
