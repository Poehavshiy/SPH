//
// Created by nikita on 30.05.16.
//

#ifndef SPHSM6_CALCULATOR_H
#define SPHSM6_CALCULATOR_H

#include "SpaceParsing.h"

struct From_cell_to_cells {
    Cell *from;
    vector<pair<Cell *, int>> to_cell_pid;

    From_cell_to_cells(Cell *f, vector<pair<Cell *, int>> &targets) {
        from = f;
        to_cell_pid = targets;
    }

    int size() {
        return to_cell_pid.size();
    }

    void replace() {
        vector<int> del;
        for (int i = 0; i < to_cell_pid.size(); ++i) {
            del.push_back(to_cell_pid[i].second);
            Particle *remove = from->get_real()->operator[](del.back());
            to_cell_pid[i].first->add_part(remove);
        }
        std::sort(del.begin(), del.end(), greater<int>());
        for (int i = 0; i < del.size(); ++i) {
            from->remove_part(del[i]);
        }
    }
};

namespace calculations {

    extern double current_time;

    extern double h;

    extern double R;

    extern double M;

    extern double k;

    extern double alpha;

    extern double beta;

    //для подсчета силы от гарничных частиц
    extern double r0;

    extern double D;

    extern double n1, n2;

    extern double deltaT;//требует глубокого переосмысления

    extern double viscous;

    double r_ij(Particle &a, Particle &b);

    //
    double w_test(double r);

    //
    double grad_w_test(double x, double y, bool direction);

    //
    double two_part_p(Particle &a, Particle &b);

    //
    double two_part_E(Particle &a, Particle &b, bool direct);

    //
    double two_part_art_visc(Particle &a, Particle &b);

    //
    double two_part_art_heat(Particle &a, Particle &b, bool direct);

    //
    double two_part_v(Particle &a, Particle &b, bool direct);

    //
    double two_part_e(Particle &a, Particle &b);

    //высчитывает часть производной скорости от действия граничной частицы
    double two_part_bforse(Particle &a, Particle &b, bool direct);

    //
    //Вот эти 3 функции считают соответствующую  производную для частицы а векторы, переддаваемые
    //в функции содержат все частицы, которые потенциально могут прореагировать с а

    /*
     * просто без задней мысли идем по всем реальным частицам и находим часть
     * от взаимодействия целевой Particle &a и all_real[i][j]
     *
     */
    double calc_p_d(Particle &a, PartPointers &all_real, PartPointers_add &all_add);

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
protected:
    int index = 0;

    double largest_V = 0;

    double dt = 0.1;


    vector<From_cell_to_cells> for_replacement;

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
    void calc_t_values(Cell &target, const int &row, const int &colum);

    //
    Particle ronge_cutt(Particle &a, int &index);

    //
    std::pair<int, int> rebaze(Particle &new_part, const int &row, const int &colum);

    virtual void calculate_derivatives();

    void calculate_final();

    void recalculate_consts(); //пересчитывает D, deltaT

    bool get_maxv_mindt(Cell &target, double &cur_max_v, double &cur_min_dt);

    void replace();

    void check_empty();

    int calc_non_empty(vector<std::pair<int, int>> &indexes_of_nonempty);

public:

    Calculator(SpaceParsing *target) {
        parsing = target;
    }

    void calculate();

};


#endif //SPHSM6_CALCULATOR_H
