//
// Created by nikita on 30.05.16.
//


#include "Calculator.h"

namespace calculations {

    double current_time = 0;

    double h = 2;

//для подсчета силы от гарничных частиц
    double r0 = 10;

    double D = 100;

    double n1 = 4;

    double n2 = 12;

    double deltaT = 0.001;//требует глубокого переосмысления

    double r_ij(Particle &a, Particle &b) {
        return sqrt(pow(a.X() - b.X(), 2) + pow(a.Y() - b.Y(), 2));
    }

//
    double w_test(double r) {
        double result = 0;
        double x = r / h;
        if (x >= 0 && x <= 1) {
            result = (1 / 3.14 * pow(h, 3)) * (1 - (3.0 / 2.0) * pow(x, 2) + (3.0 / 4.0) * pow(x, 3));
        }
        else if (x >= 1 && x <= 2) {
            result = (1 / 3.14 * pow(h, 3)) * (1.0 / 4.0) * pow(2 - x, 3);
        }
        return result;
    }

//
    double grad_w_test(double r) {
        double result = 0;
        double x = r / h;
        if (x >= 0 && x <= 1) {
            result = (1 / 3.14 * pow(h, 4)) * ((9.0 / 4.0) * pow(x, 2) - 3.0 * x);
        }
        else if (x >= 1 && x <= 2) {
            result = (1 / 3.14 * pow(h, 4)) * (3.0 / 4.0) * -pow(2 - x, 2);
        }
        return result;
    }

//
    double two_part_p(Particle &a, Particle &b) {
        double M = b.M();
        double vx = a.Vx() - b.Vx();
        double vy = a.Vy() - b.Vy();
        //скалярное произведение Vil и Wij
        double d_Wx = grad_w_test(b.X() - a.X());
        double d_Wy = grad_w_test(b.Y() - a.Y());
        double VijWij = (vx * d_Wx) + (vy * d_Wy);
        double res = M * VijWij;
        return res;
    }

//
    double two_part_v(Particle &a, Particle &b, bool direct) {
        double M = b.M();
        double r = abs(a.X() - b.X());
        if (direct == 1) r = abs(a.Y() - b.Y());
        double Wij = grad_w_test(r);
        //
        double brackets = (a.P() / pow(a.p(), 2) +
                           b.P() / pow(b.p(), 2));
        //
        double res = M * brackets * Wij;
        return res;
    }

//
    double two_part_e(Particle &a, Particle &b) {
        double vx = a.Vx() - b.Vx();
        double vy = a.Vy() - b.Vy();
        //скалярное произведение Vil и Wij
        double d_Wx = grad_w_test(b.X() - a.X());
        double d_Wy = grad_w_test(b.Y() - a.Y());
        double VijWij = (vx * d_Wx) + (vy * d_Wy);
        //
        double M = b.M();
        double brackets = (a.P() / pow(a.p(), 2) +
                           b.P() / pow(b.p(), 2));
        //
        double res = 0.5 * VijWij * M * brackets;
        return res;

    }

//высчитывает часть производной скорости от действия граничной частицы
    double two_part_bforse(Particle &a, Particle& b, bool direct) {
        double r = calculations::r_ij(a, b);
        double ratio = r0 / r;
        if (ratio > 1) return 0;
        else {
            double x = abs(a.X() - b.X());
            if (direct == 1) x = abs(a.Y() - b.Y());
            double res = D * (pow(ratio, n1) - pow(ratio, n2)) * x / r;
            res = res / a.M();
            return res;
        }
    }

    //
//Вот эти 3 функции считают соответствующую  производную для частицы а векторы, переддаваемые
//в функции содержат все частицы, которые потенциально могут прореагировать с а

/*
 * просто без задней мысли идем по всем реальным частицам и находим часть
 * от взаимодействия целевой Particle &a и all_real[i][j]
 *
 */
    double calc_p_d(Particle &a, PartPointers &all_real, PartPointers_add &all_add) {
        double res = 0;
        for (int i = 0; i < all_real.size(); ++i) {

            for (int j = 0; j < all_real[i]->size(); ++j) {
                res += calculations::two_part_p(a, *all_real[i]->operator[](j));
            }
        }
        return res;
    }

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
                    PartPointers_add &all_add, bool direct) {
        double res = 0;
        //здесь м посчитали часть производной скорости от других реальныз частиц
        for (int i = 0; i < all_real.size(); ++i) {

            for (int j = 0; j < all_real[i]->size(); ++j) {
                res += calculations::two_part_v(a, *all_real[i]->operator[](j), direct);
            }
        }
        //теперь тут надо посчитать от симетричных частиц за границей
        for (int i = 0; i < all_add.size(); ++i) {

            for (int j = 0; j < all_add[i]->size(); ++j) {
                res += calculations::two_part_v(a, all_add[i]->operator[](j), direct);
            }
        }
        res = -res;
        //а вот тут надо посчитать производную от действия граничных частиц
        for (int i = 0; i < all_bound.size(); ++i) {

            for (int j = 0; j < all_bound[i]->size(); ++j) {
                res += calculations::two_part_bforse(a, *all_bound[i]->operator[](j), direct);
            }
        }
        return res;
    }

//
    double calc_E_d(Particle &a,
                    PartPointers &all_real,
                    PartPointers_add &all_add) {
        //здесь м посчитали часть производной энергии от других реальныз частиц
        double res = 0;
        for (int i = 0; i < all_real.size(); ++i) {
            for (int j = 0; j < all_real[i]->size(); ++j) {
                res += calculations::two_part_e(a, *all_real[i]->operator[](j));
            }
        }
        return res;
    }
}


//
//
void *Calculator::get_part_around(int row, int column, int i, int type) {
    int r_row = row;
    int r_column = column;
    switch (i) {
        case 0:
            r_row--;
            if (r_row >= 0) break;
            else return NULL;
        case 1:
            r_column--;
            r_row--;
            if (r_column >= 0 && r_row >= 0) break;
            else return NULL;
        case 2:
            r_column--;
            if (r_column >= 0) break;
            else return NULL;
        case 3:
            r_column--;
            r_row++;
            if (r_column >= 0 && r_row < parsing->cells_per_y) break;
            else return NULL;
        case 4:
            r_row++;
            if (r_row < parsing->cells_per_y) break;
            else return NULL;
        case 5:
            r_column++;
            r_row++;
            if (r_row < parsing->cells_per_y && r_column < parsing->cells_per_x) break;
            else return NULL;
        case 6:
            r_column++;
            if (r_column < parsing->cells_per_x) break;
            else return NULL;
        case 7:
            r_column++;
            r_row--;
            if (r_row >= 0 && r_column < parsing->cells_per_x) break;
            else return NULL;
        default:
            return NULL;
    }

    if (type == 0) {
        return &parsing->part_groups[r_row][r_column].real_group;
    }
    else if (type == 1) return &parsing->part_groups[r_row][r_column].boundary_group;

    else return &parsing->part_groups[r_row][r_column].symetric_group;
}

//функция должна собрать все реальные частицы из  Cellов вокруг i,j Cell
PartPointers Calculator::get_real_part(int &row, int &column) {
    vector<vector<Particle *> *> result;
    result.push_back(&(parsing->part_groups[row][column].real_group));
    for (int i = 0; i < 8; ++i) {
        vector<Particle *> *cur = (vector<Particle *> *) get_part_around(row, column, i, 0);
        if (cur != NULL) {
            result.push_back(cur);
        }
    }
    return result;
}

//
PartPointers Calculator::get_bound_part(int &row, int &column) {
    vector<vector<Particle *> *> result;
    result.push_back(&(parsing->part_groups[row][column].boundary_group));
    for (int i = 0; i < 8; ++i) {
        vector<Particle *> *cur = (vector<Particle *> *) get_part_around(row, column, i, 1);
        if (cur != NULL) {
            result.push_back(cur);
        }
    }
    return result;
}

//
PartPointers_add Calculator::get_sym_part(int &row, int &column) {
    vector<vector<Particle> *> result;
    result.push_back(&(parsing->part_groups[row][column].symetric_group));
    for (int i = 0; i < 8; ++i) {
        vector<Particle> *cur = (vector<Particle> *) get_part_around(row, column, i, 2);
        if (cur != NULL) {
            result.push_back(cur);
        }
    }
    return result;
}


//вот этот метод дергается в методе calculate, который во Flow
void Calculator::calculate_derivatives() {
    /*для каждой i,j области пространства создаются вектора всех реальных, граничных и теневых
     * частиц с которыми потенциально может провзаимодействовать любая частица в этой i,j области
     */
    for (int i = 0; i < parsing->cells_per_y; ++i) {
        for (int j = 0; j < parsing->cells_per_x; ++j) {

            PartPointers all_real = get_real_part(i, j);

            PartPointers all_bound = get_bound_part(i, j);

            PartPointers_add all_add = get_sym_part(i, j);

            //для частиц каждой каждой группы высчитываются производные и запихиваются в соответствующие
            // вектрора этих производных
            calc_p_derivatives(parsing->part_groups[i][j], all_real, all_add);

            calc_vx_derivatives(parsing->part_groups[i][j], all_real, all_bound, all_add);

            calc_vy_derivatives(parsing->part_groups[i][j], all_real, all_bound, all_add);

            calc_e_derivatives(parsing->part_groups[i][j], all_real, all_add);
        }
    }
}

//чтобы посчитать производную плотности нам нужно иметь список всех реальных и ВОЗМОЖНО теневых частиц
void Calculator::calc_p_derivatives(Cell &target,
                                    PartPointers &all_real,
                                    PartPointers_add &all_add) {

    for (int i = 0; i < target.real_group.size(); ++i) {
        double res = calculations::calc_p_d(*target.real_group[i], all_real, all_add);
        p_derivatives.push_back(res);
    }
}

//
void Calculator::calc_vx_derivatives(Cell &target,
                                     PartPointers &all_real,
                                     PartPointers &all_bound,
                                     PartPointers_add &all_add) {
    for (int i = 0; i < target.real_group.size(); ++i) {
        double res = calculations::calc_V_d(*target.real_group[i], all_real, all_bound, all_add, 0);
        vx_derivatives.push_back(res);
    }
}

//
void Calculator::calc_vy_derivatives(Cell &target,
                                     PartPointers &all_real,
                                     PartPointers &all_bound,
                                     PartPointers_add &all_add) {
    for (int i = 0; i < target.real_group.size(); ++i) {
        double res = calculations::calc_V_d(*target.real_group[i], all_real, all_bound, all_add, 1);
        vy_derivatives.push_back(res);
    }
}

//
void Calculator::calc_e_derivatives(Cell &target,
                                    PartPointers &all_real,
                                    PartPointers_add &all_add) {
    for (int i = 0; i < target.real_group.size(); ++i) {
        double res = calculations::calc_E_d(*target.real_group[i], all_real, all_add);
        e_derivatives.push_back(res);
    }

}


//
/*здесь мы идем по всем Cell так же как мы шли при рассчете производных
 * и для частиц каждого Cell рассчитываем значения параметров на t+dt временном слое
 * этим занимается функция Calculator::calc_t_values(Cell &target, int &row, int &colum)
 * после того как значения параметров на новом временном слое рассчитаны, а так же
 * заполнен вектор for_replacement который хранит пары, где второй элемент указатель на Cell
 * а 1й номер-список, хранящий пары из пары индексов Cell куда нужно переместь и индекса частицы
 * в векторе rel_part у Cell*
 * Функция replace(from_second_replace &replase_inf) в цикле производит необходимые перемещения
 * */

void Calculator::calculate_final() {
    for (int i = 0; i < parsing->cells_per_y; ++i) {
        for (int j = 0; j < parsing->cells_per_x; ++j) {
            calc_t_values(parsing->part_groups[i][j], i, j);
            //заодно в этом же цикле очистим вектора виртуальных частиц
            parsing->part_groups[i][j].symetric_group.clear();
        }
    }
    //теперь переместим частицы в другие cell
    for (int i = 0; i < for_replacement.size(); ++i) {
        replace(for_replacement[i]);
    }
}
//
//
/*
 * вот это краегольно важная функция для рассчета значения на новом временном слое
 * в нее передается Cell и индексы этого Cell  в parsing->part_groups
 * индесы нужны для рассчета позиции производной для этой частицы в векторе производных
 *
 *
*/
void Calculator::calc_t_values(Cell &target, int &row, int &colum) {
    //вектор содержащий частицы вылетившие из target
    //он содержит пары: новые i,j и индекс i вылетившей частицы
    vector<second_to_first> particles_to_rebase;
    for (int i = 0; i < target.real_group.size(); ++i) {
        //индекс в одномерном векторе производных для частиц
        int index = row * parsing->cells_per_y * parsing->cells_per_x + colum * target.real_group.size() + i;

        //вернули новые параметры
        Particle new_param = ronge_cutt(*target.real_group[i], index);
        //определили, нужно ли перемещать частицу
        std::pair<int, int> new_indexes = rebaze(new_param, row, colum);
        if (new_indexes.first != -1) {
            //если нужно, то запомнили индексы куда переместить и индекс частицы в векторе
            //target.real_group которую нужно переместить
            pair<pair<int, int>, int> new_pair(new_indexes, i);
            particles_to_rebase.push_back(new_pair);
        }
        //теперь после того как мы нашли новые параметры и выяснили, нужно ли что то перемещать
        //можно старой частеце установить новые параметры
        target.real_group[i]->set_from(new_param);
    }
    //теперь для этой частицы сформирован вектор частиц, которые нужно переместить
    for_replacement.push_back(pair<vector<second_to_first>, Cell *>(particles_to_rebase, &target));
}

//пока что это просто эйлер этот метод находит параметры на +dt слое и возвращает их
Particle Calculator::ronge_cutt(Particle &a, int &index) {
    double dt = calculations::deltaT;
    double mass = a.M();
    double new_p = a.p() + dt * p_derivatives[index];
    double new_e = a.E() + dt * e_derivatives[index];
    double new_vx = a.Vx() + dt * vx_derivatives[index];
    double new_vy = a.Vy() + dt * vy_derivatives[index];
    //    double e = P / (0.4 * p);
    /*
      Particle(bool status, double p = 0,
             double P = 0, double E = 0, double Vx = 0, double Vy = 0, double M = 0);
      */
    double new_P = new_e * 0.4 * new_p;
    Particle result(0, new_p, new_P, new_e, new_vx, new_vy, mass);
    double new_x = a.X() + dt * vx_derivatives[index];
    double new_y = a.Y() + dt * vy_derivatives[index];
    Point new_pos(new_x, new_y);
    result.set_pos(new_pos);
    return result;
}

/*этот метод определяет, нужно ли перемещать частицу new_part из Cell
 * и индексами int &row, int &colum
 * и если все таки нужно то куда
 * куда возвращается в виде пары индексов
*/
std::pair<int, int> Calculator::rebaze(Particle &new_part, int &row, int &colum) {
    if (parsing->part_groups[row][row].is_inside(new_part) == true) {
        return pair<int, int>(-1, -1);
    }
    else {
        //если частица больше не внутри своей старой позиции то мы найдем индексы нового cell
        std::pair<int, int> new_indexes = parsing->find_around(row, colum, new_part);
        return new_indexes;
    }
}
//

/*
 * переместить все частицы которые нужно перемещать
 * информации о них хранится в ячейках вектора for_replacement
 */
void Calculator::replace(from_second_replace &replase_inf) {
    //typedef pair<pair<int,int>, int> second_to_first;
    //typedef pair<vector<second_to_first>, Cell*> from_second_replace;
    for (int count = 0; count < replase_inf.first.size(); ++count) {
        //индексы Cell куда мы переместим частицы из Cell*
        int i = replase_inf.first[count].first.first;
        int j = replase_inf.first[count].first.second;
        int replace_id = replase_inf.first[count].second;
        //добавили эту частицу
        parsing->part_groups[i][j].real_group.push_back(replase_inf.second->real_group[replace_id]);
        //теперь удалим этот элемент из вектора
        swap(replase_inf.second->real_group[replace_id], replase_inf.second->real_group.back());
        replase_inf.second->real_group.pop_back();
    }

}

//основной метод, вызываемый снаружи
/*
 * в нем сначала дергается calculate_derivatives();-метод, высчитывающий производные
 * затем метод calculate_final();-метод, который высчитывает новые значения для частиц на новом шаге по времени
 * calculations::curent_time+=calculations::deltaT; обновляет счетчик времени
*/
void Calculator::calculate() {

    calculate_derivatives();

    calculate_final();

    calculations::current_time += calculations::deltaT;
    // очистим все вспомогательные вектора
    p_derivatives.clear();
    vx_derivatives.clear();
    vy_derivatives.clear();
    e_derivatives.clear();
    for_replacement.clear();
    //заново нужно создать теневые частицы
    parsing->create_symetric_groups();

}
