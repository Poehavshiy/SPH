//
// Created by nikita on 30.05.16.
//


#include "Calculator.h"
#include "Flow.h"

QGraphicsScene *scene_debug;
//typedef vector<vector<Particle*>*> PartPointers;
//typedef vector<vector<Particle>*> PartPointers_add;

void draw_debug(PartPointers &all_real, PartPointers_add &all_add, PartPointers &all_bound) {
    double rad = 2;
    for (int i = 0; i < all_real.size(); ++i) {
        for (int j = 0; j < all_real[i]->size(); ++j) {
            Particle *curent = all_real[i]->operator[](j);
            scene_debug->addEllipse(curent->X() - rad, curent->Y() - rad, rad * 2.0, rad * 2.0,
                                    QPen(Qt::green), QBrush(Qt::SolidPattern));
        }
    }
    //
    for (int i = 0; i < all_add.size(); ++i) {
        for (int j = 0; j < all_add[i]->size(); ++j) {
            Particle *curent = &all_add[i]->operator[](j);
            scene_debug->addEllipse(curent->X() - rad, curent->Y() - rad, rad * 2.0, rad * 2.0,
                                    QPen(Qt::blue), QBrush(Qt::SolidPattern));
        }
    }
    //
    for (int i = 0; i < all_bound.size(); ++i) {
        for (int j = 0; j < all_bound[i]->size(); ++j) {
            Particle *curent = all_bound[i]->operator[](j);
            scene_debug->addEllipse(curent->X() - rad, curent->Y() - rad, rad * 2.0, rad * 2.0,
                                    QPen(Qt::black), QBrush(Qt::SolidPattern));
        }
    }

}


namespace calculations {

    double current_time = 0;

    double h = 25;

    double R = 8.31;

    double k = 1.2;

    double M = 0.029;

    double viscous = 0.00023;

    double alpha = 1.1;

    double beta = 1.1;

//для подсчета силы от гарничных частиц
    double r0 = 2;

    double D = 0;//0.01;//равен квадрату наибольшей скорости

    double n1 = 12;

    double n2 = 6;

    double deltaT = 0.001;//требует глубокого переосмысления

    double r_ij(Particle &a, Particle &b) {
        return sqrt(pow(a.X() - b.X(), 2) + pow(a.Y() - b.Y(), 2));
    }

//
    /*  double w_test(double r) {
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
  */
    double w_test(double r) {
        // assert(r>0);
        double result = 0;
        if (r < 2) {
            result = (15 / (7 * 3.14 * h * h)) *
                     ((2 / 3 - 9 / 8 * pow(r, 2)) + 19 / 24 * pow(r, 3) - 5 / 32 * pow(r, 4));
        }
        return result;
    }
//вооооооооооооооооот тут косяк
    //1я функция
    /*  double grad_w_test(double r) {
          double result = 0;
          double x = r / h;
          if (x >= 0 && x <= 1) {
              result = (1 / (3.14 * pow(h, 4))) * ((9.0 / 4.0) * pow(x, 2) - 3.0 * x);
          }
          else if (x >= 1 && x <= 2) {
              result = (1 / (3.14 * pow(h, 4))) * (3.0 / 4.0) * -pow(2 - x, 2);
          }
          return result;
      }
  */
    //2я функция
    /*double grad_w_test(double x, double y, bool direction) {
        double result = 0;
        double r = sqrt(x*x+y*y)/h;
        if(r<2 && r>-2) {
            if(direction == 0) {
                result = (15 / (7 * 3.14 * h * h)) * x * (-5 * pow(r, 2) + 19 * r - 72);
            }
            else
                result = (15 / (7 * 3.14 * h * h)) * y * (-5 * pow(r, 2) + 19 * r - 72);
        }
        int a=1;
        return result;
    }*/
    double grad_w_test(double x, double y, bool direction) {
        double result = 0;
        double a = 5 / (3.14 * h * h);
        double r = sqrt(x * x + y * y);
        if (r > h) return 0;
        //x=abs(x);
        //y=abs(y);
        if (direction == 0)
            result = 12 * a * x * (pow(h, 2) - 2 * h * r + pow(r, 2)) / pow(h, 4);
        else
            result = 12 * a * y * (pow(h, 2) - 2 * h * r + pow(r, 2)) / pow(h, 4);

        return result;
    }

    double grad_w_test(double x, double y, bool direction, double h_inc) {
        double result = 0;
        double a = 5 / (3.14 * h_inc * h_inc);
        double r = sqrt(x * x + y * y);
        if (r > h) return 0;
        //x=abs(x);
        //y=abs(y);
        if (direction == 0)
            result = 12 * a * x * (pow(h_inc, 2) - 2 * h_inc * r + pow(r, 2)) / pow(h_inc, 4);
        else
            result = 12 * a * y * (pow(h_inc, 2) - 2 * h_inc * r + pow(r, 2)) / pow(h_inc, 4);

        return result;
    }
//
    double two_part_p(Particle &a, Particle &b) {
        double M = b.M();
        double vx = a.Vx() - b.Vx();
        double vy = a.Vy() - b.Vy();
        //скалярное произведение Vil и Wij
        double deltaX = b.X() - a.X();
        double deltaY = b.Y() - a.Y();
        double d_Wx = grad_w_test(deltaX, deltaX, 0, a.h());
        double d_Wy = grad_w_test(deltaX, deltaY, 1, a.h());
        double VijWij = (vx * d_Wx) + (vy * d_Wy);
        double res = M * VijWij;
        return res;
    }

    double two_part_E(Particle &a, Particle &b, bool direct) {
        //direct 0 -> вычисляем Exy, если direct 1 -> вычисляем Eyx
        double mj_pj = b.M()/b.p();
        double Vj_y_Wi_x;
        double Vj_x_Wi_y;
        double deltaX = b.X() - a.X();
        double deltaY = b.Y() - a.Y();
        if(direct == 0) {
            Vj_y_Wi_x = b.Vy()*grad_w_test(deltaX,deltaY, 0, a.h());
            Vj_x_Wi_y = b.Vx()*grad_w_test(deltaX,deltaY, 1, a.h());
        }
        else {
            Vj_y_Wi_x = b.Vx()*grad_w_test(deltaX,deltaY, 1, a.h());
            Vj_x_Wi_y = b.Vy()*grad_w_test(deltaX,deltaY, 0, a.h());
        }
        double result = mj_pj*(Vj_y_Wi_x+Vj_x_Wi_y);
        return result;
    }

    //
    double two_part_art_visc(Particle &a, Particle &b) {
        double aR = sqrt(pow(a.X(), 2) + pow(a.Y(), 2));
        double bR = sqrt(pow(b.X(), 2) + pow(b.Y(), 2));
        double deltaR = aR-bR;
        //
        double av = sqrt(pow(a.Vx(), 2) + pow(a.Vy(), 2));
        double bv = sqrt(pow(b.Vx(), 2) + pow(b.Vy(), 2));
        double deltaV = aR-bR;


        double phi = a.h()* deltaR * deltaV/(pow(deltaR,2) + pow(0.1*a.h(),2));
        double result;
        if(deltaR*deltaV < 0) result = 0;
        else result = 2*(-alpha * 0.5*(a.C() + b.C())*phi + beta * pow(phi,2))/(a.p() + b.p());
        return result;
    }

    //Теплопроводность

    double two_part_art_heat(Particle &a, Particle &b, bool direct) {
       /* double deltaX = b.X() - a.X();
        double deltaY = b.Y() - a.Y();
        double divVa = b.M() *(b.Vx()*grad_w_test(deltaX, deltaY, 0)+b.Vy()*grad_w_test(deltaX, deltaY, 1))/b.p();
        double qi = alpha*h*a.p()*a.C()*abs(divVa)+beta*h*a.p()*pow(divVa,2);

        double divVb = a.M() *(b.Vx()*grad_w_test(deltaX, deltaY, 0)+b.Vy()*grad_w_test(deltaX, deltaY, 1))/b.p();
        double qj = alpha*h*b.p()*b.C()*abs(divV)+beta*h*b.p()*pow(divV,2);*/
        return 0;

    }

    //
    double two_part_v(Particle &a, Particle &b, bool direct) {
        double M = b.M();
        double deltaX = b.X() - a.X();//abs(b.X() - a.X());
        double deltaY = b.Y() - a.Y();//abs(b.Y() - a.Y());
        double Wij = grad_w_test(deltaX, deltaY, direct, a.h());

        double brackets1 = (a.P() / pow(a.p(), 2) +
                           b.P() / pow(b.p(), 2));
        //
        double res = M * brackets1 * Wij;//1е слагаемоеж
        //считаем 2е
        double brackets2 = viscous*two_part_E(a, b, direct)/pow(a.p(), 2) +
                //второе слагаемое  Ej a b поменяли
                viscous*two_part_E(b, a, direct)/pow(b.p(), 2);
        res += M * brackets2 * Wij;
        //добавим иск вязкость
        double H = two_part_art_visc(a, b);

        res+= H*Wij;

        return res;
    }

//
    double two_part_e(Particle &a, Particle &b) {
        double vx = a.Vx() - b.Vx();
        double vy = a.Vy() - b.Vy();
        //скалярное произведение Vil и Wij
        double deltaX = b.X() - a.X();//abs(b.X() - a.X());
        double deltaY = b.Y() - a.Y();//abs(b.Y() - a.Y());
        double VijWij = (vx * grad_w_test(deltaX, deltaY, 0,a.h())) + (vy * grad_w_test(deltaX, deltaY, 1,a.h()));
        //
        double M = b.M();
        double H = two_part_art_visc(a, b);
        double brackets1 = (a.P() / pow(a.p(), 2) +
                           b.P() / pow(b.p(), 2) + H);
        //
        double res = 0.5 * VijWij * M * brackets1;
        //вязкость
        double brackets2 = (viscous/(2*a.p()))*two_part_E(a,b,0) * two_part_E(a,b,1);
        res+= brackets2;

        return res;

    }

//высчитывает часть производной скорости от действия граничной частицы
    double two_part_bforse(Particle &a, Particle &b, bool direct) {
        double r = calculations::r_ij(a, b);
        double ratio = r0 / r;
        if (ratio > 1) return 0;
        else {
            double x = (a.X() - b.X());
            if (direct == 1) x = (a.Y() - b.Y());
            double res = D * (pow(ratio, n1) - pow(ratio, n2)) * x / r;
            res = res / a.M();
            return -res;
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
        //теперь от виртуальных частиц
        //теперь тут надо посчитать от симетричных частиц за границей
        for (int i = 0; i < all_add.size(); ++i) {

            for (int j = 0; j < all_add[i]->size(); ++j) {
                res += calculations::two_part_p(a, all_add[i]->operator[](j));
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
        if (for_debugin == &a) {
            draw_debug(all_real, all_add, all_bound);
            double rad = 5;
            scene_debug->addEllipse(a.X() - rad, a.Y() - rad, rad * 2.0, rad * 2.0,
                                    QPen(Qt::red), QBrush(Qt::SolidPattern));
        }
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
        //
        //теперь от виртуальных частиц
        //теперь тут надо посчитать от симетричных частиц за границей
        for (int i = 0; i < all_add.size(); ++i) {

            for (int j = 0; j < all_add[i]->size(); ++j) {
                res += calculations::two_part_e(a, all_add[i]->operator[](j));
            }
        }
        return res;
    }
}


//
//Calculator class
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
    } else if (type == 1) return &parsing->part_groups[r_row][r_column].boundary_group;

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
    //теперь переместим частицы в другие cells
    replace();
}
//
//
/*
 * вот это краеугольно важная функция для рассчета значения на новом временном слое
 * в нее передается Cell и индексы этого Cell  в parsing->part_groups
 * индесы нужны для рассчета позиции производной для этой частицы в векторе производных
*/
void Calculator::calc_t_values(Cell &target_cell, const int &row, const int &colum) {
    //Cell*-куда, ште - идекс в target
    vector<pair<Cell *, int>> to_cell_pid;
    for (int i = 0; i < target_cell.real_group.size(); ++i) {
        //вернули новые параметры
        Particle new_param = ronge_cutt(*target_cell.real_group[i], index);
        //определили, нужно ли перемещать частицу
        if (index == 3) {
            int r = 0;
        }
        std::pair<int, int> new_indexes = rebaze(new_param, row, colum);
        if (new_indexes.first != -1) {
            //если нужно, то запомнили Cell куда переместить и откуда
            //target.real_group которую нужно переместить
            pair<Cell *, int> replace_i_to(&parsing->part_groups[new_indexes.first][new_indexes.second], i);
            to_cell_pid.push_back(replace_i_to);
        }
        //теперь после того как мы нашли новые параметры и выяснили, нужно ли что то перемещать
        //можно старой частеце установить новые параметры
        target_cell.real_group[i]->set_from(new_param);
        ++index;
    }
    From_cell_to_cells new_replace_set(&target_cell, to_cell_pid);
    for_replacement.push_back(new_replace_set);
    //теперь для этого Cell сформирован вектор частиц, которые нужно переместить
}

//пока что это просто эйлер этот метод находит параметры на +dt слое и возвращает их
Particle Calculator::ronge_cutt(Particle &a, int &index) {
    double dt = calculations::deltaT;
    double mass = a.M();
    double new_p = a.p() + dt * p_derivatives[index];
    double new_e = a.E() + dt * e_derivatives[index];
    double new_vx = a.Vx() + dt * vx_derivatives[index];
    double new_vy = a.Vy() + dt * vy_derivatives[index];
    //
    //Пересчитаем D
    double new_V = pow(new_vx * new_vx + new_vy * new_vy, 0.5);
    if (new_V > largest_V) largest_V = new_V;
    //Пересчитаем delta_t
    double c = sqrt(1.4 * calculations::R * a.T() / calculations::M);
    double new_dt = calculations::h / c;
    if (new_dt < dt) dt = new_dt;
    //
    double new_P = new_e * 0.4 * new_p;
    Particle result(0, new_p, new_P, new_e, new_vx, new_vy, mass);
    double new_x = a.X() + dt * new_vx;
    double new_y = a.Y() + dt * new_vy;
    Point new_pos(new_x, new_y);
    result.set_pos(new_pos);
    return result;
}

/*этот метод определяет, нужно ли перемещать частицу new_part из Cell
 * и индексами int &row, int &colum
 * и если все таки нужно то куда
 * куда возвращается в виде пары индексов
*/
std::pair<int, int> Calculator::rebaze(Particle &new_part, const int &row, const int &colum) {
    if (parsing->part_groups[row][colum].is_inside(new_part) == true) {
        return pair<int, int>(-1, -1);
    } else {
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
void Calculator::replace() {
    for (int i = 0; i < for_replacement.size(); ++i) {
        for_replacement[i].replace();
    }
}

//
void Calculator::recalculate_consts() {
    double v_max = -100500;
    double min_dt = 100500;
    double cur_max_v, cur_min_dt;
    for (int i = 0; i < parsing->cells_per_y; ++i) {
        for (int j = 0; j < parsing->cells_per_x; ++j) {
            bool status = get_maxv_mindt(parsing->part_groups[i][j], cur_max_v, cur_min_dt);
            if (status != false) {
                if (cur_max_v > v_max) {
                    v_max = cur_max_v;
                }
                if (cur_min_dt < min_dt) {
                    min_dt = cur_min_dt;
                }
            }
        }
    }

    //непосредственно переприсвоим
    calculations::D = v_max * v_max;
    calculations::deltaT = min_dt;
}

bool Calculator::get_maxv_mindt(Cell &target, double &cur_max_vin, double &cur_min_dtin) {
    if (target.real_group.size() == 0) return false;
    else {
        double v_max = -100500;
        double min_dt = 100500;
        double cur_max_v, cur_min_dt;
        for (int i = 0; i < target.real_group.size(); ++i) {
            cur_max_v = sqrt(
                    pow(target.real_group.operator[](i)->Vx(), 2) + pow(target.real_group.operator[](i)->Vy(), 2));
            if (cur_max_v > v_max) {
                v_max = cur_max_v;
            }
            double C = sqrt(
                    (calculations::k * calculations::R / calculations::M) * target.real_group.operator[](i)->T());
            cur_min_dt = C / calculations::h;
            if (cur_min_dt < min_dt) {
                min_dt = cur_min_dt;
            }
        }
        cur_max_vin = cur_max_v;
        cur_min_dt = cur_min_dt;
        return true;
    }
}

void Calculator::check_empty() {
    for (int i = 0; i < parsing->cells_per_y; ++i) {

        for (int j = 0; j < parsing->cells_per_x; ++j) {
            parsing->part_groups[i][j].emptiness();
        }
    }
}

int Calculator::calc_non_empty(vector<std::pair<int, int>> &indexes_of_nonempty) {
    int counter = 0;
    for (int i = 0; i < parsing->cells_per_y; ++i) {
        for (int j = 0; j < parsing->cells_per_x; ++j) {
            if (parsing->part_groups[i][j].is_non_empty() == true) {
                ++counter;
                indexes_of_nonempty.push_back(pair<int, int>(i, j));
            }
        }
    }
    return counter;

}
//основной метод, вызываемый снаружи
/*
 * в нем сначала дергается calculate_derivatives();-метод, высчитывающий производные
 * затем метод calculate_final();-метод, который высчитывает новые значения для частиц на новом шаге по времени
 * calculations::curent_time+=calculations::deltaT; обновляет счетчик времени
*/
void Calculator::calculate() {

    check_empty();
    //для дебажтрования
    vector<std::pair<int, int>> indexes_of_nonempty;
    int how_many_nonempty = calc_non_empty(indexes_of_nonempty);
    //
    calculate_derivatives();
    calculate_final();
    recalculate_consts();
    int check = for_replacement.size();
    // очистим все вспомогательные вектора
    p_derivatives.clear();
    vx_derivatives.clear();
    vy_derivatives.clear();
    e_derivatives.clear();
    for_replacement.clear();
    //заново нужно создать теневые частицы
    parsing->create_symetric_groups();
    index = 0;
    calculations::D = largest_V * largest_V;
    calculations::deltaT = dt;
    calculations::current_time += calculations::deltaT;
    // cout<<calculations::deltaT<<endl;
    //  assert(calculations::deltaT<1);
}
