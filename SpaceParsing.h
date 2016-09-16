//
// Created by nikita on 27.05.16.
//

#ifndef SPHSM6_SPACEPARSING_H
#define SPHSM6_SPACEPARSING_H

#include "Particile.h"


class Cell {
    vector<Particle *> real_group;
    vector<Particle *> boundary_group;
    vector<Particle> symetric_group;
    Point left_bottom;
    Point right_top;
    bool is_contain;
    bool boundary;
    vector<Particle> *boundary_line;
    static int C;//для дебажирования юзаю

    bool is_inside(Particle &particle) {
        if (
                (particle.X() >= left_bottom.x && particle.X() < right_top.x)
                && (particle.Y() >= left_bottom.y && particle.Y() < right_top.y)
                ) {
            return true;
        }
        return false;
    }

    //boundary_cell это ячейка, которая является ближайшей граничной
    void create_sim_part(Particle *target, Cell *boundary_cell) {
        // cout << Cell::C << endl;
        // Cell::C++;
        //Новой ячейке присвоили все кастарой
        Particle new_part = *target;
        new_part.set_vx(new_part.Vx() * -1);//а скорость наоборот
        new_part.set_vy(new_part.Vy() * -1);
        //Это все процесс вычисления позиции симетричной частицы за границей
        Point A = boundary_cell->boundary_group[0]->position();
        Point B = boundary_cell->boundary_group[1]->position();
        double Dx, Dy;
        Point C = target->position();
        if (A == B) {
            Dx = A.x;
            Dy = A.y;
        } else {
            Point vector_AC(C.x - A.x, C.y - A.y);
            Point vector_AB(B.x - A.x, B.y - A.y);
            //
            double v_AC_abs = sqrt(vector_AC.x * vector_AC.x + vector_AC.y * vector_AC.y);
            double v_AB_abs = sqrt(vector_AB.x * vector_AB.x + vector_AB.y * vector_AB.y);
            //
            double cos_ab_ac = ((vector_AC.x * vector_AB.x) + (vector_AC.y * vector_AB.y)) / (v_AC_abs * v_AB_abs);
            double v_AD_abs = v_AC_abs * cos_ab_ac;
            double scale = v_AD_abs / v_AB_abs;
            Dx = A.x + (B.x - A.x) * scale;
            Dy = A.y + (B.y - A.y) * scale;
        }
        //
        Point new_pos(Dx + (Dx - C.x), Dy + (Dy - C.y));
        new_part.set_pos(new_pos);
        boundary_cell->symetric_group.push_back(new_part);
    }
    //
public:
    friend class Calculator;

    friend class SpaceParsing;

    //получение углов квадрата
    Point get_left_bottom() const {
        return left_bottom;
    }

    //
    Point get_right_top() const {
        return right_top;
    }

    //попытка вставки частицы в this Cell
    bool insert(Particle &particle) {
        bool inside = is_inside(particle);//проверка вхождения частицы внутрь
        if (inside == true && particle.boundary_status() == true) {
            boundary_group.push_back(&particle);
            boundary = true;
        } else if (inside == true && particle.boundary_status() == false) {
            real_group.push_back(&particle);
        }
        return inside;
    }

    void set_points(Point &l_b, Point &r_t) {
        left_bottom = l_b;
        right_top = r_t;
    }

    bool is_boundary() {
        return boundary;
    }

    void add_part(Particle *p) {
        real_group.push_back(p);
    }

    void create_simetric(Cell *boundary_cell) {
        //   cout << Cell::C << endl;
        //   Cell::C++;
        for (int i = 0; i < real_group.size(); ++i) {
            create_sim_part(real_group[i], boundary_cell);
            //  cout << Cell::C << endl;
            //  Cell::C++;
        }
    }

    void remove_part(int &index) {
        swap(real_group.back(), real_group[index]);
        real_group.pop_back();
    }

    const vector<Particle> *get_shadow() const {
        return &symetric_group;
    }

    const vector<Particle *> *get_real() const {
        return &real_group;
    }

    //проверить на пустоту и если пуста, сделать пустой
    void emptiness() {
        if (real_group.size() != 0) is_contain = true;
        else is_contain = false;
    }
    //внешняя проверка на пустоту
    bool is_non_empty() {
        return is_contain;
    }

    friend std::istream &operator<<(ostream &is, Cell &income) {
        is << income.left_bottom;
        is << ' ' << 800 << endl;
        is << income.right_top;
        is << ' ' << 800 << endl;
        for (int i = 0; i < income.boundary_group.size(); ++i) {
            is << income.boundary_group[i]->X() << ' ' << income.boundary_group[i]->Y() << ' ' << 0 << endl;
        }
        //
        for (int i = 0; i < income.real_group.size(); ++i) {
            is << income.real_group[i]->X() << ' ' << income.real_group[i]->Y() << ' ' << 300 << endl;
        }
        //
        for (int i = 0; i < income.symetric_group.size(); ++i) {
            is << income.symetric_group[i].X() << ' ' << income.symetric_group[i].Y() << ' ' << 500 << endl;
        }
    }

    Cell(Point &l_b, Point &r_t) {
        left_bottom = l_b;
        right_top = r_t;
        boundary = false;
    }

    Cell() {
        boundary = false;
    }
};

class SpaceParsing {
    vector<Particle> *data_ptr;
    int cells_per_x;
    int cells_per_y;
    double x_size;
    double y_size;
    vector<vector<Cell>> part_groups;//groups of particles

    double *extrude_square(vector<vector<Particle>> &boundaries) {
        double low_x = 0, top_x = 0, low_y = 0, top_y = 0;
        for (int i = 0; i < boundaries.size(); ++i) {
            for (int j = 0; j < boundaries[i].size(); ++j) {
                if (boundaries[i][j].X() < low_x) low_x = boundaries[i][j].X();
                if (boundaries[i][j].X() > top_x) top_x = boundaries[i][j].X();
                if (boundaries[i][j].Y() < low_y) low_y = boundaries[i][j].Y();
                if (boundaries[i][j].Y() > top_y) top_y = boundaries[i][j].Y();

            }
        }
        double *square = new double[4];
        square[0] = low_x;
        square[1] = top_x;
        square[2] = low_y;
        square[3] = top_y;
        return square;
    }

    //
    void build_cells(vector<vector<Particle>> &boundaries) {
        double *square = extrude_square(boundaries);
        x_size = (square[1] - square[0]) / (cells_per_x - 1);
        y_size = (square[3] - square[2]) / (cells_per_y - 1);

        for (int i = 0; i < cells_per_y; ++i) {
            for (int j = 0; j < cells_per_x; ++j) {
                double cenral_x = square[0] + j * x_size;
                double cenral_y = square[2] + i * y_size;
                Point right_top(cenral_x + x_size / 2, cenral_y + y_size / 2);
                Point left_bottom(cenral_x - x_size / 2, cenral_y - y_size / 2);
                part_groups[i][j].set_points(left_bottom, right_top);
            }
        }
        delete[] square;
    }

    //
    bool insert(Particle &particle) {//insert particle in part_groups
        bool findStatus = false;
        for (int i = 0; i < part_groups.size(); ++i) {
            for (int j = 0; j < part_groups[i].size(); ++j) {
                findStatus = part_groups[i][j].insert(particle);
                if (findStatus == true) return findStatus;
            }
        }
        return findStatus;
    }

    //
    void distribute_boundaries(vector<vector<Particle>> &boundaries) {
        int k = 0;
        for (int i = 0; i < boundaries.size(); ++i) {
            for (int j = 0; j < boundaries[i].size(); ++j) {
                //for debuging. Should always be true, since  part_groups covers all space
                bool status = insert(boundaries[i][j]);
                //cout<<i<<' '<<j<<' '<<status<<' '<<k<<endl; ++k;
            }
        }
    }

    //
    void distribute_initial(vector<Particle> &data) {
        data_ptr = &data;
        for (int i = 0; i < data.size(); ++i) {
            Particle debug = data[i];
            bool status = insert(data[i]);
//            assert (status == true);
            //  cout<<i<<endl;
            if (i == 288) {
                int a = 0;
            }
        }
    }

    //
    Cell *is_near_boundary(int row, int column) {
        //если это ячейка НЕ ГРАНИЧНАЯ

        for (int i = 0; i < 8; i+=2) {
            Cell* current = get_cell_clockwise(row, column, i);
            if (current != nullptr && current->is_boundary())
                return current;
        }
        return nullptr;
    }

    //
    void create_symetric_groups() {
        Cell *boundary_cell = nullptr;
        for (int i = 0; i < part_groups.size(); ++i) {
            for (int j = 0; j < part_groups[i].size(); ++j) {
                boundary_cell = is_near_boundary(i, j);

                if (part_groups[i][j].is_boundary() == true) {
                    part_groups[i][j].create_simetric(&part_groups[i][j]);
                } else if (part_groups[i][j].is_boundary() == false && boundary_cell != nullptr) {
                    part_groups[i][j].create_simetric(boundary_cell);
                }

            }
        }
    }

    SpaceParsing(vector<vector<double>> &geometry, int per_x, int per_y) {
        cells_per_y = per_y;
        cells_per_x = per_x;
        part_groups.resize(cells_per_y);
        for (int i = 0; i < cells_per_y; ++i) {
            part_groups[i].resize(cells_per_x);
        }

    }

    //
    void set_cells_perax(int per_x, int per_y) {
        cells_per_x = per_x;
        cells_per_y = per_y;
    }

    //для перемещения частиц
    std::pair<int, int> find_around(const int &row, const int &column, Particle &target) {
        int r = row;
        int c = column;
        if (r - 1 >= 0 && part_groups[r - 1][c].is_inside(target) == true) {//case0
            std::pair<int, int> result(r - 1, c);
            return result;
        } else if (r - 1 >= 0 && c - 1 >= 0 && part_groups[r - 1][c - 1].is_inside(target) == true) {//case1
            std::pair<int, int> result(r - 1, c - 1);
            return result;
        } else if (c - 1 >= 0 && part_groups[r][c - 1].is_inside(target) == true) { //case2
            std::pair<int, int> result(r, c - 1);
            return result;
        } else if (c - 1 >= 0 && r + 1 < cells_per_y && part_groups[r + 1][c - 1].is_inside(target) == true) { //case3
            std::pair<int, int> result(r + 1, c - 1);
            return result;
        } else if (r + 1 < cells_per_y && part_groups[r + 1][c].is_inside(target) == true) {//case4
            std::pair<int, int> result(r + 1, c);
            return result;
        } else if (r + 1 < cells_per_y && c + 1 < cells_per_x &&
                   part_groups[r + 1][c + 1].is_inside(target) == true) {//case5
            std::pair<int, int> result(r + 1, c + 1);
            return result;
        } else if (c + 1 < cells_per_x && part_groups[r][c + 1].is_inside(target) == true) {//case6
            std::pair<int, int> result(r, c + 1);
            return result;
        } else if (r - 1 >= 0 && c + 1 < cells_per_x && part_groups[r - 1][c + 1].is_inside(target) == true) {
            std::pair<int, int> result(r - 1, c + 1);
            return result;
        }
    }

    //
    Cell* get_cell_clockwise(const int &row, const int &column, int n) {
        int r = row;
        int c = column;
        Cell* result;
        int res_r, res_c;
        //эта функция не должна отрабатывать для ГРАНИЧНЫХ ячеек
        switch (n) {
            case 0:
                res_r = r + 1;
                res_c = c;
                if(res_r ==  cells_per_y) return nullptr;
                return &part_groups[res_r][res_c];
            case (1):
                res_r = r + 1;
                res_c = c + 1;
                if(res_r == cells_per_y || res_c == cells_per_x) return nullptr;
                return &part_groups[res_r][res_c];
            case (2):
                res_r = r;
                res_c = c + 1;
                if(res_c == cells_per_x ) return nullptr;
                return &part_groups[res_r][res_c];
            case (3):
                res_r = r - 1;
                res_c = c + 1;
                if(res_r == -1 || res_c == cells_per_x) return nullptr;
                return &part_groups[res_r][res_c];
            case (4):
                res_r = r - 1;
                res_c = c;
                if(res_r == -1 )return nullptr;
                return &part_groups[res_r][res_c];
            case (5):
                res_r = r - 1;
                res_c = c - 1;
                if(res_r == -1 || res_c == -1) return nullptr;
                return &part_groups[res_r][res_c];
            case (6):
                res_r = r;
                res_c = c - 1;
                if( res_c == -1 ) return nullptr;
                return &part_groups[res_r][res_c];
            case (7):
                res_r = r + 1;
                res_c = c - 1;
                if( res_r == cells_per_y || res_c == -1 ) return nullptr;
                return &part_groups[res_r][res_c];
            default:

                return nullptr;
        }

    };


public:
    friend class Calculator;

    friend class Calculator_Drawer;

    static SpaceParsing *init
            (
                    vector<vector<double>> &geometry,
                    vector<vector<Particle>> &boundaries,
                    vector<Particle> &data,
                    int per_x,
                    int per_y
            ) {
        SpaceParsing *space = new SpaceParsing(geometry, per_x, per_y);
        space->build_cells(boundaries);
        space->distribute_boundaries(boundaries);
        space->distribute_initial(data);
        space->create_symetric_groups();

        return space;
    }

    //
    friend std::istream &operator<<(ostream &is, SpaceParsing &income) {
        //
        for (int i = 0; i < income.cells_per_y; ++i) {
            for (int j = 0; j < income.cells_per_x; ++j) {
                is << income.part_groups[i][j];
            }
        }
    }

    //
    const vector<vector<Cell>> *get_parsing() {
        return &part_groups;
    }
};


#endif //SPHSM6_SPACEPARSING_H
