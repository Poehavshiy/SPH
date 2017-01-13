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

    vector<double> calc_distanse(Point target) {
        if (boundary_group.size() == 0) {
            int a = 0;
            throw 228;
        }
        double A = boundary_group[0]->Y() - boundary_group[1]->Y();
        double B = boundary_group[1]->X() - boundary_group[0]->X();
        double C = boundary_group[0]->X() * boundary_group[1]->Y() -
                   boundary_group[1]->X() * boundary_group[0]->Y();
        double distance = fabs((A * target.x + B * target.y + C) / sqrt(A * A + B * B));
        return vector<double>({A, B, C, distance});
    }


    //boundary_cell это ячейка, которая является ближайшей граничной
    void create_sim_part(Particle *target, Cell *boundary_cell) {

        vector<double> result = calc_distanse(target->position());
        double A = result[0];
        double B = result[1];
        double C = result[2];
        double D = -(A * target->Y() - B * target->X());
        double Y1 = -(D * A + C * B) / (B * B + A * A);
        //if (Double.IsInfinity(Y1)) Y1 = 0
        double X1 = (A * Y1 + D) / B;
        if (fabs(B) < 0.000000000001)
            X1 = -C / A;
        //if (Double.IsInfinity(X1)) X1 = 0;
        double X2 = 2 * X1 - target->X();
        double Y2 = 2 * Y1 - target->Y();
        Point new_pos(X2, Y2);
        //создание частицы
        Particle new_part = *target;
        new_part.set_vx(new_part.Vx() * -1);//а скорость наоборот
        new_part.set_vy(new_part.Vy() * -1);

        new_part.set_pos(new_pos);
        boundary_cell->symetric_group.push_back(new_part);

    }
    //

    //
public:
    friend class Calculator;

    friend class SpaceParsing;

    bool is_inside(Particle &particle) {
        if (
                (particle.X() >= left_bottom.x && particle.X() < right_top.x)
                && (particle.Y() >= left_bottom.y && particle.Y() < right_top.y)
                ) {
            return true;
        }
        return false;
    }

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
        //если вошла,и частица является граничной, силовой
        if (inside == true && particle.boundary_status() == true) {
            boundary_group.push_back(&particle);
            boundary = true;
        }
            //если входит и является реальной
        else if (inside == true && particle.boundary_status() == false) {
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

    void create_simetric(Cell *boundary_cell) {
        for (int i = 0; i < real_group.size(); ++i) {
            create_sim_part(real_group[i], boundary_cell);
        }
    }

    void remove_part(int &index) {
        swap(real_group[real_group.size() - 1], real_group[index]);
        real_group.pop_back();
    }

    const vector<Particle> *get_shadow() const {
        return &symetric_group;
    }

    const vector<Particle *> *get_real() const {
        return &real_group;
    }

    const vector<Particle *> *get_bound() const {
        return &boundary_group;
    }

    //проверить на пустоту и если пуста, сделать пустой
    void emptiness() {
        if (real_group.size() != 0) is_contain = true;
        else is_contain = false;
    }

    void clear_shadow() {
        symetric_group.clear();
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

struct Include_to {
    vector<Cell *> to;
    vector<int> for_including;
};

class SpaceParsing {
    vector<Particle> *data_ptr;
    int cells_per_x;
    int cells_per_y;
    double x_size;
    double y_size;
    vector<vector<Cell>> part_groups;//groups of particles
    vector<vector<Include_to>> replacement;

    double *extrude_square(vector<vector<Particle>> &boundaries) {
        double low_x = 100500, top_x = -1, low_y = 100500, top_y = -1;
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

    //Построение Cell ов
    void build_cells(vector<vector<Particle>> &boundaries) {
        double *square = extrude_square(boundaries);
        x_size = ((square[1] - square[0]) / ((double) cells_per_x)) * 1.0001;
        y_size = ((square[3] - square[2]) / ((double) cells_per_y)) * 1.0001;

        for (int i = 0; i < cells_per_y; ++i) {
            for (int j = 0; j < cells_per_x; ++j) {
                double cenral_x = square[0] + j * x_size + x_size / 2.0;
                double cenral_y = square[2] + i * y_size + y_size / 2.0;
                Point right_top(cenral_x + x_size / 2.0, cenral_y + y_size / 2.0);
                Point left_bottom(cenral_x - x_size / 2.0, cenral_y - y_size / 2.0);
                part_groups[i][j].set_points(left_bottom, right_top);
            }
        }
        delete[] square;
    }

    //начальная вставка частицы для построения начального распределения
    bool insert(Particle &particle) {
        bool findStatus = false;
        for (int i = 0; i < part_groups.size(); ++i) {
            for (int j = 0; j < part_groups[i].size(); ++j) {
                findStatus = part_groups[i][j].insert(particle);
                if (findStatus == true) {
                    return findStatus;
                }
            }
        }
        return findStatus;
    }

    //Рспределение граничных частиц (те, что дают силу)
    void distribute_boundaries(vector<vector<Particle>> &boundaries) {
        int k = 0;
        for (int i = 0; i < boundaries.size(); ++i) {
            for (int j = 0; j < boundaries[i].size(); ++j) {
                //for debuging. Should always be true, since  part_groups vector covers all space
                bool status = insert(boundaries[i][j]);
                //cout<<i<<' '<<j<<' '<<status<<' '<<k<<endl; ++k;
            }
        }
        int a = 1;
    }

    //распределим начальные частицы
    void distribute_initial(vector<Particle> &data) {
        data_ptr = &data;
        for (int i = 0; i < data.size(); ++i) {
            bool status = insert(data[i]);
//            assert (status == true);
            //  cout<<i<<endl;
        }
    }

    //
    Cell *is_near_boundary(int row, int column) {
        //если это ячейка НЕ ГРАНИЧНАЯ

        for (int i = 0; i < 8; i += 2) {
            Cell *current = get_cell_clockwise(row, column, i);
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
                if (i == 0 && j == 0) {
                    int a = 1;
                }
                if (part_groups[i][j].is_boundary() == true) {
                    part_groups[i][j].create_simetric(&part_groups[i][j]);
                } else if (part_groups[i][j].is_boundary() == false && boundary_cell != nullptr) {
                    part_groups[i][j].create_simetric(boundary_cell);
                }

            }
        }
    }

    //конструктор, который тупо объявляет двумерный масиив Cell
    SpaceParsing(vector<vector<double>> &geometry, int per_x, int per_y) {
        cells_per_y = per_y;
        cells_per_x = per_x;
        part_groups.resize(cells_per_y);
        replacement.resize(cells_per_y);
        for (int i = 0; i < cells_per_y; ++i) {
            part_groups[i].resize(cells_per_x);
            replacement[i].resize(cells_per_x);
        }

    }

    //для перемещения частиц- проход по 8-окрестности row column и поиск, куда там входит target
    Cell *find_around(const int &row, const int &column, Particle &target) {

        for (int i = 0; i < 8; ++i) {
            Cell *current = get_cell_clockwise(row, column, i);
            if (current != nullptr)
                if (current->is_inside(target))
                    return current;
        }
        return nullptr;
    }

    //
    Cell *get_cell_clockwise(const int &row, const int &column, int n) {
        int r = row;
        int c = column;
        Cell *result;
        int res_r, res_c;
        //эта функция не должна отрабатывать для ГРАНИЧНЫХ ячеек
        switch (n) {
            case 0:
                res_r = r + 1;
                res_c = c;
                if (res_r == cells_per_y) return nullptr;
                return &part_groups[res_r][res_c];
            case (1):
                res_r = r + 1;
                res_c = c + 1;
                if (res_r == cells_per_y || res_c == cells_per_x) return nullptr;
                return &part_groups[res_r][res_c];
            case (2):
                res_r = r;
                res_c = c + 1;
                if (res_c == cells_per_x) return nullptr;
                return &part_groups[res_r][res_c];
            case (3):
                res_r = r - 1;
                res_c = c + 1;
                if (res_r == -1 || res_c == cells_per_x) return nullptr;
                return &part_groups[res_r][res_c];
            case (4):
                res_r = r - 1;
                res_c = c;
                if (res_r == -1)return nullptr;
                return &part_groups[res_r][res_c];
            case (5):
                res_r = r - 1;
                res_c = c - 1;
                if (res_r == -1 || res_c == -1) return nullptr;
                return &part_groups[res_r][res_c];
            case (6):
                res_r = r;
                res_c = c - 1;
                if (res_c == -1) return nullptr;
                return &part_groups[res_r][res_c];
            case (7):
                res_r = r + 1;
                res_c = c - 1;
                if (res_r == cells_per_y || res_c == -1) return nullptr;
                return &part_groups[res_r][res_c];
            default:

                return nullptr;
        }

    };

    void fill_for_ij(int &row, int &column, Include_to &for_ij) {

        vector<Particle *> &group_ref = part_groups[row][column].real_group;
        int size = group_ref.size();
        for (int i = 0; i < size; ++i) {
            if (!part_groups[row][column].is_inside(*group_ref[i])) {
                Cell *to = find_around(row, column, *group_ref[i]);
                if (to != nullptr) {
                    for_ij.to.push_back(to);
                    for_ij.for_including.push_back(i);
                }
            }
        }
    }

    void fill_replacement() {

        for (int i = 0; i < cells_per_y; ++i) {
            for (int j = 0; j < cells_per_x; ++j) {
                Include_to for_ij;
                fill_for_ij(i, j, for_ij);
                replacement[i][j].for_including = std::move(for_ij.for_including);
                replacement[i][j].to = std::move(for_ij.to);
            }
        }
    }

    //
    void acctualy_replace() {
        //сначала добавим
        for (int i = 0; i < cells_per_y; ++i) {
            for (int j = 0; j < cells_per_x; ++j) {
                //если есть что переместить из этой ячейки
                if (replacement[i][j].to.size() != 0) {
                    Include_to &current = replacement[i][j];
                    for (int k = 0; k < current.to.size(); ++k) {
                        //переместить из i,j Cell частичку с к-атой позиции в k-й Cell
                        current.to[k]->real_group.push_back(part_groups[i][j].real_group[current.for_including[k]]);
                    }
                }
            }
        }
        //теперь удалим ie индексы
        for (int i = 0; i < cells_per_y; ++i) {
            for (int j = 0; j < cells_per_x; ++j) {
                if (replacement[i][j].to.size() != 0) {

                    Include_to &current = replacement[i][j];
                    for (int k = current.to.size() - 1; k >= 0; --k) {
                        //удалим  из i,j Cell частички с к-атой позиций
                        part_groups[i][j].remove_part(current.for_including[k]);
                        int qwerty = 1;
                    }
                }
            }
        }
    }

    void clear_symetric_groups() {
        //теперь удалим ie индексы
        for (int i = 0; i < cells_per_y; ++i) {
            for (int j = 0; j < cells_per_x; ++j) {
                part_groups[i][j].clear_shadow();
            }
        }
    }

public:
    friend class Calculator;

    friend class Calculator_Drawer;

    double h;

    void replace() {
        fill_replacement();
        acctualy_replace();
    }

    static SpaceParsing *init
            (
                    vector<vector<double>> &geometry,
                    vector<vector<Particle>> &boundaries,
                    vector<Particle> &data,
                    int per_x,
                    int per_y,
                    double h
            ) {
        SpaceParsing *space = new SpaceParsing(geometry, per_x, per_y);
        space->build_cells(boundaries);
        space->distribute_boundaries(boundaries);
        space->distribute_initial(data);
        space->h = h;
        return space;
    }

    //
    const vector<vector<Cell>> *get_parsing() {
        return &part_groups;
    }


};


#endif //SPHSM6_SPACEPARSING_H
