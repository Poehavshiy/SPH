//
// Created by nikita on 15.06.16.
//

#include "Flow_drawer.h"
#include "myqgraphicsview.h"

Flow_Drawer::Flow_Drawer(const string &boundaryFile, const string &initFile) :
        Flow(boundaryFile, initFile) {
    calculator = new Calculator(s_distribution);
}

void Flow_Drawer::draw_boundary(QGraphicsScene *scene) {
    const vector<vector<Particle>> *data = reinterpret_cast<const vector<vector<Particle>> *>(this->get_part(0));
    double rad = 2;
    for (int i = 0; i < data->size(); ++i) {
        for (int j = 0; j < data->operator[](i).size(); ++j) {
            Particle *curent = const_cast<Particle *>(&data->operator[](i)[j]);
            scene->addEllipse(curent->X() - rad, curent->Y() - rad, rad * 2.0, rad * 2.0,
                              QPen(Qt::black), QBrush(Qt::SolidPattern));
        }
    }
}

void Flow_Drawer::draw_data(QGraphicsScene *scene) {
    const vector<Particle> *data = reinterpret_cast<const vector<Particle> *>(this->get_part(1));
    //   const vector<Particle> *data = (vector<Particle> *)this->get_part(1);
    double rad = 2;
    for (int i = 0; i < data->size(); ++i) {
        const Particle *curent = &data->operator[](i);
        scene->addEllipse(curent->X() - rad, curent->Y() - rad, rad * 2.0, rad * 2.0,
                          QPen(Qt::green), QBrush(Qt::SolidPattern));
    }
}

void Flow_Drawer::draw_grid(QGraphicsScene *scene) {
    const vector<vector<Cell>> *data = reinterpret_cast<const vector<vector<Cell>> *>(this->get_part(2));
    double rad = 2;
    for (int i = 0; i < data->size(); ++i) {
        for (int j = 0; j < data->operator[](i).size(); ++j) {
            const Cell *cur = &data->operator[](i)[j];
            const Point l_b = cur->get_left_bottom();
            const Point r_t = cur->get_right_top();
            //
            int colour[3] = {0, 0, 0};
            colour[1] = 128;
            QColor myC(colour[0], colour[1], colour[2]);
            QPointF pt1(l_b.x, l_b.y);
            QPointF pt2(l_b.x, r_t.y);
            QPointF pt3(r_t.x, r_t.y);
            QPointF pt4(r_t.x, l_b.y);
            //
            QLineF lt1(pt1, pt2);
            QLineF lt2(pt2, pt3);
            QLineF lt3(pt3, pt4);
            QLineF lt4(pt4, pt1);
            //
            scene->addLine(lt1, QPen(myC));
            scene->addLine(lt2, QPen(myC));
            scene->addLine(lt3, QPen(myC));
            scene->addLine(lt4, QPen(myC));

            //
        }
    }
}

void Flow_Drawer::draw_shadow(QGraphicsScene *scene) {
    const vector<vector<Cell>> *data = reinterpret_cast<const vector<vector<Cell>> *>(this->get_part(2));
    double rad = 2;
    int colour[3] = {0, 0, 0};
    for (int i = 0; i < data->size(); ++i) {
        for (int j = 0; j < data->operator[](i).size(); ++j) {
            const Cell *current = &data->operator[](i)[j];
            const vector<Particle> *curent_shadow = current->get_shadow();
            for (int i = 0; i < curent_shadow->size(); ++i) {
                scene->addEllipse(curent_shadow->operator[](i).X() - rad,
                                  curent_shadow->operator[](i).Y() - rad,
                                  rad * 2.0, rad * 2.0, QPen(QColor(colour[0], colour[1], colour[2])),
                                  QBrush(Qt::SolidPattern));
            }
        }
    }
}

void Flow_Drawer::draw_data_bycells(QGraphicsScene *scene) {
    const vector<vector<Cell>> *data = reinterpret_cast<const vector<vector<Cell>> *>(this->get_part(2));
    double rad = 2;
    int c[6][3] = {
            {255, 0,   0},
            {0,   255, 0},
            {0,   0,   255},
            {255,   255,   0},
            {0,   255,   255},
            {255,   0,   255}
    };
    QColor colors[6];
    for (int i = 0; i < 6; ++i) {
        colors[i].setRgb(c[i][0], c[i][1], c[i][2]);
    }
    int colour_counter = 0;
    for (int i = 0; i < data->size(); ++i) {
        for (int j = 0; j < data->operator[](i).size(); ++j) {
            if (colour_counter == 6) colour_counter = 0;
            const Cell *current = &data->operator[](i)[j];
            const vector<Particle *> *curent_real = current->get_real();
            for (int i = 0; i < curent_real->size(); ++i) {
                scene->addEllipse(curent_real->operator[](i)->X() - rad,
                                  curent_real->operator[](i)->Y() - rad,
                                  rad * 2.0, rad * 2.0, colors[colour_counter], QBrush(Qt::SolidPattern));
            }
            if(curent_real->size()!=0) ++colour_counter;
        }
    }
}

void Flow_Drawer::calculate_step(QGraphicsScene *scene) {
    draw_boundary(scene);
    draw_grid(scene);
    draw_shadow(scene);
    //
    calculator->calculate();
    draw_data_bycells(scene);
    // cout<<data.size();
}