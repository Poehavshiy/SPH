//
// Created by nikita on 15.06.16.
//

#include "Flow_drawer.h"
#include "myqgraphicsview.h"

Flow_Drawer::COLOUR Flow_Drawer::GetColour(double v, double vmin, double vmax) {
    COLOUR c = {1.0, 1.0, 1.0}; // white
    double dv;

    if (v < vmin)
        v = vmin;
    if (v > vmax)
        v = vmax;
    dv = vmax - vmin;

    if (v < (vmin + 0.25 * dv)) {
        c.r = 0;
        c.g = 4 * (v - vmin) / dv;
    } else if (v < (vmin + 0.5 * dv)) {
        c.r = 0;
        c.b = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
    } else if (v < (vmin + 0.75 * dv)) {
        c.r = 4 * (v - vmin - 0.5 * dv) / dv;
        c.b = 0;
    } else {
        c.g = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
        c.b = 0;
    }

    return (c);
}


Flow_Drawer::Flow_Drawer(const string &boundaryFile, const string &initFile) :
        Flow(boundaryFile, initFile),
        frames({50, 100, 150, 200}), times({0.025 / 1000, 0.05 / 1000, 0.075 / 1000, 0.1 / 1000}) {
    calculator = new Calculator(s_distribution, smooth_length);
    /* cof  = 3;
     cof1 = 200;*/
    cof = 1;
    cofx = 0;
    cofy = 0;

};

void Flow_Drawer::draw_boundary(QGraphicsScene *scene) {
    const vector<vector<Particle>> *draw = reinterpret_cast<const vector<vector<Particle>> *>(this->get_part(0));
    double rad = 2;
    for (int i = 0; i < draw->size(); ++i) {
        for (int j = 0; j < draw->operator[](i).size(); ++j) {
            Particle *curent = const_cast<Particle *>(&draw->operator[](i)[j]);
            scene->addEllipse(curent->X() * cof + cofx - rad, curent->Y() * cof - rad + cofy, rad * 2.0, rad * 2.0,
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
            QPointF pt1(l_b.x * cof + cofx, l_b.y * cof + cofy);
            QPointF pt2(l_b.x * cof + cofx, r_t.y * cof + cofy);
            QPointF pt3(r_t.x * cof + cofx, r_t.y * cof + cofy);
            QPointF pt4(r_t.x * cof + cofx, l_b.y * cof + cofy);
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
    const vector<vector<Cell>> *draw = reinterpret_cast<const vector<vector<Cell>> *>(this->get_part(2));
    double rad = 2;
    for (int i = 0; i < draw->size(); ++i) {
        for (int j = 0; j < draw->operator[](i).size(); ++j) {
            const Cell *current = &draw->operator[](i)[j];
            const vector<Particle> *curent_shadow = current->get_shadow();
            for (int k = 0; k < curent_shadow->size(); ++k) {
                double value = curent_shadow->operator[](k).P();

                COLOUR currentC = GetColour(value, 0, maxP);
                QColor QTcurrentC;
                QTcurrentC.setRgb(255 * currentC.r, 255 * currentC.g, 255 * currentC.b);
                scene->addEllipse(curent_shadow->operator[](k).X() * cof + cofx - rad,
                                  curent_shadow->operator[](k).Y() * cof + cofy - rad,
                                  rad * 2.0, rad * 2.0, QPen(Qt::black), QBrush(Qt::black));
            }
        }
    }
}

void Flow_Drawer::draw_data_bycells(QGraphicsScene *scene) {
    const vector<vector<Cell>> *data = reinterpret_cast<const vector<vector<Cell>> *>(this->get_part(2));
    double rad = 2;
    int c[8][3] = {
            {255, 0,   0},
            {0,   255, 0},
            {0,   0,   255},
            {255, 255, 0},
            {0,   255, 255},
            {255, 0,   255},
            {128, 128, 0},
            {128, 0,   128}
    };
    QColor colors[8];
    for (int i = 0; i < 8; ++i) {
        colors[i].setRgb(c[i][0], c[i][1], c[i][2]);
    }
    int colour_counter = 0;
    for (int i = 0; i < data->size(); ++i) {
        for (int j = 0; j < data->operator[](i).size(); ++j) {
            if (colour_counter == 8) colour_counter = 0;
            const Cell *current = &data->operator[](i)[j];
            const vector<Particle *> *curent_real = current->get_real();
            for (int i = 0; i < curent_real->size(); ++i) {
                scene->addEllipse(curent_real->operator[](i)->X() - rad,
                                  curent_real->operator[](i)->Y() - rad,
                                  rad * 2.0, rad * 2.0, colors[colour_counter], QBrush(Qt::SolidPattern));
            }
            if (curent_real->size() != 0) ++colour_counter;
        }
    }
}

void Flow_Drawer::draw_dataP(QGraphicsScene *scene) {
    double rad = 1;
    for (int i = 0; i < data.size(); ++i) {
        double value = data[i].P();
        COLOUR currentC = GetColour(value, 0, calculator->get_maxP());
        QColor QTcurrentC;
        QTcurrentC.setRgb(255 * currentC.r, 255 * currentC.g, 255 * currentC.b);
        scene->addEllipse(data[i].X() * cof - rad + cofx,
                          data[i].Y() * cof - rad + cofy,
                          rad * 2.0, rad * 2.0, QTcurrentC, QBrush(Qt::SolidPattern));
    }
}

void Flow_Drawer::show_information(QGraphicsScene *scene) {
    set_text("iteration=", scene, pair<int, int>(-50, -150), calculator->get_iteration());
    set_text("time=", scene, pair<int, int>(-50, -130), calculator->get_time());
    set_text("maxP=", scene, pair<int, int>(-50, -110), calculator->get_maxP());
    set_text("minP=", scene, pair<int, int>(-50, -90), calculator->get_minP());
    set_text("maxV=", scene, pair<int, int>(50, -150), calculator->get_maxV());
    set_text("theoreticalV=", scene, pair<int, int>(50, -130), calculator->get_V_theory());
}

void Flow_Drawer::set_text(QString &&text, QGraphicsScene *scene, pair<int, int> &&position, double &&value) {
    QString test_to_show = text + QString::number(value);
    QGraphicsTextItem *text_item = scene->addText(test_to_show);
    text_item->setPos(position.first, position.second);
}

void Flow_Drawer::write_py_data(int target_itertion) {
    std::ofstream outfile(mathplot_path + "mathplot" + to_string(target_itertion) + ".txt");
    outfile << "begin.\n";
    for (int i = 0; i < data.size(); ++i) {
        double x = data[i].X();
        double y = data[i].Y();
        double value = data[i].P();
        outfile << x << " " << y << " " << value << '\n';
    }
    outfile << "end.\n";
}

void Flow_Drawer::calculate_step(QGraphicsScene *scene) {
    calculator->calculate();
    draw_dataP(scene);
    //draw_data_bycells(scene);

    //  draw_boundary(scene);
    draw_grid(scene);
    // draw_shadow(scene);
    //

    show_information(scene);

    if (calculator->get_iteration() == frames[python_iter]) {
        ++python_iter;
        write_py_data(python_iter);
    }



    // cout<<data.size();
}