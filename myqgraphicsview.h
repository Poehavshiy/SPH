#ifndef MyQGraphicsView_H
#define MyQGraphicsView_H

#include "Flow_drawer.h"
#include "Visualisator.h"

class MyQGraphicsView : public QGraphicsView {
Q_OBJECT
public:
    explicit MyQGraphicsView(QWidget *parent = 0);

    ~MyQGraphicsView();

signals:

public slots:

    //drawing
    void mousePressEvent(QMouseEvent *e);

    void ClearWindow();

    void mouseDoubleClickEvent(QMouseEvent *e);

private:
    string calc_case = "2";
    string file1 = "/home/nikita/SPH/txt/testB_" + calc_case + ".txt";
    string file2 = "/home/nikita/SPH/txt/testI_" +calc_case + ".txt";
    string path_tovis = "/home/nikita/SPHS/csv";
    Visualisator* vis;
    QGraphicsScene *scene;
    Flow_Drawer *flow;
    //
    int counter_help = 0;
public:
    void draw();

    void clear();


};

#endif
