#ifndef MyQGraphicsView_H
#define MyQGraphicsView_H

#include "Flow_drawer.h"

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
    string file1 = "/home/nikita/SPHSm6/txt/testB_3.txt";
    string file2 = "/home/nikita/SPHSm6/txt/testI_3.txt";
    QGraphicsScene *scene;
    Flow_Drawer *flow;
    //
    int counter_help = 0;
public:
    void draw();

    void clear();


};

#endif
