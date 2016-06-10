#include "myqgraphicsview.h"

//from QT
MyQGraphicsView::MyQGraphicsView(QWidget *parent) :
        QGraphicsView(parent) {
    //boudCounter=0;
    scene = new QGraphicsScene();
    this->setSceneRect(-100, -100, 600, 600);
    this->setScene(scene);
    flow = new Flow_Drawer(file1, file2);
}

//
MyQGraphicsView::~MyQGraphicsView() {
    delete flow;
    delete scene;
}

//
void MyQGraphicsView::clear() {
    scene->clear();
}
//

void MyQGraphicsView::mousePressEvent(QMouseEvent *e) {
    //flow->calculate(scene);

}

//
void MyQGraphicsView::mouseDoubleClickEvent(QMouseEvent *e) {

}

//
void MyQGraphicsView::draw() {
 /*   ++i;
    double rad = 20;
    scene->addEllipse(i - rad, i - rad, rad * 2.0, rad * 2.0,
                      QPen(Qt::green), QBrush(Qt::SolidPattern));*/
    if(calculations::current_time<1) flow->calculate_step(scene);

}

//
void  MyQGraphicsView::ClearWindow() {
    scene->clear();
}