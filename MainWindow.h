#ifndef MainWindow_H
#define MainWindow_H

#include "myqgraphicsview.h"

class MainWindow : public QMainWindow {
Q_OBJECT
public:
    MainWindow(QWidget *parent = 0);

    ~MainWindow();


public slots:
    void start();
private:
    QPushButton* go;
    QTimer* timer;
    MyQGraphicsView *view;
    //menus
};

#endif
