#include "MainWindow.h"

MainWindow::MainWindow(QWidget *parent)
        : QMainWindow(parent) {

    setWindowTitle(tr("Something Title"));
    view = new MyQGraphicsView;//creating grafic scence inside of main window
    this->setGeometry(100, 100, 900, 900);
    this->setCentralWidget(view);
    go = new QPushButton("start", this);
    this->setMenuWidget(go);
   // connect(go, SIGNAL (released()), this, SLOT (start()));

    // основа рассчета

    timer = new QTimer;
    connect(timer, SIGNAL (timeout()), this, SLOT (start()));
    timer->start(100);

    //
}

void MainWindow::start() {
  //  timer = new QTimer;
 //   cout<<"go";
 //   timer -> start(1000* 20);
 //   while(timer->remainingTime() != 0) {
//        if (timer->remainingTime() % 1000 == 0) {
           // view->clear();
    view->clear();
    view->draw();
  //      }
 //   }
}

MainWindow::~MainWindow() {
}
