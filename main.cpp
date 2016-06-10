<<<<<<< HEAD
#include <QApplication>
#include "MainWindow.h"

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);

	MainWindow myW;
	myW.show();
	return a.exec();
}
=======
#include "Flow.h"


int main() {

    string file1 = "/home/nikita/SPHSm6/testB_1.txt";
    string file2 = "/home/nikita/SPHSm6/testI_1.txt";
    Flow A(file1, file2);
    A.calculate();

    return 0;
}
>>>>>>> 883cc5708c435c9115b864d6fdf26606dffc1703
