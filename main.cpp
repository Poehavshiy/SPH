#include "Flow.h"


int main() {

    string file1 = "/home/nikita/SPHSm6/testB_1.txt";
    string file2 = "/home/nikita/SPHSm6/testI_1.txt";
    Flow A(file1, file2);
    A.calculate();

    return 0;
}