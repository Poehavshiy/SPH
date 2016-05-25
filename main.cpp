#include "Flow.h"

int main() {

    string file1 = "/home/nikita/SPHSm6/testB.txt";
    string file2 = "/home/nikita/SPHSm6/testI.txt";
    Flow A(file1, file2);
    A.calculate();
    return 0;
}