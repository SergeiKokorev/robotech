#include <iostream>
#include <vector>

#include "Quaternion.h"
#include "Vector.h"


int main(){

    std::cout << "Quaternion default constructor\n";
    QRobotech<double>* q = new QRobotech();
    std::cout << "Constructed successfully!\n";

    return 0;
};
