#include <iostream>
#include <vector>
#include "Vector.h"


template <class T>
void printVector(const VectorRobotech<T>& vec)
{
    std::cout << "X=" << vec.getX() << "\tY=" << vec.getY() << "\tZ=" << vec.getZ() << "\n";
    std::cout << "Magnitude=" << vec.magnitude() << "\n";
}


int main() {

    std::cout << "Constructors testing\n";
    std::cout << "Default constructor\n";
    VectorRobotech<double> defVector;
    printVector(defVector);

    std::cout << "Vector from 1D linear array\n";
    double arr[3]{2.0, 1.5, 3.5};
    VectorRobotech linearVector(arr);
    printVector(linearVector);

    std::cout << "Vector from std::vector\n";
    std::vector<double> vec{1.5, 2.5, 3.2};
    VectorRobotech vectorVector(vec);
    printVector(vectorVector);

    std::cout << "Copy constructor\n";
    VectorRobotech copyVector(linearVector);
    printVector(copyVector);

    std::cout << "Constructor from numbers\n";
    double x = 1.2;
    double y = 1.4;
    double z = 1.8;
    VectorRobotech numVector(x, y, z);
    printVector(numVector);

    std::cout << "\n#############################################\n";
    std::cout << "Test overload operators\n";
    std::cout << "Additions:\n";
    VectorRobotech v1(2.5, 3.0, 4.2);
    VectorRobotech v2(3.1, 4.2, 1.0);
    std::cout << "First vector\t";
    printVector(v1);
    std::cout << "Second vector\t";
    printVector(v2);

    return 0;
}
