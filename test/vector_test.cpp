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
    double arr[3] = {2.0, 1.5, 3.5};
    double *pArr = arr;
    VectorRobotech<double> linearVector(pArr);
    pArr = nullptr;
    delete[] pArr;
    printVector(linearVector);
    linearVector.norm();
    std::cout << "Normalized vector:\n";
    printVector(linearVector);

    std::cout << "Vector from std::vector\n";
    std::vector<double> *vec = new std::vector<double> {1.2, 3.2, 4.5};
    VectorRobotech<double> vectorVector(vec);
    printVector(vectorVector);

    std::cout << "Copy constructor\n";
    VectorRobotech<double> copyVector(linearVector);
    printVector(copyVector);

    std::cout << "Constructor from numbers\n";
    double x = 1.2;
    double y = 1.4;
    double z = 1.8;
    VectorRobotech<double> numVector(x, y, z);
    printVector(numVector);

    std::cout << "\n#############################################\n";
    std::cout << "Test setter methods\n";
    VectorRobotech<double> vs1(1.2, 3.2, 1.1);
    std::cout << "Initial vector\n";
    printVector(vs1);
    std::cout << "setX(4.0), setY(1.5), setZ(2.5)\n";
    vs1.setX(4.0).setY(1.5).setZ(2.5);
    printVector(vs1);

    std::cout << "\n#############################################\n";
    std::cout << "Test overload operators\n";
    std::cout << "Additions:\n";
    VectorRobotech<double> v1(2.5, 3.0, 4.2);
    VectorRobotech<double> v2(3.1, 4.2, 1.0);
    std::cout << "First vector\t";
    printVector(v1);
    std::cout << "Second vector\t";
    printVector(v2);
    printVector(v1 + v2);

    std::cout << "Substraction:\n";
    std::cout << "First - Second\n";
    printVector(v1 - v2);
    std::cout << "Second - First\n";
    printVector(v2 - v1);

    std::cout << "Cross product\n";
    std::cout << "First *  Second\n";
    printVector(v1 * v2);
        std::cout << "Second *  First\n";
    printVector(v2 * v1);

    return 0;
}
