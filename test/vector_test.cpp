#include <iostream>
#include <vector>
#include "Vector.h"
#include "LinearAlgebra.h"


template <class T>
void printVector(const VectorRobotech<T>& vec)
{
    std::cout << "X=" << vec.getX() << "\tY=" << vec.getY() << "\tZ=" << vec.getZ() << "\n";
    std::cout << "Magnitude=" << vec.magnitude() << "\n";
}

template <class T>
void printMatrix(rtMatrix2D<T>& m)
{
    int nRows = m.getNumRows();
    int nCols = m.getNumCols();

    for (int row=0; row<nRows; ++row)
    {
        for (int col=0; col<nCols; ++col)
        {
            std::cout << m.getElement(row, col) << " ";
        }
        std::cout << "\n";
    }
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

    std::cout << "Outer product\n";
    std::cout << "First vector v1\n";
    printVector(v1);
    std::cout << "Second vector v2\n";
    printVector(v2);
    std::cout << "v1.outer(v2)\n";
    rtMatrix2D<double> res = v1.outer(v2);
    printMatrix(res);

    std::cout << "v2.outer(v1)\n";
    res = v2.outer(v1);
    printMatrix(res);

    std::cout << "Wedge product of two vectors\n";
    std::cout << "The first one v1\n";
    printVector(v1);
    std::cout << "The second one v2\n";
    printVector(v2);
    std::cout << "v1.wedge(v2)\n";
    res = v1.wedge(v2);
    printMatrix(res);
    std::cout << "v2.wedge(v1)\n";
    res = v2.wedge(v1);
    printMatrix(res);

    return 0;
}
