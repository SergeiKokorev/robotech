#include <iostream>
#include "LinearAlgebra.h"


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


int main(){

    double A[12] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
    double *pA = A;
    double B[8] = {8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
    double *pB = B;
    double C[12] = {12.0, 11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
    double *pC = C;
    double E[9] = {2.0, 3.0, 1.0, 6.0, 8.0, 1.0, 3.0, 6.0, 7.0};
    double *pE = E;
    double scalar = 1.0;

    rtMatrix2D<double> matA(4, 3, pA);
    rtMatrix2D<double> matB(4, 2, pB);
    rtMatrix2D<double> matC(4, 3, pC);
    rtMatrix2D<double> matD(matA);
    rtMatrix2D<double> squareMat(3, 3, pE);

    /* Copy constructor */
    std::cout << "Copy constructor. matA has been copied to matD" << "\n";
    // Method printMatrix() is the private method and needs to be switched to public method
    printMatrix(matD);

    /* Test element access */
    std::cout << "############## Element access method ##########" << "\n";
    std::cout << "matA, element [0][0]\t" << matA.getElement(0, 0) << "\n";
    std::cout << "matA, element [0][1]\t" << matA.getElement(0, 1) << "\n";
    std::cout << "matA, element [0][2]\t" << matA.getElement(0, 2) << "\n";
    std::cout << "matA, element [1][0]\t" << matA.getElement(1, 0) << "\n";
    std::cout << "matA, element [1][1]\t" << matA.getElement(1, 1) << "\n";
    std::cout << "matA, element [1][2]\t" << matA.getElement(1, 2) << "\n";
    std::cout << "matA, element [2][0]\t" << matA.getElement(2, 0) << "\n";
    std::cout << "matA, element [2][1]\t" << matA.getElement(2, 1) << "\n";
    std::cout << "matA, element [2][2]\t" << matA.getElement(2, 2) << "\n";
    std::cout << "matA, element [3][0]\t" << matA.getElement(3, 0) << "\n";
    std::cout << "matA, element [3][1]\t" << matA.getElement(3, 1) << "\n";
    std::cout << "matA, element [3][2]\t" << matA.getElement(3, 2) << "\n";
    std::cout << "###############################################" << "\n";
    
    /* Test setter method */
    std::cout << "############## Element setter method ##########" << "\n";
    matA.setElement(0, 1, 3.0);
    std::cout << "Set element matA[0][1] to 3.0. Element equal\t" << matA.getElement(0, 1) << "\n";
    std::cout << "###############################################" << "\n";
    matA.setElement(0, 1, 2.0);

    /* Resize matrix*/
    std::cout << "############ Resize matrix matD ###############" << "\n";
    std::cout << "Initial matrix" << "\n";
    printMatrix(matD);
    std::cout << "Resized matrix" << "\n";
    matD.resize(3, 3);
    printMatrix(matD);

    /* Set matrix to Identity*/
    std::cout << "### Set matrix squareMat to Identity matrix ###" << "\n";
    std::cout << "Initial matrix" << "\n";
    printMatrix(squareMat);
    rtMatrix2D<double> copySquareMat(squareMat);
    squareMat.setToIdentity();
    std::cout << "Matrix after set to Identity matrix" << "\n";
    printMatrix(squareMat);

    /* Tets operators */
    std::cout << "################ Test operators ###############" << "\n";
    std::cout << "Initial matrices" << "\n";
    std::cout << "matA, matB and matC" << "\n";
    printMatrix(matA);
    printMatrix(matB);
    printMatrix(matC);
    rtMatrix2D<double> resMat(4, 3);
    resMat = matA + scalar;
    std::cout << "After right scalar addition\n";
    printMatrix(resMat);
    std::cout << "After left scalar addition\n";
    resMat = scalar + matA;
    printMatrix(resMat);
    std::cout << "Addition matA to matC. Result matrix:" << "\n";
    resMat = matA + matC;
    printMatrix(resMat);
    std::cout << "matA - scalar\n";
    resMat = matA - scalar;
    printMatrix(resMat);
    std::cout << "scalar - matA\n";
    resMat = scalar - matA;
    printMatrix(resMat);
    std::cout << "matA - matc\n";
    resMat = matA - matC;
    printMatrix(resMat);
    std::cout << "scalar * matA\n";
    resMat = scalar * matA;
    printMatrix(resMat);
    std::cout << "matA * scalar\n";
    resMat = matA * scalar;
    printMatrix(resMat);
    std::cout << "matA * matC (matrix multiplication)\n";
    resMat = matA * matC;
    printMatrix(resMat);
    std::cout << "matA.dot(matC) = scalar\n";
    double resScalar = matA.dot(matC);
    std::cout << resScalar << "\n";

    std::cout << "############# Inverse Matrix ################\n";
    std::cout << "Initial Matrix\n";
    printMatrix(copySquareMat);

    if (copySquareMat.inverse())
    {
        std::cout << "Inverse matrix\n";
        printMatrix(copySquareMat);
    }
    else
    {
        std::cout << "Matrix doen't has inverse.\n";
        printMatrix(copySquareMat);
    }

    return 0;
}
