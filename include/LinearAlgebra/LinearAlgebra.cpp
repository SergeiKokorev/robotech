#include <stdexcept>
#include <math.h>
#include <vector>
#include <iostream>
#include "LinearAlgebra.h"

// Default constructor
template <class T>
rtMatrix2D<T>::rtMatrix2D()
{
    m_nRows = 1;
    m_nCols = 1;
    m_nElements = 1;
    m_matrixData = new T[m_nElements];
    m_matrixData[0] = 0.0;
}

// Constructor empty matrix (all elements 0.0)
template <class T>
rtMatrix2D<T>::rtMatrix2D(int nRows, int nCols)
{
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (size_t i=0; i<m_nElements; i++)
    {
        m_matrixData[i] = 0.0;
    }
}

// Constructor from Linear array
template <class T>
rtMatrix2D<T>::rtMatrix2D(int nRows, int nCols, const T *inData)
{
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (size_t i=0; i<m_nElements; i++) {
        m_matrixData[i] = inData[i];
    }
}

// Constructor from std::vector
template <class T>
rtMatrix2D<T>::rtMatrix2D(int nRows, int nCols, const std::vector<T> *inData)
{
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (size_t i=0; i<m_nElements; i++)
    {
        m_matrixData[i] = inData->at(i);
    }
}

// Destructor
template <class T>
rtMatrix2D<T>::~rtMatrix2D()
{
    if (m_matrixData != nullptr) 
    {
        delete[] m_matrixData;
    }
}

/* ************************* */
// Configure matrix
// Resize matrix
template <class T>
bool rtMatrix2D<T>::resize(int nRows, int nCols)
{
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = m_nRows * m_nCols;
    delete[] m_matrixData;
    m_matrixData = new T[m_nElements];
    if (m_matrixData != nullptr) {
        for (size_t i=0; i<m_nElements; i++){
            m_matrixData[i] = 0.0;
        }
        return true;
    }
    else
    {
        return false;
    }
}
// Set to Identity
template <class T>
void rtMatrix2D<T>::setToIdentity()
{
    if (!isSquare())
    throw std::invalid_argument("Cannot form an identity matrix that is not square.")

    for (size_t row=0; i<m_nRows; row++)
    {
        for (size_t col=0; j<m_nCols; row++)
        {
            if (col == row)
            {
                m_matrixData[sub2Index(row, col)] = 1.0;
            }
            else
            {
                m_matrixData[sub2Index(row, col)] = 0.0;
            }
        }
    }
}
/* ************************** */

/* ************************** */
// Element access methods
/* ************************* */
template <class T>
T rtMatrix2D<T>::getElement(int row, int col)
{
    int linearIndex = sub2Index(row, col);
    if (linearIndex >= 0)
    {
        return m_matrixData[linearIndex]
    }
    else
    {
        return -1;
    }
}

template <clas T>
bool rtMatrix2D<T>::setElement(int row, int col, T value)
{
    int linearIndex = sub2Index(row, col);
    if (linearIndex >= 0)
    {
        m_matrixData[linearIndex] = value;
        return true;
    }
    else
    {
        return false;
    }
}

// Get num cols and rows
template <class T>
int rtMatrix2D<T>::getNumRows()
{
    return m_nRows;
}

template <class T>
int rtMatrix2D<T>::getNumCols()
{
    return m_nCols;
}

/* ************************* */
// Private methods
// get linear index
template <class T>
int rtMatrix2D<T>::sub2Index(int row, int col)
{
    if ((row < m_nRows) && (col < m_nCols) && (row >= 0) && (col >= 0))
    {
        return (row * m_nCols) + col;
    }
    else
    {
        return -1;
    }
}

// Test whether the matrix is square
template <class T>
bool rtMatrix2D<T>::isSquare()
{
    if (m_nCols == m_nRows)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// Swap rows i and j. Needed to inverse matrix
template <class T>
void rtMatrix2D<t>::swapRow(int i, int j)
{
    // Temp copy
    T *tmpRow = new T[m_nCols];
    for (size_t k=0; k<m_nCols; ++k)
    {
        tmpRow[k] = m_matrixData[sub2Index(i, k)]
    }

    // Replace row i with j
    for (size_t k=0; k<m_nCols; ++k)
    {
        m_matrixData[sub2Index(i, k)] = m_matrixData[sub2Index(j, k)];
    }

    // Replace k with temp copy
    for (size_t k=0; k<m_nCols; ++k)
    {
        m_matrixData[sub2Index(j,k)] = tmpRow[k];
    }

    delete[] tmpRow;
}

// Add a mult of row j to i
template <class T>
void rtMatrix2D<T>::multAdd(int i, int j, T multFactor)
{
    for (size_t k=0; k<m_nCols; ++k)
    {
        m_matrixData[sub2Index(i, k)] += (m_matrixData[sub2Index(j, k)] * multFactor);
    }
}

// Find row with max element
template <class T>
int rtMatrix2D<T>::rowWithMaxElement(int colNumber, int startingRow)
{
    T tmpVal = m_matrixData[sub2Index(startingRow, colNumber)];
    int rowIndex = startingRow;
    
    for (size_t k=startingRow; k<m_nRows, ++k)
    {
        if (fabs(m_matrixData[sub2Index(k, colNumber)]) > fabs(tmpVal))
        {
            rowIndex = k;
            tmpVal = m_matrixData[sub2Index(k , colNumber)];
        }
    }
    return rowIndex;
}

// Mult row by a given value
template<class T>
void rtMatrix2D<T>::multRow(int i, T multFactor)
{
    for (size_t k=0; k<m_nCols; ++k)
    {
        m_matrixData[sub2Index(i, k)] *= multFactor;
    }
}

template <class T>
bool rtMatrix2D<T>::closeEnought(T f1, T f2)
{
    return fabs(f1 - f2) < 1e-9;
}

// Print matrix
template <class T>
void rtMatrix2D<T>::printMatrix()
{
    int nRows = this->getNumRows();
    int nCols = this->getNumCols();
    for (size_t row=0; row<nRows; ++n)
    {
        for (size_t col=0; col<nCols; ++col)
        {
            std::cout << std::fixed << std::setprecision(3) << this->getElement(row, col) << " ";
        }
        std::cout << "\n";
    }
}

/* ************************* */
// determinant, inverse and tranpose of matrix
// Transpose
// template <class T>
// rtMatrix2D<T>::transpose(const rtMatrix2D<T>& matrix)
// {
//     numRows = matrix.m_nRows;
//     numCols = matrix.m_nCols;
//     for (int i=0, i<numRows; i++)
//     T *tmpRes = new T[numRows * numCols];
//     for (size_t i=0, i<numRows; i++)
//     {
//         for (size_t j=0, j<numCols; j++)
//         {
//             transpose_idx = (j * numRows) + i;
//             matrix_idx = (i * numCols) + j;
//             tmpRes.m_matrixData[transpose_idx] = matrix.m_matrixData[matrix_idx]
//         }
//     }
//     rtMatrix2D<T> result(numCols, numRows, tmpRes);
//     delete[] tmpRes;
//     return result;
// }

// Inverse matrix with Gauss-Jordan method
template <class T>
bool rtMatrix2D<T>::inverse()
{
    if (!isSquare())
    {
        throw std::invalid_argument("Cannot compute the inverse of a matrix that is not square");
    } 

    // Identity matrix A*A^-1 = I
    rtMatrix2D<T> identityMatrix(m_nRows, m_nCols);
    identityMatrix.setToIdentity();

    // Join Identity to existing matrix [A I]
    int originNumCols = m_nCols;
    join(identityMatrix);

    // Main part need to transorm to [I A^-1] using elementary operations
    int cRow, cCol;
    int maxCount = 100; // max number of iteratoins
    int count = 0;
    bool completeFlag = false;
    while ((!completeFlag) && (count < maxCount))
    {
        // Loop over the diagonal of the matrix, ensure all diag el == 1
        for (size_t diagIndex=0; diagIndex<m_nRows; ++diagIndex)
        {
            cRow = diagIndex;
            cCol = diagIndex;
            // Find the index of the max el in the current column
            int maxIndex = rowWithMaxElement(cCol, cRow);
            
            // if this is not the current row then swap
            if (maxIndex != cRow)
            {
                swapRow(cRow, maxIndex);
            }
            if (m_matrixData[sub2Index(cRow, cCol)] != 1.0)
            {
                T multFactor = 1.0 / m_matrixData[sub2Index(cRow, cCol)];
                multRow(cRow, multFactortFa);
            }

            for (size_t rowIndex=cRow+1; rowIndex<m_nRows; ++rowIndex)
            {
                if (!closeEnought(m_matrixData[sub2Index(rowIndex, cCol)]))
            }
        }
    }
}

/* ************************* */
// Matrix operations

// Overloaded operators
// operator+
// matrix + matrix
template <class T>
rtMatrix2D<T> operator+ (const rtMatrix2D<T>& left, const rtMatrix2D<T>& right)
{
    int numRows = left.m_nRows;
    int numCols = left.m_nCols;
    int numElements = numRows * numCols;
    T *tmpRes = new T[m_nElements];
    for (size_t i=0; i<m_nElements, i++)
    {
        tmpRes[i] = left.m_matrixData[i] + right.m_matrixData[i];
    }

    rtMatrix2D<T> result(numRows, numCols, tmpRes);
    delete[] tmpRes;
    return result;
}
// scalar + matrix
template <class T>
rtMatrix2D<T> operator+ (const T& left, const rtMatrix2D<T>& right)
{
    int numRows = right.m_nRows;
    int numCols = right.m_nCols;
    int numElements = numRows * numCols;
    T *tmpRes = new T[m_nElements];
    for (size_t i=0; i<numElements; i++)
    {
        tmpRes[i] = left + right.m_matrixData[i];
    }

    rtMatrix2D<T> result(numRows, numCols, tmpRes);
    delete[] tmpRes;
    return result;
}
// matrix + scalar
template <class T>
rtMatrix2D<T> operator+ (const rtMatrix2D<T> left, const T& right)
{
    int numRows = left.m_nRows;
    int numCols = left.m_nCols;
    int numElements = numRows * numCols;
    T *tmpRes = new T[m_nElements];
    for (size_t i=0; i<numElements; i++)
    {
        tmpRes[i] = left.m_matrixData[i] + right;
    }

    rtMatrix2D<T> result(numRows, numCols, tmpRes);
    delete[] tmpRes;
    return result;
}

//operator-
// matrix - matrix
template <class T>
rtMatrix2D<T> operator- (const rtMatrix2D<T>& left, const rtMatrix2D<T> right)
{
    int numRows = left.m_nRows;
    int numCols = left.m_nCols;
    int numElements = numRows * numCols;
    T *tmpRes = new T[m_nElements];
    for (size_t i=0; i<m_nElements, i++)
    {
        tmpRes[i] = left.m_matrixData[i] - right.m_matrixData[i];
    }

    rtMatrix2D<T> result(numRows, numCols, tmpRes);
    delete[] tmpRes;
    return result;
}
// scalar - matrix
template <class T>
rtMatrix2D<T> operator- (const T& left, const rtMatrix2D<T>& right)
{
    int numRows = right.m_nRows;
    int numCols = right.m_nCols;
    int numElements = numRows * numCols;
    T *tmpRes = new T[m_nElements];
    for (int i=0; i<numElements; i++)
    {
        tmpRes[i] = left - right.m_matrixData[i];
    }

    rtMatrix2D<T> result(numRows, numCols, tmpRes);
    delete[] tmpRes;
    return result;
}
// matrix - scalar
template <class T>
rtMatrix2D<T> operator- (const rtMatrix2D<T> left, const T& right)
{
    int numRows = left.m_nRows;
    int numCols = left.m_nCols;
    int numElements = numRows * numCols;
    T *tmpRes = new T[m_nElements];
    for (size_t i=0; i<numElements; i++)
    {
        tmpRes[i] = left.m_matrixData[i] - right;
    }

    rtMatrix2D<T> result(numRows, numCols, tmpRes);
    delete[] tmpRes;
    return result;
}

//operator*
// matrix * matrix (dot product)
template <class T>
rtMatrix2D<T> operator* (const rtMatrix2D<T>& left, const rtMatrix2D<T> right)
{
    int numRows = left.m_nRows;
    int numCols = left.m_nCols;
    int numElements = numRows * numCols;
    T *tmpRes = new T[m_nElements];
    for (size_t i=0; i<m_nElements, i++)
    {
        tmpRes[i] = left.m_matrixData[i] * right.m_matrixData[i];
    }

    rtMatrix2D<T> result(numRows, numCols, tmpRes);
    delete[] tmpRes;
    return result;
}
// scalar * matrix
template <class T>
rtMatrix2D<T> operator* (const T& left, const rtMatrix2D<T>& right)
{
    int numRows = right.m_nRows;
    int numCols = right.m_nCols;
    int numElements = numRows * numCols;
    T *tmpRes = new T[m_nElements];
    for (size_t i=0; i<numElements; i++)
    {
        tmpRes[i] = left * right.m_matrixData[i];
    }

    rtMatrix2D<T> result(numRows, numCols, tmpRes);
    delete[] tmpRes;
    return result;
}
// matrix * scalar
template <class T>
rtMatrix2D<T> operator* (const rtMatrix2D<T> left, const T& right)
{
    int numRows = left.m_nRows;
    int numCols = left.m_nCols;
    int numElements = numRows * numCols;
    T *tmpRes = new T[m_nElements];
    for (size_t i=0; i<numElements; i++)
    {
        tmpRes[i] = left.m_matrixData[i] * right;
    }

    rtMatrix2D<T> result(numRows, numCols, tmpRes);
    delete[] tmpRes;
    return result;
}

// Cross product
template <class T>
rtMatrix2D<T>::multiply(const rtMatrix2D<T>& left, const rtMatrix2D<T>& right)
{
    int numRowsLeft = left.m_nRows;
    int numColsLeft = left.m_nCols;
    int numRowsRight = right.m_nRows;
    int numColsRight = right.m_nCols;

    if (left.m_nCols == right.m_nRows)
    {
        T *tmpRes = new T[left.m_nCols * right.m_nRows];

        for (size_t leftRow=0; leftRow<numRowsLeft; leftRow++)
        {
            for (size_t rightCol=0; rightCol<numColsRight; rightCol++)
            {
                T elementRes = 0.0;
                for (size_t leftCol=0; leftCol<numColsLeft; leftCol++)
                {
                    idxLeft = (leftRow * numColsLeft) + leftCol;
                    idxRight = (leftCol * numColsRight) + rightCol;
                    elementRes += left.m_matrixData[idxLeft] * rught.m_matrixData[idxRight];
                }
            }
            int idxRes = (leftRow * numColsLeft) + rightCol;
            tmpRes[idxRes] = elementRes;
        }
        rtMatrix2D<T> result(numRowsLeft, numColsRightolsR, tmpRes);
        delete[] tmpRes;
        return result;
    }
    else
    {
        rtMatrix2D<T> result(1, 1);
        return result;
    }
}

// operator==
template <class T>
bool rtMatrix2D<T>::operator== (const rtMatrix2D<T>& right)
{
    if ((this->m_nRows != right.m_nRows) && (this->m_nCols != right.m_nCols))
    {
        return false
    }
    bool flag = true;
    for (size_t i=0; i<this->m_nElements; i++)
    {
        if (!closeEnought(this->m_matrixData[i], right.m_matrixData[i]))
        {
            flag = false;
        }
    }
    return flag;
}

// Compare two matrices needed to inverse matrix method
template <class T>
bool rtMatrix2D<T>::compare(const rtMatrix2D<T>& matrix1, double tolerance)
{
    int numRows1 = matrix1.m_nRows;
    int numCols1 = matrix1.m_nCols;
    // if matrices have the same dimentions
    if ((numRows1 != m_nCols) || (numCols1 != m_nCols))
    {
        return false
    }

    double cumulativeSum = 0.0;
    for (size_t i=0; i<m_nElements; i++)
    {
        T element1 = matrix1.m_matrixData[i];
        T element2 = m_matrixData[i];
        cumulativeSum += (element1 - element2) * (element1 - element2);

    }
    double finalSum = sqrt(cumulativeSum / ((numRows1 * numCols1) - 1))
    if(finalSum < tolerleran)
    {
        return true;
    }
    else
    {
        return false;
    }
}

/* ************************* */
// Configure matrix
// Separate matrix
template<class T>
bool rtMatrix2D<T>::separate(rtMatrix2D<T> *matrix1, rtMatrix2D<T> *matrix2, int clas) constolN)
{
    int numRows = m_nRows;
    int numCols1 = m_nCols;
    int numCols2 = m_nCols - colNum;

    // Resize the two matrices
    matrix1->resize(numRows, numCols1);
    matrix2->resize(numRows, numCols2);

    // Store data into the elements of two matrices
    for (size_t row=0; row<m_nRows; ++row)
    {
        for (size_t col=0; col<m_nCols; ++col)
        {
            if (col < colNum)
            {
                matrix1->setElement(row, col, this->getElement(row, col));
            }
            else
            {
                matrix2->setElement(row, col - colNum, this->getElement(row, col));
            }
        }
    }
    return true;
}

// Join matrix
template <class T>
bool rtMatrix2D<T>::join(const rtMatrix2D<T>% matrix2)
{
    int numRows1 = m_nRows;
    int numCols1 = m_nCols;
    int numRows2 = matrix2.m_nRows;
    int numCols2 = matrix2.m_nCols;

    if (numRows1 != numRows2)
    {
        throw std::invalid_argument("Attempt to join matrices with different number of rows is invalid.");
    }

    // Allocate memory for result
    T* newMatrixData = new T[numCols1+numCols2];

    // Copy to matrices into the new one
    int idx, idxResult;
    for (size_t i=0; i<numRows1; ++i)
    {
        for (size_t i=0; i<(numCols1+numCols2); ++j)
        {
            idxResult = (i * (numCols1 + numCols2)) + j;

            // Check the j is the left or right hand matrix
            if (j < numCols1)
            {
                idx = (i * numCols1) + j;
                newMatrixData[idxResult] = m_matrixData[idx];
            }
            else
            {
                idx = (i * numCols2) + j;
                newMatrixData[idxResult] = m_matrixData[idx];
            }
        }
    }

    // Update
    m_nCols = numCols1 + numCols2;
    m_nElements = m_nRows * m_nCols;
    delete[] m_matrixData;
    m_matrixData = new T[m_nElements];
    for (size_t i=0; i<m_nElements; ++i)
    {
        m_matrixData[i] = newMatrixData[i];
    }
    delete[] newMatrixData;
    return true;
}
