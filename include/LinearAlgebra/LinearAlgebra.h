#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#include <stdexcept>
#include <math.h>
#include <vector>

template <class T>
class rtMatrix2D
{
    public:
        // Definthe various constructor
        rtMatrix2D();
        rtMatrix2D(int nRows, int nCols);
        // Constructor from linear 1D array
        rtMatrix2D(int nRows, int nCols, const T *inData);
        // Copy cstructor
        rtMatrix2D(rtMatrix2D<T> const& inMatrix);
        // Constructor from vector
        rtMatrix2D(int nRows, int nCols, const std::vector<T> *inData);

        // Destructor
        ~rtMatrix2D() = default;

        // Matrix cross product
        rtMatrix2D<T> multiply(const rtMatrix2D<T>& right) const;

        // Iverse matrix

        // Manipulate methods
        // Transpose matrix and inverse matrix
        bool transpose();
        bool inverse();

        // Configure methods
        bool resize(int numRows, int numCols);
        void setToIdentity();

        // Element access
        T getElement(int row, int col);
        bool setElement(int row, int col, T value);
        int getNumRows() const;
        int getNumCols() const;

        // Overload 
        bool operator== (const rtMatrix2D<T>& right);
        bool compare(const rtMatrix2D<T>& matrix1, double tolerance);

        // Overload +, - and * operators
        template <class U> friend rtMatrix2D<U> operator+ (const rtMatrix2D<U>& left, const rtMatrix2D<U>& right);
        template <class U> friend rtMatrix2D<U> operator+ (const U& left, const rtMatrix2D<U>& right);
        template <class U> friend rtMatrix2D<U> operator+ (const rtMatrix2D<U>& left, const U& right);

        template <class U> friend rtMatrix2D<U> operator- (const rtMatrix2D<U>& left, const rtMatrix2D<U>& right);
        template <class U> friend rtMatrix2D<U> operator- (const U& left, const rtMatrix2D<U>& right);
        template <class U> friend rtMatrix2D<U> operator- (const rtMatrix2D<U>& left, const U& right);

        template <class U> friend rtMatrix2D<U> operator* (const rtMatrix2D<U>& left, const rtMatrix2D<U>& right);
        template <class U> friend rtMatrix2D<U> operator* (const U& left, const rtMatrix2D<U>& right);
        template <class U> friend rtMatrix2D<U> operator* (const rtMatrix2D<U>& left, const U& right);

        T dot(const rtMatrix2D<T>& rhs);

    // Temporary switch to public for testing matrix
    private:
        int sub2Index(int row, int col);

        // Needed to find an inverse Matrix
        bool separate(rtMatrix2D<T> *matrix1, rtMatrix2D<T> *matrix2, int colNum);
        bool isSquare();
        bool closeEnought(T f1, T f2);
        void swapRow(int i, int j);
        void multAdd(int i, int j, T multFactor);
        void multRow(int i, T multFactor);
        bool join(const rtMatrix2D<T>& matrix2);
        int rowWithMaxElement(int colNumber, int startingRow);
        void printMatrix();

    private:
        T *m_matrixData;
        int m_nRows, m_nCols, m_nElements;

};

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
    for (size_t i=0; i<m_nElements; i++) 
    {
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

// Copy constructor
template <class T>
rtMatrix2D<T>::rtMatrix2D(rtMatrix2D<T> const& inMatrix)
{
    m_nRows = inMatrix.m_nRows;
    m_nCols = inMatrix.m_nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (size_t i=0; i<m_nElements; ++i)
    {
        m_matrixData[i] = inMatrix.m_matrixData[i];
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
        for (size_t i=0; i<m_nElements; ++i){
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
        { throw std::invalid_argument("Cannot form an identity matrix that is not square."); }

    for (size_t row=0; row<m_nRows; ++row)
    {
        for (size_t col=0; col<m_nCols; ++col)
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
        return m_matrixData[linearIndex];
    }
    else
    {
        return -1;
    }
}

// Setter method
template <class T>
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
int rtMatrix2D<T>::getNumRows() const
{
    return m_nRows;
}

template <class T>
int rtMatrix2D<T>::getNumCols() const
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

// Join two matrices
template <class T>
bool rtMatrix2D<T>::join(const rtMatrix2D<T>& matrix2)
{
    int numRows1 = m_nRows;
    int numRows2 =matrix2.m_nRows;
    int numCols1 = m_nCols;
    int numCols2 = matrix2.m_nCols;

    if (numRows1 != numRows2)
    {
        return false;
    }
    T* newMatrixData = new T[numRows1*(numCols1+numCols2)];

    int linearIndex, resultLinearIndex;
    for (size_t i=0; i<numRows1; ++i)
    {
        for (size_t j=0; j<(numCols1+numCols2); ++j)
        {
            resultLinearIndex = (i * (numCols1+numCols2)) + j;

            if (j < numCols1)
            {
                linearIndex = (i * numCols1) + j;
                newMatrixData[resultLinearIndex] = m_matrixData[linearIndex];
            }
            else
            {
                linearIndex = (i * numCols2) + (j - numCols1);
                newMatrixData[resultLinearIndex] = matrix2.m_matrixData[linearIndex];
            }
        }
    }
    m_nCols = numCols1 + numCols2;
    m_nElements = m_nRows * m_nCols;
    delete[] m_matrixData;
    m_matrixData = new T[m_nElements];
    for (size_t i=0; i<m_nElements; ++i)
    {
        m_matrixData[i] = newMatrixData[i];
    }
    return true;
}

// Swap rows i and j. Needed to inverse matrix
template <class T>
void rtMatrix2D<T>::swapRow(int i, int j)
{
    // Temp copy
    T *tmpRow = new T[m_nCols];
    for (size_t k=0; k<m_nCols; ++k)
    {
        tmpRow[k] = m_matrixData[sub2Index(i, k)];
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
    
    for (size_t k=startingRow; k<m_nRows; ++k)
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
    return fabs(f1 - f2) < 1e-6;
}

template <class T>
bool rtMatrix2D<T>::separate(rtMatrix2D<T>* matrix1, rtMatrix2D<T>* matrix2, int colNum)
{
    int numRows = m_nRows;
    int numCols1 = colNum;
    int numCols2 = m_nCols - colNum;

    matrix1->resize(numRows, numCols1);
    matrix2->resize(numRows, numCols2);
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
                matrix2->setElement(row, col-colNum, this->getElement(row, col));
            }
        }
    }
    return true;
}

// Print matrix
template <class T>
void rtMatrix2D<T>::printMatrix()
{
    int nRows = this->getNumRows();
    int nCols = this->getNumCols();

    for (size_t row=0; row<nRows; ++row)
    {
        for (size_t col=0; col<nCols; ++col)
        {
            std::cout << this->getElement(row, col) << " ";
        }
        std::cout << "\n";
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
    T *tmpRes = new T[numElements];
    for (size_t i=0; i<numElements; i++)
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

    T *tmpRes = new T[numElements];
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
rtMatrix2D<T> operator+ (const rtMatrix2D<T>& left, const T& right)
{
    int numRows = left.m_nRows;
    int numCols = left.m_nCols;
    int numElements = numRows * numCols;
    T *tmpRes = new T[numElements];
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
rtMatrix2D<T> operator- (const rtMatrix2D<T>& left, const rtMatrix2D<T>& right)
{
    int numRows = left.m_nRows;
    int numCols = left.m_nCols;
    int numElements = numRows * numCols;
    T *tmpRes = new T[numElements];
    for (size_t i=0; i<numElements; ++i)
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
    T *tmpRes = new T[numElements];
    for (int i=0; i<numElements; ++i)
    {
        tmpRes[i] = left - right.m_matrixData[i];
    }

    rtMatrix2D<T> result(numRows, numCols, tmpRes);
    delete[] tmpRes;
    return result;
}
// matrix - scalar
template <class T>
rtMatrix2D<T> operator- (const rtMatrix2D<T>& left, const T& right)
{
    int numRows = left.m_nRows;
    int numCols = left.m_nCols;
    int numElements = numRows * numCols;
    T *tmpRes = new T[numElements];
    for (size_t i=0; i<numElements; ++i)
    {
        tmpRes[i] = left.m_matrixData[i] - right;
    }

    rtMatrix2D<T> result(numRows, numCols, tmpRes);
    delete[] tmpRes;
    return result;
}

//operator*
// matrix * matrix (matrix multilication)
template <class T>
rtMatrix2D<T> operator* (const rtMatrix2D<T>& lhs, const rtMatrix2D<T>& rhs)
{
    int r_numRows = rhs.m_nRows;
    int r_numCols = rhs.m_nCols;
    int l_numRows = lhs.m_nRows;
    int l_numCols = lhs.m_nCols;
    T *tmpRes = new T[lhs.m_nRows * rhs.m_nCols];

    for (size_t lhsRow=0; lhsRow<l_numRows; ++lhsRow)
    {
        for (size_t rhsCol=0; rhsCol<r_numCols; ++rhsCol)
        {
            T elementResult = 0.0;
            for (size_t lhsCol=0; lhsCol<r_numCols; ++lhsCol)
            {
                int lhsLinearIndex = (lhsRow * l_numCols) + lhsCol;
                int rhsLinearIndex = (lhsCol * r_numRows) + rhsCol;
                elementResult += (lhs.m_matrixData[lhsLinearIndex] * rhs.m_matrixData[rhsLinearIndex]);

            }
            int resultLinearIndex = (lhsRow * r_numCols) + rhsCol;
            tmpRes[resultLinearIndex] = elementResult;
        }
    }
    rtMatrix2D result(l_numRows, r_numCols, tmpRes);
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
    T *tmpRes = new T[numElements];
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
rtMatrix2D<T> operator* (const rtMatrix2D<T>& left, const T& right)
{
    int numRows = left.m_nRows;
    int numCols = left.m_nCols;
    int numElements = numRows * numCols;
    T *tmpRes = new T[numElements];
    for (size_t i=0; i<numElements; i++)
    {
        tmpRes[i] = left.m_matrixData[i] * right;
    }

    rtMatrix2D<T> result(numRows, numCols, tmpRes);
    delete[] tmpRes;
    return result;
}

// Dot product (scalar)
template <class T>
T rtMatrix2D<T>::dot(const rtMatrix2D<T>& right)
{
    int numRows = right.m_nRows;
    int numCols = right.m_nCols;
    int numElements = numRows * numCols;

    T result = 0.0;
    for (size_t i=0; i<numElements; ++i)
    {
        result += m_matrixData[i] * right.m_matrixData[i];
    }
    return result;
}

// operator==
template <class T>
bool rtMatrix2D<T>::operator== (const rtMatrix2D<T>& right)
{
    if ((this->m_nRows != right.m_nRows) && (this->m_nCols != right.m_nCols))
    {
        return false;
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
        return false;
    }

    double cumulativeSum = 0.0;
    for (size_t i=0; i<m_nElements; i++)
    {
        T element1 = matrix1.m_matrixData[i];
        T element2 = m_matrixData[i];
        cumulativeSum += (element1 - element2) * (element1 - element2);

    }
    double finalSum = sqrt(cumulativeSum / ((numRows1 * numCols1) - 1));
    if (finalSum < tolerance)
    {
        return true;
    }
    else
    {
        return false;
    }
}

/* Inverse matrix */
template <class T>
bool rtMatrix2D<T>::inverse()
{
    if (not this->isSquare())
    {
        return false;
    }
    rtMatrix2D<T> identityMatrix(m_nRows, m_nCols);
    identityMatrix.setToIdentity();
    int originalNumCols = m_nCols;
    join(identityMatrix);
    std::cout << "Joint matrix\n";

    int cRow, cCol;
    int maxCount = 100;
    int count = 0;
    bool completeFlag = false;
    while ((!completeFlag) && (count < maxCount))
    {
        for (size_t diagIndex=0; diagIndex<m_nRows; ++diagIndex)
        {
            cRow = diagIndex;
            cCol = diagIndex;

            int maxIndex = rowWithMaxElement(cCol, cRow);
            if (maxIndex != cRow)
            {
                swapRow(cRow, maxIndex);
            }
            if (m_matrixData[sub2Index(cRow, cCol)] != 1.0)
            {
                T multFactor = 1.0 / m_matrixData[sub2Index(cRow, cCol)];
                multRow(cRow, multFactor);
            }

            for (size_t rowIndex=cRow+1; rowIndex<m_nRows; ++rowIndex)
            {
                if (!closeEnought(m_matrixData[sub2Index(rowIndex, cCol)], 0.0))
                {
                    int rowOneIndex = cCol;
                    T currentElementValue = m_matrixData[sub2Index(rowIndex, cCol)];
                    T rowOneValue = m_matrixData[sub2Index(rowOneIndex, cCol)];
                    if (!closeEnought(rowOneValue, 0.0))
                    {
                        T correctionFactor = -(currentElementValue / rowOneValue);
                        multAdd(rowIndex, rowOneIndex, correctionFactor);
                    }
                }
            }
            for (size_t colIndex=cCol+1; colIndex<originalNumCols; ++colIndex)
            {
                if (!closeEnought(m_matrixData[sub2Index(cRow, colIndex)], 0.0))
                {
                    int rowOneIndex = colIndex;
                    T currentElementValue = m_matrixData[sub2Index(cRow, colIndex)];
                    T rowOneValue = m_matrixData[sub2Index(rowOneIndex, colIndex)];
                    if (!closeEnought(rowOneValue, 0.0))
                    {
                        T correctionFactor = -(currentElementValue / rowOneValue);
                        multAdd(cRow, rowOneIndex, correctionFactor);
                    }
                }
            }
        }
        rtMatrix2D<T> leftHalf;
        rtMatrix2D<T> rightHalf;
        this->separate(&leftHalf, &rightHalf, originalNumCols);
        if (leftHalf == identityMatrix)
        {
            completeFlag = true;
            m_nCols = originalNumCols;
            std::cout << originalNumCols << " Origin \n";
            m_nElements = m_nRows * m_nCols;
            delete[] m_matrixData;
            m_matrixData = new T[m_nElements];
            for (size_t i=0; i<m_nElements; ++i)
            {
                m_matrixData[i] = rightHalf.m_matrixData[i];
            }
        }
        count++;
    }
    return completeFlag;
}


#endif // LINEARALGEBRA_H
