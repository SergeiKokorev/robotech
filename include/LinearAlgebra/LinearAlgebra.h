#ifndef LINEARALGEBRA_H
#define LINEARGEBRA_H

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
        rtMatrix2D(const rtMatrix2D<T>& inMatrix);
        rtMatrix2D(int nRows, int nCols, const std::vector<T> *inData);

        // Destructor
        ~rtMatrix2D();

        // Matrix cross product
        template <class U> rtMatrix2D<U> multiply(const rtMatrix2D<U>& left, const rtMatrix2D<U>& right);

        // Iverse matrix

        // Manipulate methods
        // Transpose matrix
        bool transpose();
        bool inverse();

        // Configure methods
        bool resize(int numRows, int num Cols);
        void setToIdentity();

        // Element access
        T getElement(int row, int col);
        bool setElement(int row, int col, T value);
        int getNumRows();
        int getNumCols();

        // Overload 
        bool operator== (const rtMatrix2D<T>& right);
        bool compare(const rtMatrix2D<T>& matrix1, double tolerance);

        // Overload +, - and * operators
        template <class U> friend rtMatrix2D<U> operator+ (const rtMatrix2D<U>& left, const rtMatrix2D<U>&right);
        template <class U> friend rtMatrix2D<U> operator+ (const U& left, const rtMatrix2D<U>& right);
        template <class U> friend rtMatrix2D<U> operator+ (const rtMatrix2D<U>& left, const U right);

        template <class U> friend rtMatrix2D<U> operator- (const rtMatrix2D<U>& left, const rtMatrix2D<U>& right);
        template <class U> friend rtMatrix2D<U> operator- (const U& left, const rtMatrix2D<U>& right);
        template <class U> friend rtMatrix2D<U> operator- (const rtMatrix2D<U>& left, const U& right);

        template <class U> friend rtMatrix2D<U> operator* (const rtMatrix2D<U>& left, const rtMatrix2D<U>& right);
        template <class U> friend rtMatrix2D<U> operator* (const U& left, const rtMatrix2D<U>& right);
        template <class U> friend rtMatrix2D<U> operator* (const rtMatrix2D<U>& left, const U& right);

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
        int m_nRows, m_nCols, m_nElements

};

#endif LINEARALGEBRA_H
