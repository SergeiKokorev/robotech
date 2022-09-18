#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>
#include <vector>
#include <iostream>

#include "LinearAlgebra.h"

template <class T>
class VectorRobotech
{
    public:
        // Constructors
        VectorRobotech();
        // From 1D linear array
        VectorRobotech(const T *inData);
        // Constructor from vector
        VectorRobotech(const std::vector<T> *inData);
        // Copy constructor
        VectorRobotech(const VectorRobotech<T>& V);
        // Constructor from numbers
        VectorRobotech(const T& x, const T& y, const T& z);
        // Destructor
        ~VectorRobotech();

        // Elements access method
        T getElement(size_t i) const;
        T getX() const;
        T getY() const;
        T getZ() const;

        // Element setter method
        VectorRobotech<T> setElement(size_t i, const T& e);
        VectorRobotech<T> setX(const T& x);
        VectorRobotech<T> setY(const T& y);
        VectorRobotech<T> setZ(const T& z);

        // Overload operators mehtods
        // operator+
        template <class U> friend VectorRobotech<U> operator+ (const VectorRobotech<U>& lhs, const VectorRobotech<U>& rhs);
        // operator-
        template <class U> friend VectorRobotech<U> operator- (const VectorRobotech<U>& lhs, const VectorRobotech<U>& rhs);
        // operator* vector (cross) product
        template <class U> friend VectorRobotech<U> operator* (const VectorRobotech<U>& lhs, const VectorRobotech<U>& rhs);
        // Scalar (dot) product
        template <class U> friend U operator^ (const VectorRobotech<U>& lhs, const VectorRobotech<U>& rhs);

        // Other vector methods
        // Outer product
        rtMatrix2D<T> outer(const VectorRobotech<T>& rhs) const;
        rtMatrix2D<T> wedge(const VectorRobotech<T>& rhs) const;


        // Vector magnitude
        T magnitude() const;
        // Unit (normalized) vector
        bool norm();

    protected:
        T m_x;
        T m_y;
        T m_z;
        T *m_vData;
};

/* ******************************
        Constructors
****************************** */
// Default constructor
template <class T>
VectorRobotech<T>::VectorRobotech() : m_x(1.0), m_y(1.0), m_z(1.0)
{
    m_vData = new T[3];
    m_vData[0] = m_x;
    m_vData[1] = m_y;
    m_vData[2] = m_z;
}
// From 1D linear array
template <class T>
VectorRobotech<T>::VectorRobotech(const T *inData)
{
    m_vData = new T[3];
    for (size_t i=0; i<3; ++i)
    {
        m_vData[i] = inData[i];
    }
    m_x = m_vData[0];
    m_y = m_vData[1];
    m_z = m_vData[2];
}
// From vector
template <class T>
VectorRobotech<T>::VectorRobotech(const std::vector<T> *inData)
{
    m_vData = new T[3];
    for (size_t i=0; i<3; ++i)
    {
        m_vData[i] = inData->at(i);
    }
    m_x = m_vData[0];
    m_y = m_vData[1];
    m_z = m_vData[2];
}
// From nums
template <class T>
VectorRobotech<T>::VectorRobotech(const T& x, const T& y, const T& z)
{
    m_x = x;
    m_y = y;
    m_z = z;
    m_vData = new T[3];
    m_vData[0] = m_x;
    m_vData[1] = m_y;
    m_vData[2] = m_z;
}
// Copy constructor
template <class T>
VectorRobotech<T>::VectorRobotech(const VectorRobotech<T>& V)
{
    m_x = V.m_x;
    m_y = V.m_y;
    m_z = V.m_z;
    m_vData = new T[3];
    for (size_t i=0; i<3; ++i)
    {
        m_vData[i] = V.m_vData[i];
    }
}

// Destructor
template <class T>
VectorRobotech<T>::~VectorRobotech()
{
    if (m_vData != nullptr)
    {
        delete[] m_vData;
    }
}

/* ****************************
        Element access
**************************** */
template <class T>
T VectorRobotech<T>::getElement(size_t i) const
{
    return m_vData[i];
}
template <class T>
T VectorRobotech<T>::getX() const
{
    return m_x;
}
template <class T>
T VectorRobotech<T>::getY() const
{
    return m_y;
}
template <class T>
T VectorRobotech<T>::getZ() const
{
    return m_z;
}

/* ****************************
        Setter methods
**************************** */
template <class T>
VectorRobotech<T> VectorRobotech<T>::setElement(size_t i, const T& value)
{
    m_vData[i] = value;
    return *this;
}
template <class T>
VectorRobotech<T> VectorRobotech<T>::setX(const T& x)
{
    m_x = x;
    m_vData[0] = m_x;
    return *this;
}
template <class T>
VectorRobotech<T> VectorRobotech<T>::setY(const T& y)
{
    m_y = y;
    m_vData[1] = m_y;
    return *this;
}
template <class T>
VectorRobotech<T> VectorRobotech<T>::setZ(const T& z)
{
    m_z = z;
    m_vData[2] = m_z;
    return *this;
}

/* ******************************
        Overload operators
****************************** */
// Overload operator+
template <class T>
VectorRobotech<T> operator+ (const VectorRobotech<T>& lhs, const VectorRobotech<T>& rhs)
{
    VectorRobotech<T> result(
        lhs.m_x + rhs.m_x,
        lhs.m_y + rhs.m_y,
        lhs.m_z + rhs.m_z
    );
    return result;
}
// Overload operator-
template <class T>
VectorRobotech<T> operator- (const VectorRobotech<T>& lhs, const VectorRobotech<T>& rhs)
{
    VectorRobotech<T> result(
        lhs.m_x - rhs.m_x,
        lhs.m_y - rhs.m_y,
        lhs.m_z - rhs.m_z
    );
    return result;
}
// Overload operator* cross (vector) product
template <class T>
VectorRobotech<T> operator* (const VectorRobotech<T>& lhs, const VectorRobotech<T>& rhs)
{
    T x = lhs.m_y * rhs.m_z - lhs.m_z * rhs.m_y;
    T y = lhs.m_x * rhs.m_z - lhs.m_z * rhs.m_x;
    T z = lhs.m_x * rhs.m_y - lhs.m_y * rhs.m_x;
    VectorRobotech<T> result(x, y, z);
    return result;
}
// Overload operator^ dot (scalar) product
template <class T>
T operator^ (const VectorRobotech<T>& lhs, const VectorRobotech<T>& rhs)
{
    return lhs.m_x * rhs.m_x + lhs.m_y * rhs.m_y + lhs.m_z * rhs.m_z;
}

// Other vector methods
// Outer product. Return nxn matrix
template <class T>
rtMatrix2D<T> VectorRobotech<T>::outer(const VectorRobotech<T>& rhs) const
{
    int nRows = 3;
    int nCols = 3;
    rtMatrix2D<T> result(nRows, nCols);

    for (size_t row=0; row<3; ++row)
    {
        T elementRes = 0.0;        
        for (size_t col=0; col<3; ++col)
        {
            T value = this->m_vData[row] * (-1) * rhs.m_vData[col];
            result.setElement(row, col, value);
        }
    }
    return result;
}

// Wedge product as v1.outer(v2) - v2.outer(v1)
template <class T>
rtMatrix2D<T> VectorRobotech<T>::wedge(const VectorRobotech<T>& rhs) const
{
    return this->outer(rhs) - rhs.outer(*this);
}


/* ******************************
        Magnitude
****************************** */
template <class T>
T VectorRobotech<T>::magnitude() const
{
    return std::sqrt(std::pow(m_x, 2) + std::pow(m_y, 2) + std::pow(m_z, 2));
}

/* ******************************
        Normalized vector
****************************** */
template <class T>
bool VectorRobotech<T>::norm()
{
    T magnitude = this->magnitude();
    if(magnitude == 0){
        return false;
    } else {
    m_x = m_x / magnitude;
    m_y = m_y / magnitude;
    m_z = m_z / magnitude;        
    m_vData[0] = m_x;
    m_vData[1] = m_y;
    m_vData[2] = m_z;
    }
    return true;
}

#endif // VECTOR_H
