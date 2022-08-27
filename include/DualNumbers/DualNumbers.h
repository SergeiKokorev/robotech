#ifndef DUALNUMBERS_H
#define DUALNUMBERS_H

#include <vector>
#include <cmath>

template <class T>
class DNumRobotech
{
    public:
        // Constructors
        DNumRobotech();
        // Construct from linear 1D
        DNumRobotech(const T *inData);
        // Construct from two nums
        DNumRobotech(const T& real, const T& dual);
        // Copy constructor
        DNumRobotech(const DNumRobotech<T>& DN);
        // Construct from vector
        DNumRobotech(const std::vector<T> *inData);

        // Destructor
        ~DNumRobotech();

        // Element access
        T real() const;
        T dual() const;

        // Overloaded operators
        // Overload operator+
        template <class U> friend DNumRobotech<U> operator+ (const DNumRobotech<U>& lhs, const DNumRobotech<U>& rhs);
        template <class U> friend DNumRobotech<U> operator+ (const DNumRobotech<U>& lhs, const U& rhs);
        template <class U> friend DNumRobotech<U> operator+ (const U& lhs, const DNumRobotech<U>& rhs);
        // Overload operator-
        template <class U> friend DNumRobotech<U> operator- (const DNumRobotech<U>& lhs, const DNumRobotech<U>& rhs);
        template <class U> friend DNumRobotech<U> operator- (const DNumRobotech<U>& lhs, const U& rhs);
        template <class U> friend DNumRobotech<U> operator- (const U& lhs, const DNumRobotech<U>& rhs);
        // Overload operator*
        template <class U> friend DNumRobotech<U> operator* (const DNumRobotech<U>& lhs, const DNumRobotech<U>& rhs);
        template <class U> friend DNumRobotech<U> operator* (const DNumRobotech<U>& lhs, const U& rhs);
        template <class U> friend DNumRobotech<U> operator* (const U& lhs, const DNumRobotech<U>& rhs);
        // Overload operator/
        template <class U> friend DNumRobotech<U> operator/ (const DNumRobotech<U>& lhs, const DNumRobotech<U>& rhs);
        template <class U> friend DNumRobotech<U> operator/ (const DNumRobotech<U>& lhs, const U& rhs);
        template <class U> friend DNumRobotech<U> operator/ (const U& lhs, const DNumRobotech<U>& rhs);

        bool operator== (const DNumRobotech<T>& rhs);

        // Conjugate DN
        DNumRobotech<T> conjugate() const;

        private:
            bool closeEnought(T f1, T f2);

        private:
            T *m_dnData;
            T m_real;
            T m_dual;
};

/*
********************************
    Constructos and Destructor
*******************************
*/
// Default
template <class T>
DNumRobotech<T>::DNumRobotech() : m_real(0.0), m_dual(0.0) 
{
    m_dnData = new T[2];
    m_dnData[0] = m_real;
    m_dnData[1] = m_dual;
}
// From two nums
template <class T>
DNumRobotech<T>::DNumRobotech(const T& real, const T& dual)
{
    m_real = real;
    m_dual = dual;
    m_dnData = new T[2];
    m_dnData[0] = m_real;
    m_dnData[1] = m_dual;
}
// From 1D linear array
template <class T>
DNumRobotech<T>::DNumRobotech(const T *inData)
{
    m_real = inData[0];
    m_dual = inData[1];
    m_dnData = new T[2];
    m_dnData[0] = m_real;
    m_dnData[1] = m_dual;
}
// Copy constructor
template <class T>
DNumRobotech<T>::DNumRobotech(const DNumRobotech<T>& DN)
{
    m_real = DN.real();
    m_dual = DN.dual();
    m_dnData = new T[2];
    m_dnData[0] = m_real;
    m_dnData[1] = m_dual;
}
// From  vector
template <class T>
DNumRobotech<T>::DNumRobotech(const std::vector<T> *inData)
{
    m_real = inData->at(0);
    m_dual = inData->at(1);
    m_dnData = new T[2];
    m_dnData[0] = m_real;
    m_dnData[1] = m_dual;
}
// Destructor
template <class T>
DNumRobotech<T>::~DNumRobotech()
{
    if (m_dnData != nullptr){
        delete[] m_dnData;
    }
}

/*
********************************
        ELEMENT ACCESS
********************************
*/
template <class T>
T DNumRobotech<T>::real() const
{
    return m_dnData[0];
}
template <class T>
T DNumRobotech<T>::dual() const
{
    return m_dnData[1];
}

/*
********************************
        Utility method
********************************
*/
template <class T>
bool DNumRobotech<T>::closeEnought(T f1, T f2)
{
    return (f1 - f2) < 1e-9;
}

/*
********************************
        Conjugate
********************************
*/
template <class T>
DNumRobotech<T> DNumRobotech<T>::conjugate() const
{
    T real = m_real;
    T dual = (-1) * m_dual;
    DNumRobotech conj(real, dual);
    return conj;
}

/*
********************************
        Overloaded operators
********************************
*/
/* Overload operator+ */
// DN + DN
template <class T>
DNumRobotech<T> operator+ (const DNumRobotech<T>& lhs, const DNumRobotech<T>& rhs)
{
    T real = lhs.real() + rhs.real();
    T dual = lhs.m_dual + rhs.m_dual;
    DNumRobotech<T> result(real, dual);
    return result;
}
// DN + scalar
template <class T>
DNumRobotech<T> operator+ (const DNumRobotech<T>& lhs, const T& rhs)
{
    T real = lhs.m_real + rhs;
    T dual = lhs.m_dual;
    DNumRobotech<T> result(real, dual);
    return result;
}
// Scalar + DN
template <class T>
DNumRobotech<T> operator+ (const T& lhs, const DNumRobotech<T>& rhs)
{
    T real = rhs.m_real + lhs;
    T dual = rhs.m_dual;
    DNumRobotech<T> result(real, dual);
    return result;
}

/* Overload operator- */
// DN - DN
template <class T>
DNumRobotech<T> operator- (const DNumRobotech<T>& lhs, const DNumRobotech<T>& rhs)
{
    T real = lhs.m_real - rhs.m_real;
    T dual = lhs.m_dual - rhs.m_dual;
    DNumRobotech<T> result(real, dual);
    return result;
}
// DN - scalar
template <class T>
DNumRobotech<T> operator- (const DNumRobotech<T>& lhs, const T& rhs)
{
    T real = lhs.m_real - rhs;
    T dual = lhs.m_dual;
    DNumRobotech<T> result(real, dual);
    return result;
}
// Scalar - DN
template <class T>
DNumRobotech<T> operator- (const T& lhs, const DNumRobotech<T>& rhs)
{
    T real = lhs - rhs.m_real;
    T dual = rhs.m_dual;
    DNumRobotech<T> result(real, dual);
    return result;
}

/* Overload operator* */
// DN * DN
template <class T>
DNumRobotech<T> operator* (const DNumRobotech<T>& lhs, const DNumRobotech<T>& rhs)
{
    T real = lhs.m_real * rhs.m_real;
    T dual = lhs.m_real * rhs.m_dual + lhs.m_dual * rhs.m_real;
    DNumRobotech<T> result(real, dual);
    return result;
}
// DN * scalar
template <class T>
DNumRobotech<T> operator* (const DNumRobotech<T>& lhs, const T& rhs)
{
    T real = lhs.m_real * rhs;
    T dual = lhs.m_dual * rhs;
    DNumRobotech<T> result(real, dual);
    return result;
}
// Scalar * DN
template <class T>
DNumRobotech<T> operator* (const T& lhs, const DNumRobotech<T>& rhs)
{
    T real = rhs.m_real * lhs;
    T dual = rhs.m_dual * lhs;
    DNumRobotech<T> result(real, dual);
    return result;
}

/* Overload operator/ */
// DN / DN
template <class T>
DNumRobotech<T> operator/ (const DNumRobotech<T>& lhs, const DNumRobotech<T>& rhs)
{
    T real = lhs.m_real / rhs.m_real;
    T dual = (lhs.m_dual * rhs.m_real - lhs.m_real * rhs.m_dual) / std::pow(rhs.m_real, 2);
    DNumRobotech<T> result(real, dual);
    return result;
}
// DN / scalar
template <class T>
DNumRobotech<T> operator/ (const DNumRobotech<T>& lhs, const T& rhs)
{
    T real = lhs.m_real / rhs;
    T dual = lhs.m_dual / rhs;
    DNumRobotech<T> result(real, dual);
    return result;
}
// Scalar / DN
template <class T>
DNumRobotech<T> operator/ (const T& lhs, const DNumRobotech<T>& rhs)
{
    T real = lhs / rhs.m_real;
    T dual = lhs * rhs.m_dual / std::pow(rhs.m_real, 2);
    DNumRobotech<T> result(real, dual);
    return result;
}

/* Overload operator== */
template <class T>
bool DNumRobotech<T>::operator== (const DNumRobotech<T>& rhs)
{
    bool flag = true;
    for (size_t i=0; i<2; ++i)
    {
        if (!closeEnought(this->m_dnData[i], rhs.m_dnData[i]))
        {
            flag = false;
        }
    }
    return flag;
}

#endif // DUALNUMBERS_H
