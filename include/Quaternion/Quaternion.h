#ifndef QUTERNION_H
#define QUTERNION_H

#include <vector>
#include <cmath>
#include <iostream>

#include "Vector.h"


template <class T>
class QRobotech
{
    public:
        // Constructors
        QRobotech();
        // Constructors from four nums
        // q1, q2, q3 will be transformed to class Vector
        QRobotech(const T& q0, const T& q1, const T& q2, const T& q3);
        // Constructor from a num and std::vector
        // std::vector will be translated to class Vector
        QRobotech(const T& q0, const std::vector<T>& q);
        // Constructor from a num and class Vector
        QRobotech(const T& q0, const VectorRobotech<T>& q);
        // Copy constructir
        QRobotech(const QRobotech<T>& Q);
        
        // Destructor
        ~QRobotech();

        // Element access
        // Getter methods
        T Re();
        VectorRobotech<T> Im();
        // Argument of a quternion
        T arg();
        // Versor of a quternion
        QRobotech<T> vers();

        // Setter methods
        bool setRe();
        bool setIm();

        // Overloaded operators
        template <class U> friend QRobotech<U> operator+ (const QRobotech<U>& lhs, const QRobotech<U>& rhs);
        template <class U> friend QRobotech<U> operator- (const QRobotech<U>& lhs, const QRobotech<U>& rhs);
        // Simple multiplications
        template <class U> friend QRobotech<U> operator^ (const QRobotech<U>& lhs, const QRobotech<U>& rhs);
        template <class U> friend QRobotech<U> operator^ (const QRobotech<U>& lhs, const U& rhs);
        template <class U> friend QRobotech<U> operator^ (const U& lhs, const QRobotech<U>& rhs);
        // Quternion multiplication
        template <class U> friend QRobotech<U> operator* (const QRobotech<U>& lhs, const QRobotech<U>& rhs);

        // Scalar product
        QRobotech<T> sDot();
        // Vector product
        QRobotech<T> vDot();

        // Quaternion manipulations
        // Quaternion inversion
        QRobotech<T> inverse();
        // Conjugated quternion
        QRobotech<T> conj();
        // Matrix representation
        // rtMatrix2D<T> matrix()
        // Module
        T mod();
        // Module of a vector part
        T vModule();
        // Module of a scalar part
        T sModule();
        // Normalized quternion
        QRobotech<T> norm();



    private:
        T m_q0;
        VectorRobotech<T> *m_vector;
        T *m_qData;
};

/* *******************************
    Constructors & Desctructor
*********************************/
// Default
template <class T>
QRobotech<T>::Qrobotech()
{
    m_q0 = 1.0;
    m_vector = new VectorRobotech(0.0, 0.0, 0.0);
    m_qData = new T[2];
    m_qData[0] = m_q0;
    m_qData[1] = m_vector;
}

// Destructor
template <class T>
QRobotech<T>::~QRobotech()
{
    if(m_vector != nullprt)
    {
        delete[] m_vector;
        std::cout << "m_vector deleted successfully\n";
    }
    if(m_qData != nullptr)
    {
        delete[] m_qData;
        std::cout << "m_qData deleted successfully\n";
    }
}

#endif // QUTERNION_H
