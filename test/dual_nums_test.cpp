#include <iostream>
#include <vector>
#include "DualNumber.h"

template <class T>
void printDNums(const DNumRobotech<T> dnum)
{
    std::cout << "Real part:\t" << dnum.real() << "\tDual part:\t" << dnum.dual() << "\n";
}

int main() {

    std::cout << "Constructors testing\n";

    DNumRobotech<double> dNum;
    std::cout << "Default constructor:\t" << "Real: " << dNum.real() << "\tDual: " << dNum.dual() << "\n";

    std::vector<double> v{1.5, 2.5};
    DNumRobotech vecDNum(v);
    std::cout << "Construct from vector:\tReal: " << vecDNum.real() << "\tDual: " << vecDNum.dual() << "\n";

    double linear[2]{4.5, 3.0};
    DNumRobotech linearDNum(linear);
    std::cout << "Construct from 1D linear:\tReal: " << linearDNum.real() << "\tDual: " << linearDNum.dual() << "\n";

    double f1{2.3};
    double f2{3.4};
    DNumRobotech twoDNum(f1,f2);
    std::cout << "Construct from two nums:\tReal: " << twoDNum.real() << "\tDual: " << twoDNum.dual() << "\n";

    DNumRobotech copyDNum(twoDNum);
    std::cout << "Copy constructor:\tReal: " << copyDNum.real() << "\tDual: " << copyDNum.dual() << "\n";

    std::cout << "\n#############################################\n";
    double scalar = 1.5;
    std::cout << "Test overloaded operators" << "\n";

    std::cout <<"Additions:\n";
    std::cout << "first:\t";
    printDNums(vecDNum);
    std::cout << "second:\t";
    printDNums(twoDNum);
    printDNums(vecDNum + twoDNum);
    std::cout << "Scalar additions\n" << "Scalar:\t" << scalar << "\n";
    std::cout << "Dual number:\t";
    printDNums(vecDNum);
    std::cout << "Right hand scalar" << "\n";
    printDNums(scalar + vecDNum);
    std::cout << "Left hand scalar:\t" << "\n";
    printDNums(vecDNum + scalar);

    std::cout <<"Substructions:\n";
    std::cout << "first:\t";
    printDNums(vecDNum);
    std::cout << "second:\t";
    printDNums(twoDNum);
    printDNums(vecDNum - twoDNum);
    std::cout << "Scalar substruct\n" << "Scalar:\t" << scalar << "\n";
    std::cout << "Dual number:\t";
    printDNums(twoDNum);
    std::cout << "Left hand scalar" << "\n";
    printDNums(scalar - twoDNum);
    std::cout << "Right hand scalar:\t" << "\n";
    printDNums(twoDNum - scalar);

    std::cout <<"Product:\n";
    std::cout << "first:\t";
    printDNums(vecDNum);
    std::cout << "second:\t";
    printDNums(twoDNum);
    printDNums(vecDNum * twoDNum);
    std::cout << "Scalar product\n" << "Scalar:\t" << scalar << "\n";
    std::cout << "Dual number:\t";
    printDNums(twoDNum);
    std::cout << "Left hand scalar" << "\n";
    printDNums(scalar * twoDNum);
    std::cout << "Right hand scalar:\t" << "\n";
    printDNums(twoDNum * scalar);

    std::cout <<"Divisions:\n";
    std::cout << "first:\t";
    printDNums(vecDNum);
    std::cout << "second:\t";
    printDNums(twoDNum);
    printDNums(vecDNum / twoDNum);
    std::cout << "Scalar devisions\n" << "Scalar:\t" << scalar << "\n";
    std::cout << "Dual number:\t";
    printDNums(twoDNum);
    std::cout << "Left hand scalar" << "\n";
    printDNums(scalar / twoDNum);
    std::cout << "Right hand scalar:\t" << "\n";
    printDNums(twoDNum / scalar);

    std::cout << "\n#######################################\n";
    std::cout << "Conjugate testing\n";
    std::cout << "Dual number:\t";
    printDNums(twoDNum);
    std::cout << "Conjugated:\t";
    printDNums(twoDNum.conjugate());

    std::cout << "\nCompare operator testing\n";
    std::cout << "Dual numbers:\t";
    printDNums(twoDNum);
    printDNums(vecDNum);
    std::cout << "twoDNum == vecDNum:\t" << (twoDNum == vecDNum) << "\n";
    std::cout << "vecDNum == twoDNum:\t" << (twoDNum == vecDNum) << "\n";
    std::cout << "twoDNum == copyDNum:\t" << (twoDNum == copyDNum) << "\n";
    std::cout << "copyDNum == twoDNum:\t" << (twoDNum == copyDNum) << "\n";

    return 0;
}
