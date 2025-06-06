#include <iostream>
//#include "E:\Programming projects\C++\Quasichemical model\lib\eigen-3.4.0\eigen-3.4.0\Eigen\Dense"
//#include "lib/eigen-3.4.0/eigen-3.4.0/Eigen/Dense"
#include "Eigen/Dense"


using namespace Eigen;

int main()
{
    std::cout << "Hello world!" << std::endl;
    Matrix<float,3,3> matrixA;
    matrixA.setZero();
    std::cout << matrixA << std::endl;
}