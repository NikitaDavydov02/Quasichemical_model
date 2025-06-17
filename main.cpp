#include <iostream>
//#include "E:\Programming projects\C++\Quasichemical model\lib\eigen-3.4.0\eigen-3.4.0\Eigen\Dense"
//#include "lib/eigen-3.4.0/eigen-3.4.0/Eigen/Dense"
#include "Eigen/Dense"
#include "Model.h"
#include <functional>

using namespace Eigen;

void newthon(std::function<double(double)> func, double* res, double initial_guess = 0.0, double tolerance = 1e-6, int max_iterations = 1000)
{
    // Placeholder for Newton's method implementation
    // This function will be implemented in the future
    std::cout << "Newton's method placeholder called." << std::endl;
    double x = initial_guess; 
    for(int i=0;i<max_iterations;i++)
    {
        double fx = func(x);
        if (std::abs(fx) < tolerance) {
            std::cout << "Converged to: " << x << std::endl;
            if (res) {
                *res = x; // Store the result if res is not null
            }
            return;
        }
        // Update x based on the derivative (not implemented here)
        double dx= 0.01; // Placeholder for derivative step
        double f_df = func(x + dx) - fx; // Placeholder for derivative calculation
        double derivative = f_df / dx; // Placeholder for derivative
        double step = fx / derivative; // Newton's step
        x -= step; // Update x
    }
    return;
}
/**
 * @brief main function.
 * @param none
 * @return 0 if succeed
 * @throws none
 */
int main()
{
    std::cout << "Hello world!" << std::endl;
    Matrix<double,3,3> matrixA;
    matrixA.setZero();
    matrixA(0, 0) = 1.0;
    matrixA(1, 1) = 2.0;   
    matrixA(2, 2) = 2.0;
    std::cout << matrixA << std::endl;
    Matrix3d m = matrixA.inverse();
    std::cout << m << std::endl;
    ///
    Model model;
    model.initialized = false;
    model.Init();
    std::cout << model.initialized << std::endl;

    double x;
    newthon([&model](double x) { return model.function(x); },&x, 0.0, 1e-6, 1000);
    std::cout << "Result from Newton's method: " << x << std::endl;
}
