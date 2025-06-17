#include <iostream>
//#include "E:\Programming projects\C++\Quasichemical model\lib\eigen-3.4.0\eigen-3.4.0\Eigen\Dense"
//#include "lib/eigen-3.4.0/eigen-3.4.0/Eigen/Dense"
#include "Eigen/Dense"
#include "Model.h"
#include <functional>

using namespace Eigen;

void newthon_single(std::function<double(double)> func, double* res, double initial_guess = 0.0, double tolerance = 1e-6, int max_iterations = 1000)
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
void newthon(std::function<void(Eigen::VectorXd&,Eigen::VectorXd&,int)> func, double* res, int n, double* initial_guess = NULL, double tolerance = 1e-6, int max_iterations = 1000)
{
    // Placeholder for Newton's method implementation
    // This function will be implemented in the future
    std::cout << "Newton's method placeholder called." << std::endl;
    if(initial_guess == NULL) {
        std::cerr << "Error: initial_guess is NULL." << std::endl;
        return; // Handle the case where initial_guess is NULL
    }
    VectorXd xs(n);
    VectorXd fs(n);
    VectorXd f_dfs(n);
    MatrixXd J(n, n);
    for(int i=0;i<n;i++)
    {
        xs(i) = initial_guess[i]; // Initialize xs with initial_guess
    }

    for(int i=0;i<max_iterations;i++)
    {
        func(xs, fs, n);
        //Check convergence
        double delta = fs.dot(fs); // Dot product to calculate the sum of squares
       
        if (std::abs(delta) < tolerance) {
            std::cout << "Converged" << std::endl;
            for(int j=0;j<n;j++) {
                if (res) {
                    res[j] = xs(j) ; // Store the result if res is not null
                }
            }
            return;
        }
        // Update x based on the derivative (not implemented here)
        double dx= 0.01; // Placeholder for derivative step
        //Calculate Jacobian matrix with Eigen library
        
        for(int j=0;j<n;j++) {
            xs(j)  += dx; // Increment xs for Jacobian calculation
            func(xs, f_dfs, n); // Calculate function values at incremented xs
            for(int k=0;k<n;k++) {
                J(k, j) = (f_dfs(k) - fs(k)) / dx; // Calculate Jacobian column
            }
            xs(j) -= dx; // Reset xs after increment
        }
        //std::cout << "Jacobian matrix J: " << std::endl;
        //std::cout << J << std::endl; // Print Jacobian matrix
        MatrixXd J_inv = J.inverse(); // Inverse of Jacobian matrix
        std::cout << "Inverse Jacobian matrix J_inv: " << std::endl;
        std::cout << J_inv << std::endl; // Print inverse Jacobian matrix
       
        VectorXd step_vec = J_inv * fs; // Calculate step vector
        xs-= step_vec; // Update xs with step vector
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
    newthon_single([&model](double x) { return model.function(x); },&x, 0.0, 1e-6, 1000);
    std::cout << "Result from single Newton's method: " << x << std::endl;

    double* xss = new double[2];
    double* initial_guess = new double[2];
    int n = 3; // Number of variables
    for(int i = 0; i < n; i++) {
        initial_guess[i] = 0.0001; // Initial guess for multi-variable Newton's method
    }
    newthon([&model](Eigen::VectorXd& xs, Eigen::VectorXd& ys, int n) { model.function(xs, ys, n); }, xss, n, initial_guess, 1e-6, 1000);
    
    std::cout << "Result from multi-variable Newton's method: ";
    for(int i = 0; i < n; i++) {
        std::cout << xss[i] << " ";
    }
}
