#include <iostream>
//#include "E:\Programming projects\C++\Quasichemical model\lib\eigen-3.4.0\eigen-3.4.0\Eigen\Dense"
//#include "lib/eigen-3.4.0/eigen-3.4.0/Eigen/Dense"
#include "Eigen/Dense"
#include "Model.h"
#include "scheutjens_model.h"
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
        std::cout << "Iteration: " << i << std::endl;
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
        //std::cout << "Inverse Jacobian matrix J_inv: " << std::endl;
        //std::cout << J_inv << std::endl; // Print inverse Jacobian matrix
       
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
    std::cout <<"***********************"<< std::endl;
    std::cout <<"***********************"<< std::endl;
    std::cout << "SCHEUTJENS MODEL TEST" << std::endl;
    std::cout <<"***********************"<< std::endl;
    std::cout <<"***********************"<< std::endl; 
    std::cout << "Press any button to continue..."<< std::endl; 
    int s=0;
    std::cin >> s;

    Scheutjens_Model poly_model;
    poly_model.Init();

    std::cout<< "Scheutjens model initialized: " << poly_model.initialized << std::endl;
    int m = poly_model.M;
    std::cout << "M: " << m << std::endl;
    std::cout << "r: " << poly_model.r << std::endl;

    double* X_initial = new double[m];
    double* X_res = new double[m];
    for(int i = 0; i < m; i++) {
        X_initial[i] = log(poly_model.phi_bulk[0]/(1-poly_model.phi_bulk[0])); // Initial guess for Scheutjens model
    }

    newthon([&poly_model](Eigen::VectorXd& xs, Eigen::VectorXd& ys, int n) { poly_model.function(xs, ys, n); }, X_res, m, X_initial, 1e-8, 1000);
    
    std::cout << "Scheutjens model results: ";
    for(int i = 0; i < m; i++) {
        double phi_i = exp(X_res[i]) / (1 + exp(X_res[i])); // Convert results back to phi_i
        std::cout << phi_i << " " << std::endl;
    }

    Eigen::VectorXd p_i = poly_model.p_matrix.col(0);
    std::cout << "p_i: " << std::endl;
    std::cout << p_i<< std::endl;
}
