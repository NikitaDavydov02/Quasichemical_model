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
void newthon(std::function<void(Eigen::VectorXd&,Eigen::VectorXd&,int)> func, double* res, int n, double  break_coeff=1.0, double* initial_guess = NULL, double tolerance = 1e-6, int max_iterations = 1000)
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
    
    //<COMMENT>
    //xs << 0.40815,	0.17409	,0.069353	,0.039993,	0.02751,	0.021089,	0.017388,	0.015084,	0.013569,	0.012513,	0.011798,	0.011267	,0.010878,	0.0105,	0.010374,	0.010214	,0.010096,	0.01001,	0.0099499,	0.0099087;
    //</COMMENT>

    for(int i=0;i<max_iterations;i++)
    {
        std::cout << "Iteration: " << i << std::endl;
        func(xs, fs, n);
        //Check convergence
        
        double delta = fs.dot(fs); // Dot product to calculate the sum of squares

        std::cout<< std::endl;
        std::cout << "Delta f: " << fs << std::endl;
        std::cout<< std::endl;

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
        xs-= step_vec*break_coeff; // Update xs with step vector

        std::cout << "Updated xs: " << xs.transpose() << std::endl; // Print updated xs
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
    
    Scheutjens_Model poly_model;
    poly_model.Init();

    std::cout<< "Scheutjens model initialized: " << poly_model.initialized << std::endl;
    int m = poly_model.M;
    std::cout << "M: " << m << std::endl;
    std::cout << "r: " << poly_model.r << std::endl;

    std::cout << "Press any button to continue..."<< std::endl; 
    int s=0;
    std::cin >> s;

    double* X_initial = new double[m];
    double* X_res = new double[m];
    for(int i = 0; i < m; i++) {
        X_initial[i] = log(poly_model.phi_bulk[0]/(1-poly_model.phi_bulk[0])); // Initial guess for Scheutjens model
    }

    double break_coeff=0.1;
    newthon([&poly_model](Eigen::VectorXd& xs, Eigen::VectorXd& ys, int n) { poly_model.function(xs, ys, n); }, X_res, m, break_coeff, X_initial, 1e-10, 1000);
    
    std::cout << "Scheutjens model results: " << std::endl;
    for(int i = 0; i < m; i++) {
        double phi_i = exp(X_res[i]) / (1 + exp(X_res[i])); // Convert results back to phi_i
        std::cout << phi_i << " " << std::endl;
    }

    std::cout << std::endl;
    Eigen::VectorXd p_i = poly_model.p_matrix.col(0);
    std::cout << "p_i: " << std::endl;
    std::cout << p_i<< std::endl;

    std::cout << "Convergence check: " << std::endl;
    std::cout << "Press any button to continue..." << std::endl;
    std::cin >> s;

    Eigen::VectorXd X_res_vec(m);
    Eigen::VectorXd ln_conv(m);
    for(int i=0;i<m;i++)
    {
        X_res_vec(i)=X_res[i];
    }  

    poly_model.function(X_res_vec,ln_conv,m);
    std::cout << "ln_conv: " << std::endl;
    std::cout << ln_conv << std::endl;
    std::cout << "p_matrix: " << std::endl;
    std::cout << poly_model.p_matrix << std::endl;
    std::cout << "phi_i_dashed: " << std::endl;
    std::cout << poly_model.phi_i_dashed << std::endl;
    std::cout << "phi_i: " << std::endl;
    std::cout << poly_model.phi_i << std::endl;
}
