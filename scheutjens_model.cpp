#include "scheutjens_model.h"
#include <iostream>


void Scheutjens_Model::Init() {
    // Initialization logic here
    z=10;       
    r=10;       // Chain length
    M=3*r;      // Number of layers
    chi = 1.2;  // Flory parameter
    chi_s=0.0;  // Adsorption parameter
    phi_bulk[0]=0.01;
    phi_bulk[1]=1-phi_bulk[0];

    lamda_matrix.setZero(M,M);
    p_matrix.setZero(M,r);
    w_matrix.setZero(M,M);

    phi_i.setConstant(M, phi_bulk[0]); // Initialize phi_i with the bulk polymer fraction
    phi_solv_i.setConstant(M, phi_bulk[1]); // Initialize phi_solv_i with the bulk solvent fraction
    avr_phi_solv_i.setZero(M);
    avr_phi_i.setZero(M);
    phi_i_dashed.setZero(M);
    
    double lamda_0 = 4.0/6.0;
    double lamda_1 = 1.0/6.0;

    // Initialize the lambda matrix with some values, for example:
    for(int i=0; i<M; i++) {
        for(int j=0; j<M; j++) {
            if(i == j)
                lamda_matrix(i, j) = lamda_0; // Diagonal elements
            else if(j==i+1 || j==i-1)
                lamda_matrix(i, j) = lamda_1; //Neighbour elements
            else
                lamda_matrix(i, j) = 0.0; // Other elements 
        }
    }
    std::cout << "Lambda matrix initialized:" << std::endl;
    Print_Matrix(lamda_matrix);

    lnP_bulk =chi*(phi_bulk[0]-phi_bulk[1])+log(phi_bulk[1]);//[TO DO] check if this is correct

    initialized = true; // Set initialized to true after initialization
}
void Scheutjens_Model::function( Eigen::VectorXd& xs,  Eigen::VectorXd& ys, int n)
{
    //xs = ln(phi_i/(1-phi_i))
    //ys = ln(phi_i_dashed/phi_i)

    for(int i=0; i<M; i++)
    {
        phi_i(i) = exp(xs(i)) / (1 + exp(xs(i))); // Convert xs to phi_i
        phi_solv_i(i) = 1 - phi_i(i); // Assuming solvent fraction is 1 - polymer fraction
    }

    Print_Vector(phi_i);
    Print_Vector(phi_solv_i);

    Calculate_avr_phi_i();

    Print_Vector(avr_phi_i);
    Print_Vector(avr_phi_solv_i);


    //<calculate lnPi>
    for(int i=0; i<M; i++)
    {
       // double lnP_i = chi*(avr_phi_i[i]-avr_phi_solv_i[i])+log(avr_phi_solv_i[i]);//[TO DO] check if this is correct
        double lnP_i = chi*(avr_phi_i[i]-avr_phi_solv_i[i])+log(avr_phi_solv_i[i]);

        if(i==0)
            lnP_i += chi_s;
        double P_i = exp(lnP_i);
        double P_bulk = exp(lnP_bulk);
        double p_i = P_i / P_bulk; // Calculate p_i
        p_matrix(i, 0) = p_i; // Store in p_matrix
        for(int j=0; j<M; j++)
        {
            w_matrix(i, j) = lamda_matrix(i, j) * p_i; // Calculate w_matrix
        }
    }
    Print_Matrix(p_matrix);
    //</calculate lnPi>

    //<calculate p_matrix>
    Eigen::VectorXd p_col_prev = p_matrix.col(0);
    Eigen::VectorXd p_col_new;
    for(int i=1; i<r; i++)
    {
        p_col_new = w_matrix * p_col_prev; // Calculate new column of p_matrix
        p_matrix.col(i) = p_col_new; // Store in p_matrix
        p_col_prev = p_col_new; // Update for next iteration
    }
    Print_Matrix(p_matrix);
    //</calcualte p_matrix>

    //<calculate phi_i_dashed>
    for(int i=0;i<M;i++)
    {
        double sum=0;
        for(int s=0;s<r;s++)
        {
            double p_1 = p_matrix(i, s);
            double p_2 = p_matrix(i,r-1-s);
            sum += p_matrix(i, s) * p_matrix(i,r-1-s);
            // sum += p_matrix(i, s) * p_matrix(i,r-s+1);
        }
        phi_i_dashed(i)=phi_bulk[0]*sum/(r*p_matrix(i,0)); // Calculate phi_i_dashed
    }
    Print_Vector(phi_i_dashed);
    //</calculate phi_i_dashed>

    //<calculate ys>
    for(int i=0; i<M; i++)
    {
        ys(i) = log(phi_i_dashed(i) / phi_i(i)); // Calculate ys
    }
    Print_Vector(ys);
    //</calculate ys>
}
void Scheutjens_Model::Calculate_avr_phi_i()
{
    avr_phi_i = lamda_matrix * phi_i;
    avr_phi_solv_i = lamda_matrix * phi_solv_i;
}
void Scheutjens_Model::Print_Vector(const Eigen::VectorXd& vec) {
    //return; // Commented out to avoid excessive output
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << vec << std::endl;
}
void Scheutjens_Model::Print_Matrix(const Eigen::MatrixXd& mat) {
   // return; // Commented out to avoid excessive output
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << mat << std::endl;
}