#include "scheutjens_model.h"
#include <iostream>


void Scheutjens_Model::Init() {
    // Initialization logic here
      
    r=1000;       // Chain length
    M=20;      // Number of layers
    chi = 0.0;  // Flory parameter
    chi_s= 0.0;  // Adsorption parameter
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
    
    //lamda_0 = 4.0/6.0;
    //lamda_1 = 1.0/6.0;

    lamda_0=0.5;
    lamda_1=0.25;

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

    lnP_bulk =chi*(phi_bulk[0]-phi_bulk[1])+log(phi_bulk[0]);//[TO DO] check if this is correct

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


    Calculate_ln_p_i();
    Print_Matrix(p_matrix);

    Calculate_p_matrix();

    Calculate_phi_i_dashed();

    //<calculate ys>
    for(int i=0; i<M; i++)
    {
        ys(i) = log(phi_i_dashed(i) / phi_i(i)); // Calculate ys
    }
    Print_Vector(ys);
    //</calculate ys>
}
void Scheutjens_Model::Calculate_ln_p_i()
{
    for(int i=0; i<M; i++)
    {
       // double lnP_i = chi*(avr_phi_i[i]-avr_phi_solv_i[i])+log(avr_phi_solv_i[i]);//[TO DO] check if this is correct
        double lnP_i = chi*(avr_phi_i[i]-avr_phi_solv_i[i])+log(avr_phi_i[i]);

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
}
void Scheutjens_Model::Calculate_phi_i_dashed()
{
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
}
void Scheutjens_Model::Calculate_p_matrix()
{
    /*<OLD>Eigen::VectorXd p_col_prev = p_matrix.col(0);
    Eigen::VectorXd p_col_new;
    for(int i=1; i<r; i++)
    {
        p_col_new = w_matrix * p_col_prev; // Calculate new column of p_matrix
        p_matrix.col(i) = p_col_new; // Store in p_matrix
        p_col_prev = p_col_new; // Update for next iteration
    }
    Print_Matrix(p_matrix);
    </OLD>*/

    //New implementation
    Eigen::VectorXd p_col_prev = p_matrix.col(0);
    Eigen::VectorXd p_col_new;
    int last_ascent_column_index = 0;
    if(r%2==0)
        last_ascent_column_index = r/2; // Last column index for ascending part
    else
        last_ascent_column_index = (r+1)/2; // Last column index for ascending part

    for(int s=1;s<r;s++)
    {
        if(s<=last_ascent_column_index)
        {
            int new_col_size = M+s; // Size of the new column vector
            p_col_new.setZero(new_col_size);
            for(int i=0;i<new_col_size;i++)
            {
                if(i==0)
                    p_col_new(i)=lamda_0*p_col_prev(i)+lamda_1*p_col_prev(i+1); // First layer
                else if(i==new_col_size-1)
                    p_col_new(i)=lamda_1*p_col_prev(i-1)+lamda_0*1 + lamda_1*1;
                else if(i==new_col_size-2)
                    p_col_new(i)=lamda_1*p_col_prev(i-1)+lamda_0*p_col_prev(i) + lamda_1*1;
                else
                    p_col_new(i)=lamda_1*p_col_prev(i-1)+lamda_0*p_col_prev(i)+lamda_1*p_col_prev(i+1);
            }
            p_matrix.col(s)=p_col_new.segment(0, M); // Store in p_matrix
            p_col_prev= p_col_new;
        }
        else
        {
            int distance_till_end = r-1-s; // Distance to the end of a trigangle
            int new_col_size = M+distance_till_end; // Size of the new column vector
            p_col_new.setZero(new_col_size);
            for(int i=0;i<new_col_size;i++)
            {
                if(i==0)
                    p_col_new(i)=lamda_0*p_col_prev(i)+lamda_1*p_col_prev(i+1); // First layer
                else
                    p_col_new(i)=lamda_1*p_col_prev(i-1)+lamda_0*p_col_prev(i)+lamda_1*p_col_prev(i+1);
            }
            p_matrix.col(s)=p_col_new.segment(0, M); // Store in p_matrix
            p_col_prev = p_col_new;
        }
    }
}
void Scheutjens_Model::Calculate_avr_phi_i()
{
    avr_phi_i = lamda_matrix * phi_i;
    avr_phi_solv_i = lamda_matrix * phi_solv_i;
}
void Scheutjens_Model::Print_Vector(const Eigen::VectorXd& vec) {
    return; // Commented out to avoid excessive output
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << vec << std::endl;
}
void Scheutjens_Model::Print_Matrix(const Eigen::MatrixXd& mat) {
    return; // Commented out to avoid excessive output
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << mat << std::endl;
}