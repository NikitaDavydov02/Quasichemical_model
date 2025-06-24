#include "Eigen/Dense"

class Scheutjens_Model{
    public:
    
    void Init();
    void function(Eigen::VectorXd&  xs,   Eigen::VectorXd&  ys, int n);
    void Calculate_avr_phi_i();
    void Print_Vector(const Eigen::VectorXd& vec);
    void Print_Matrix(const Eigen::MatrixXd& mat);

    bool initialized;
    int M;
    int z;
    int r;
    double chi;   //Flory parameter
    double chi_s; //Adsorption_parameter
    double phi_bulk[2]; //bulk molar fractions [0] - poly, [1] - solvent
    Eigen::VectorXd phi_i; //molar fractions of polymer in each layer
    Eigen::VectorXd phi_solv_i;
    Eigen::VectorXd phi_i_dashed;
    double lnP_bulk;

    Eigen::VectorXd  avr_phi_i; //average molar fraction of polymer in each layer
    Eigen::VectorXd  avr_phi_solv_i; //average molar fraction of solvent in each layer

    Eigen::MatrixXd lamda_matrix;
    Eigen::MatrixXd p_matrix;
    Eigen::MatrixXd w_matrix;

};