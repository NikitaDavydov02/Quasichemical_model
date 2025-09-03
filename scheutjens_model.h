#include "Eigen/Dense"

class Scheutjens_Model{
    public:
    
    void Init();
    void InitFromFile();
    void function(Eigen::VectorXd&  xs,   Eigen::VectorXd&  ys, int n);
    void Calculate_avr_phi_i();
    void Print_Vector(const Eigen::VectorXd& vec);
    void Print_Matrix(const Eigen::MatrixXd& mat);
    void Calculate_p_matrix();
    void Calculate_phi_i_dashed();
    void Calculate_ln_p_i();

    int max_iter;
    double required_error;
    double newthon_break_coeff;

    bool initialized;
    bool new_method;
    bool output;
    int M;
    int r;
    double chi;   //Flory parameter
    double chi_s; //Adsorption_parameter
    double phi_bulk[2]; //bulk molar fractions [0] - poly, [1] - solvent
    Eigen::VectorXd phi_i; //molar fractions of polymer in each layer
    Eigen::VectorXd phi_solv_i;
    Eigen::VectorXd phi_i_dashed;
    double lnP_bulk;
    double lamda_0;
    double lamda_1;

    Eigen::VectorXd  avr_phi_i; //average molar fraction of polymer in each layer
    Eigen::VectorXd  avr_phi_solv_i; //average molar fraction of solvent in each layer

    Eigen::MatrixXd lamda_matrix;
    Eigen::MatrixXd p_matrix;
    Eigen::MatrixXd w_matrix;

};