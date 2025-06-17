#include "Eigen/Dense"

class Model{
    public:
    bool initialized;
    void Init();
    double function(double x);
    //void function(double* xs, double* ys, int n);
    void function(Eigen::VectorXd&  xs,   Eigen::VectorXd&  ys, int n);
};