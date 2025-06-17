#include "Model.h"

void Model::Init() {
    // Initialization logic here
    initialized = true; // Set initialized to true after initialization
}
double Model::function(double x)
{
    return x*x-1;
}
/*void Model::function(double* xs, double* ys, int n)
{
    ys[0]=xs[0]*xs[0]+xs[1]*xs[1]-4;
    ys[1]=xs[0]*xs[0]*xs[1]-1;
}*/
void Model::function( Eigen::VectorXd& xs,  Eigen::VectorXd& ys, int n)
{
    //ys(0)=xs(0)*xs(0)+xs(1)*xs(1)-4;
    //ys(1)=xs(0)*xs(0)*xs(1)-1;
    int x0=1;
    int y0=2;
    int z0=3;
    ys(0) = xs(0) * xs(0) + xs(1) * xs(1) - 4;
    ys(1) = (xs(0)-4) *(xs(0)-4)  + xs(1) * xs(1) - 4;
    if(n==3)
    {
        //ys(0) = (xs(2)-4-z0) *(xs(2)-4-z0)  + (xs(0)-x0) * (xs(0)-x0)+ (xs(1)-y0) * (xs(1)-y0) - 4;
        //ys(1) = (xs(0)-4-x0) *(xs(0)-4-x0)  + (xs(1)-y0) *(xs(1)-y0)+ (xs(2)-z0) *(xs(2)-z0) - 4;
        //ys(2) = (xs(1)-4-y0) *(xs(1)-4-y0)  + (xs(0)-x0) * (xs(0)-x0)+ (xs(2)-z0) * (xs(2)-z0) - 4;
        ys(0) = (xs(2)) *(xs(2))  + (xs(0)) * (xs(0))+ (xs(1)) * (xs(1)) - 25;
        ys(1)= xs(0)-3;
        ys(2)=xs(1)-4;
    }
}