#include "Model.h"

void Model::Init() {
    // Initialization logic here
    initialized = true; // Set initialized to true after initialization
}
double Model::function(double x)
{
    return x*x-1;
}