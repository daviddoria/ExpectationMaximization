#ifndef GAUSSIANND_H
#define GAUSSIANND_H

/*
http://www.cs.mcgill.ca/~cs644/Godfried/2005/Fall/boutin/gaussian.html

p(x) = 1/((2 PI)^(d/2) SIGMA^(1/2)) e^{-1/2 (x-u)^T SIGMA^(-1)  (x - u)}

MLE: 

u = (1/N) sum_{n=1}^N (x_n)
SIGMA = (1/N) sum_{n=1}^N (x_n - u)(x_n - u)^T
*/

#include "Model.h"

class GaussianND : public Model
{
  public:
    virtual double Evaluate(Eigen::VectorXd x);

    virtual double WeightedEvaluate(Eigen::VectorXd x);
};

#endif
