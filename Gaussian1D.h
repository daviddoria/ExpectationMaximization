#ifndef GAUSSIAN1D_H
#define GAUSSIAN1D_H

/*

p(x) = 1/(sigma sqrt(2 pi)) e^{-(x-u)^2 / (2 sigma^2)}

MLE: 

u = (1/N) sum_{n=1}^N (x_n)
sigma^2 = (1/N) sum_{n=1}^N (x_n - u)^2
*/

#include "GaussianND.h"

#include <vtkMath.h>

#include <cmath>
#include <vector>

class Gaussian1D : public GaussianND
{
  public:
    
    Gaussian1D();
    
    double GetMean(){return this->Mean(0);}
    void SetMean(double m){this->Mean(0) = m;}

    double GetVariance(){return this->Variance(0,0);}
    void SetVariance(double v){this->Variance(0,0) = v;}

    using GaussianND::Evaluate;
    double Evaluate(double x);

    using GaussianND::WeightedEvaluate;
    double WeightedEvaluate(double x)
    {
      vnl_vector<double> xvec(1);
      xvec(0) = x;
      return WeightedEvaluate(xvec);
    }

};

#endif