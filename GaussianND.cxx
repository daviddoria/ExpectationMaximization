/*
http://www.cs.mcgill.ca/~cs644/Godfried/2005/Fall/boutin/gaussian.html

p(x) = 1/((2 PI)^(d/2) SIGMA^(1/2)) e^{-1/2 (x-u)^T SIGMA^(-1)  (x - u)}

MLE: 

u = (1/N) sum_{n=1}^N (x_n)
SIGMA = (1/N) sum_{n=1}^N (x_n - u)(x_n - u)^T
*/

#include <vtkMath.h>
#include <vtkSmartPointer.h>
#include <vtkDenseArray.h>

#include <cmath>
#include <vector>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_determinant.h>
#include <vnl/vnl_inverse.h>

#include "GaussianND.h"

double GaussianND::Evaluate(vnl_vector<double> x)
{
  // p(x) = 1/((2 PI)^(d/2) |SIGMA|^(1/2)) e^{-1/2 (x-u)^T SIGMA^(-1)  (x - u)}
  double denom = pow(2*vtkMath::Pi(), this->Dimensionality/2.0)*sqrt(vnl_determinant(this->Variance));
  return (1./denom) * exp(-1.0/2.0 * dot_product((x-this->Mean), vnl_inverse(this->Variance) * (x-this->Mean)));
}

double GaussianND::WeightedEvaluate(vnl_vector<double> x)
{
  return this->MixingCoefficient * Evaluate(x);
}
