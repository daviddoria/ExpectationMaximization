#include <cmath>
#include <iostream>
#include <vector>

#include "GaussianModel.h"

GaussianModel::GaussianModel(const unsigned int dimensionality)
{
  this->Dimensionality = dimensionality;
  this->Variance = Eigen::MatrixXd::Identity(dimensionality, dimensionality);

  this->Mean = Eigen::VectorXd::Zero(dimensionality);
}

double GaussianModel::Evaluate(const Eigen::VectorXd& x) const
{
  // p(x) = 1/((2 PI)^(d/2) |SIGMA|^(1/2)) e^{-1/2 (x-u)^T SIGMA^(-1) (x - u)}
  double denom = pow(2.0*M_PI, this->Dimensionality/2.0)*sqrt(this->Variance.determinant());
  return (1.0/denom) * exp(-1.0/2.0 * (x-this->Mean).transpose() * this->Variance.inverse() * (x-this->Mean));
}

double GaussianModel::WeightedEvaluate(const Eigen::VectorXd& x) const
{
  return this->MixingCoefficient * Evaluate(x);
}

void GaussianModel::Print() const
{
    std::cout << "Mean: " << this->Mean << std::endl;
    std::cout << "Variance: " << this->Variance << std::endl;
}
