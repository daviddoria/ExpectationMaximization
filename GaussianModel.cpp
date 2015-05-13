/*
Copyright (C) 2015 David Doria, daviddoria@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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
