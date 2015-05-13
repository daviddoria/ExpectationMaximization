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

#ifndef GAUSSIANMODEL_H
#define GAUSSIANMODEL_H

/*
This is a subclass of Model that implements the Gaussian distribution.

For more information, see http://www.cs.mcgill.ca/~cs644/Godfried/2005/Fall/boutin/gaussian.html

p(x) = 1/((2 PI)^(d/2) SIGMA^(1/2)) e^{-1/2 (x-u)^T SIGMA^(-1)  (x - u)}

Maximum Likelihood Estimates:

u = (1/N) sum_{n=1}^N (x_n)
SIGMA = (1/N) sum_{n=1}^N (x_n - u)(x_n - u)^T
*/

#include "Model.h"

class GaussianModel : public Model
{
  public:
    GaussianModel(const unsigned int dimensionality);

    virtual double Evaluate(const Eigen::VectorXd& x) const;

    virtual double WeightedEvaluate(const Eigen::VectorXd& x) const;

    void Print() const;

    virtual GaussianModel* Clone()
    {
        GaussianModel* newObject = new GaussianModel(this->Dimensionality);

        newObject->SetMean(this->Mean);
        newObject->SetVariance(this->Variance);
        newObject->SetMixingCoefficient(this->MixingCoefficient);

        return newObject;
    }
};

#endif
