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

#ifndef MODEL_H
#define MODEL_H

#include <cmath>
#include <vector>

#include <Eigen/Dense>

class Model
{
  public:
    /** The dimensionality of the vector space. */
    void SetDimensionality(const unsigned int dim);
    unsigned int GetDimensionality();
    
    /** The mean vector of the model. */
    Eigen::VectorXd GetMean() const;
    void SetMean(const Eigen::VectorXd& m);

    /** The covariance matrix of the model. */
    Eigen::MatrixXd GetVariance() const;
    void SetVariance(const Eigen::MatrixXd& v);

    /** The mixing coefficient (the weight of this distribution relative to other distributions). */
    double GetMixingCoefficient() const;
    void SetMixingCoefficient(const double m);

    /** Compute the probability of a point according to this model. */
    virtual double Evaluate(const Eigen::VectorXd& x) const = 0;

    /** Compute the probability of a point according to this model, weighted by its mixing coefficient. */
    virtual double WeightedEvaluate(const Eigen::VectorXd& x) const = 0;

    /** Display properties of the model. */
    virtual void Print() const = 0;

    virtual Model* Clone() = 0;

  protected:

    Eigen::VectorXd Mean;

    Eigen::MatrixXd Variance;
    
    unsigned int Dimensionality;
    
    double MixingCoefficient;
};

#endif
