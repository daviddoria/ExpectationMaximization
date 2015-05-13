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

#include "Model.h"

void Model::SetDimensionality(const unsigned int dim)
{
  this->Dimensionality = dim;
}

unsigned int Model::GetDimensionality()
{
  return this->Dimensionality;
}

Eigen::VectorXd Model::GetMean() const
{
  return this->Mean;
}

void Model::SetMean(const Eigen::VectorXd& m)
{
  this->Mean = m;
}

Eigen::MatrixXd Model::GetVariance() const
{
  return this->Variance;
}

void Model::SetVariance(const Eigen::MatrixXd& v)
{
  this->Variance = v;
}

double Model::GetMixingCoefficient() const
{
  return this->MixingCoefficient;
}

void Model::SetMixingCoefficient(double m)
{
  this->MixingCoefficient = m;
}
