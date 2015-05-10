#include "Model.h"

void Model::SetDiagonalCovariance(const Eigen::VectorXd& diagonal)
{
  this->Variance.diagonal() = diagonal;
}

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

void Model::SetMean(Eigen::VectorXd m)
{
  this->Mean = m;
}

Eigen::MatrixXd Model::GetVariance() const
{
  return this->Variance;
}

void Model::SetVariance(Eigen::MatrixXd v)
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
