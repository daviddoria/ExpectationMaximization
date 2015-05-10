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
