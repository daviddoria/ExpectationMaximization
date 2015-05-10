#include "Model.h"

void Model::Init()
{
  // Initialize mean
  this->Mean.resize(this->Dimensionality);
  for(unsigned int i = 0; i < this->Dimensionality; i++)
  {
    this->Mean(i) = 0;
  }
  
  // Initialize variance
  this->Variance.resize(this->Dimensionality, this->Dimensionality);
  Eigen::VectorXd diagonal(this->Dimensionality);
  for(unsigned int i = 0; i < this->Dimensionality; i++)
  {
    diagonal(i) = i;
  }
    
  this->SetDiagonalCovariance(diagonal);
}

void Model::SetDiagonalCovariance(const Eigen::VectorXd& diagonal)
{
  this->Variance.diagonal() = diagonal;
//  for(unsigned int i = 0; i < this->Dimensionality; i++)
//  {
//    for(unsigned int j = 0; j < this->Dimensionality; j++)
//    {
//      if(i == j)
//      {
//        this->Variance(i,j) = diagonal[i];
//      }
//      else
//      {
//        this->Variance(i,j) = 0;
//      }
//    }
//  }
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
