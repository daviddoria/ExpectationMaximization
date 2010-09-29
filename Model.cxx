#include "Model.h"

void Model::Init()
{
  // Initialize mean
  this->Mean.set_size(this->Dimensionality);
  for(unsigned int i = 0; i < this->Dimensionality; i++)
    {
    this->Mean(i) = 0;
    }
  
  // Initialize variance
  this->Variance.set_size(this->Dimensionality, this->Dimensionality);
  std::vector<double> diagonal;
  for(unsigned int i = 0; i < this->Dimensionality; i++)
    {
    diagonal.push_back(i);
    }
    
  this->SetDiagonalCovariance(diagonal);
}

void Model::SetDiagonalCovariance(std::vector<double> diagonal)
{
  for(unsigned int i = 0; i < this->Dimensionality; i++)
    {
    for(unsigned int j = 0; j < this->Dimensionality; j++)
      {
      if(i == j)
	{
	this->Variance(i,j) = diagonal[i];
	}
      else
	{
	this->Variance(i,j) = 0;
	}
      }
    }
}

void Model::SetDimensionality(unsigned int dim)
{
  this->Dimensionality = dim;
}

unsigned int Model::GetDimensionality()
{
  return this->Dimensionality;
}

vnl_vector<double> Model::GetMean()
{
  return this->Mean;
}

void Model::SetMean(vnl_vector<double> m)
{
  this->Mean = m;
}

vnl_matrix<double> Model::GetVariance()
{
  return this->Variance;
}

void Model::SetVariance(vnl_matrix<double> v)
{
  this->Variance = v;
}


double Model::GetMixingCoefficient()
{
  return this->MixingCoefficient;
}

void Model::SetMixingCoefficient(double m)
{
  this->MixingCoefficient = m;
}