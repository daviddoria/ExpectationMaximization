#include "Model.h"

void Model::Init()
{
  //this->Variance->Resize(this->Dimensionality, this->Dimensionality);
	
  this->Variance.set_size(this->Dimensionality, this->Dimensionality);
  this->Mean.set_size(this->Dimensionality);
  
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
    //this->Variance->SetValue(i,i,diagonal[i]);
    this->Variance(i,i) = diagonal[i];
    }
}

void Model::SetDimensionality(int dim)
{
  this->Dimensionality = dim;
}

int Model::GetDimensionality()
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