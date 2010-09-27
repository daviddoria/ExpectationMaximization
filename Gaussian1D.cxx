#include "Gaussian1D.h"

Gaussian1D::Gaussian1D()
{
  this->Variance.set_size(1,1);
  this->Variance(0,0) = 1.0;
  
  this->Mean.set_size(1);
  this->Mean(0) = 0.0;
  
  this->Dimensionality = 1;
}

void OutputModelInfo(std::vector<Model> models)
{
  for(unsigned int i = 0; i < models.size(); i++)
    {
    std::cout << "Model " << i << " : Mean = " << models[i].GetMean()
              << " Variance = " << models[i].GetVariance()
              << " Mixing coefficient = " << models[i].GetMixingCoefficient() << std::endl;
    }
}

double Gaussian1D::Evaluate(double x)
{
  /*
  vnl_vector<double> xvec(1);
  xvec(0) = x;
  return Evaluate(xvec);
  */
  
  return (1./sqrt(2*vtkMath::Pi()*this->Variance(0,0))) * exp(-(pow(x-this->Mean(0),2))/(2.*this->Variance(0,0)));
}