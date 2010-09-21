#ifndef MODEL_H
#define MODEL_H

#include <vtkMath.h>

#include <cmath>
#include <vector>

class Model
{
  public:
    double GetMean(){return this->Mean;}
    void SetMean(double m){this->Mean = m;}

    double GetVariance(){return this->Variance;}
    void SetVariance(double v){this->Variance = v;}

    double GetStandardDeviation(){return sqrt(this->Variance);}
    void SetStandardDeviation(double s){this->Variance = pow(s,2);}

    double GetMixingCoefficient(){return this->MixingCoefficient;}
    void SetMixingCoefficient(double m){this->MixingCoefficient = m;}

    double Evaluate(double x)
    {
      return (1./sqrt(2*vtkMath::Pi()*this->Variance)) * exp(-(pow(x-this->Mean,2))/(2.*this->Variance));
    }

    double WeightedEvaluate(double x)
    {
      return (this->MixingCoefficient/sqrt(2*vtkMath::Pi()*this->Variance)) * exp(-(pow(x-this->Mean,2))/(2.*this->Variance));
    }
  private:
    double Mean;
    double Variance;
    double MixingCoefficient;
};

void PlotModels(std::vector<Model> models, double range[2]);
void PlotModels(std::vector<Model> models, std::vector<double> points);
void OutputModelInfo(std::vector<Model> models);

#endif