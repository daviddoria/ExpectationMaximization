#ifndef MODEL_H
#define MODEL_H

#include <vtkMath.h>
#include <vtkSmartPointer.h>
#include <vtkDenseArray.h>

#include <cmath>
#include <vector>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_determinant.h>
#include <vnl/vnl_inverse.h>

class Model
{
  public:
    
    void Init();
    
    void SetDiagonalCovariance(std::vector<double> diagonal);
    
    void SetDimensionality(int dim);
    int GetDimensionality();
    
    vnl_vector<double> GetMean();
    void SetMean(vnl_vector<double> m);

    vnl_matrix<double> GetVariance();
    void SetVariance(vnl_matrix<double> v);

    double GetMixingCoefficient();
    void SetMixingCoefficient(double m);

    virtual double Evaluate(vnl_vector<double> x) = 0;

    virtual double WeightedEvaluate(vnl_vector<double> x) = 0;
    
  protected:

    vnl_vector<double> Mean;

    vnl_matrix<double> Variance;
    
    int Dimensionality;
    
    double MixingCoefficient;
   
   
};

#endif