#ifndef MODEL_H
#define MODEL_H

#include <cmath>
#include <vector>

#include <Eigen/Dense>

class Model
{
  public:
    
    void Init();
    
    void SetDiagonalCovariance(std::vector<double> diagonal);
    
    void SetDimensionality(unsigned int dim);
    unsigned int GetDimensionality();
    
    Eigen::VectorXd GetMean();
    void SetMean(Eigen::VectorXd m);

    Eigen::MatrixXd GetVariance();
    void SetVariance(Eigen::MatrixXd v);

    double GetMixingCoefficient();
    void SetMixingCoefficient(double m);

    virtual double Evaluate(Eigen::VectorXd x) = 0;

    virtual double WeightedEvaluate(Eigen::VectorXd x) = 0;
    
  protected:

    Eigen::VectorXd Mean;

    Eigen::MatrixXd Variance;
    
    unsigned int Dimensionality;
    
    double MixingCoefficient;
   
   
};

#endif
