#ifndef MODEL_H
#define MODEL_H

#include <cmath>
#include <vector>

#include <Eigen/Dense>

class Model
{
  public:

    void SetDiagonalCovariance(const Eigen::VectorXd& diagonal);
    
    void SetDimensionality(const unsigned int dim);
    unsigned int GetDimensionality();
    
    Eigen::VectorXd GetMean() const;
    void SetMean(Eigen::VectorXd m);

    Eigen::MatrixXd GetVariance() const;
    void SetVariance(Eigen::MatrixXd v);

    double GetMixingCoefficient() const;
    void SetMixingCoefficient(const double m);

    virtual double Evaluate(const Eigen::VectorXd& x) const = 0;

    virtual double WeightedEvaluate(const Eigen::VectorXd& x) const = 0;

    virtual void Print() const = 0;
    
  protected:

    Eigen::VectorXd Mean;

    Eigen::MatrixXd Variance;
    
    unsigned int Dimensionality;
    
    double MixingCoefficient;
};

#endif
