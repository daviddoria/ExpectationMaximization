#ifndef MODEL_H
#define MODEL_H

#include <cmath>
#include <vector>

#include <Eigen/Dense>

class Model
{
  public:
    /** The dimensionality of the vector space. */
    void SetDimensionality(const unsigned int dim);
    unsigned int GetDimensionality();
    
    /** The mean vector of the model. */
    Eigen::VectorXd GetMean() const;
    void SetMean(const Eigen::VectorXd& m);

    /** The covariance matrix of the model. */
    Eigen::MatrixXd GetVariance() const;
    void SetVariance(const Eigen::MatrixXd& v);

    /** The mixing coefficient (the weight of this distribution relative to other distributions). */
    double GetMixingCoefficient() const;
    void SetMixingCoefficient(const double m);

    /** Compute the probability of a point according to this model. */
    virtual double Evaluate(const Eigen::VectorXd& x) const = 0;

    /** Compute the probability of a point according to this model, weighted by its mixing coefficient. */
    virtual double WeightedEvaluate(const Eigen::VectorXd& x) const = 0;

    /** Display properties of the model. */
    virtual void Print() const = 0;

    virtual Model* Clone() = 0;

  protected:

    Eigen::VectorXd Mean;

    Eigen::MatrixXd Variance;
    
    unsigned int Dimensionality;
    
    double MixingCoefficient;
};

#endif
