#ifndef MIXTUREMODEL_H
#define MIXTUREMODEL_H

#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include "Model.h"

class MixtureModel
{
  public:

    void SetModels(const std::vector<Model*>& models)
    {
        this->Models = models;
    }

    /** Compute the probability of a point according to each model, weighted by its mixing coefficient. */
    float WeightedEvaluate(const Eigen::VectorXd& x) const
    {
        float mixtureSum = 0;
        for(unsigned int i = 0; i < this->Models.size(); ++i)
        {
            mixtureSum += this->Models[i]->WeightedEvaluate(x);
        }
        return mixtureSum;
    }

  protected:

    std::vector<Model*> Models;
};

#endif
