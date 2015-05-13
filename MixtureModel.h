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

    unsigned int GetNumberOfModels() const
    {
        return this->Models.size();
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

    Model* GetModel(const unsigned int i) const
    {
        return Models[i];
    }

    void CopyFrom(const MixtureModel& mixtureModel)
    {
        this->Models.resize(mixtureModel.GetNumberOfModels());

        for(unsigned int i = 0; i < this->Models.size(); ++i)
        {
            this->Models[i] = mixtureModel.GetModel(i)->Clone();
        }
    }

    unsigned int GetDimensionality() const
    {
        return this->Models[0]->GetDimensionality();
    }

  protected:

    std::vector<Model*> Models;
};

#endif
