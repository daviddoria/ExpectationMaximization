/*
Copyright (C) 2015 David Doria, daviddoria@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MIXTUREMODEL_H
#define MIXTUREMODEL_H

#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include "Helpers/Helpers.h"

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

        if(mixtureSum == 0)
        {
            std::cout << "mixtureSum = 0!" << std::endl;
        }

        if(Helpers::IsNaN(mixtureSum))
        {
            std::cout << "mixtureSum is NaN!" << std::endl;
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
