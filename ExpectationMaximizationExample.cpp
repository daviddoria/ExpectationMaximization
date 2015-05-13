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

#include "ExpectationMaximization.h"
#include "GaussianModel.h"

#include <iostream>
#include <random>
#include <chrono>

Eigen::MatrixXd GenerateData(const unsigned int numberOfSamplesPerGroup, const unsigned int dimensionality);
long int GetSeed(const bool random);

int main(int, char*[])
{
    int dimensionality = 1;

    // Generate some data
    Eigen::MatrixXd data = GenerateData(1000, dimensionality);

    // Initialize the model
    std::vector<Model*> models(2);

    for(unsigned int i = 0; i < models.size(); i++)
    {
      Model* model = new GaussianModel(dimensionality);
      models[i] = model;
    }

    MixtureModel mixtureModel;
    mixtureModel.SetModels(models);

    ExpectationMaximization expectationMaximization;
    expectationMaximization.SetData(data);
    expectationMaximization.SetMixtureModel(mixtureModel);
    expectationMaximization.SetMinChange(1e-5);
    expectationMaximization.SetMaxIterations(100);

    expectationMaximization.Compute();

    std::cout << "Final models:" << std::endl;
    MixtureModel finalModel = expectationMaximization.GetMixtureModel();
    for(unsigned int i = 0; i < finalModel.GetNumberOfModels(); ++i)
    {
      finalModel.GetModel(i)->Print();
    }

  return EXIT_SUCCESS;
}

Eigen::MatrixXd GenerateData(const unsigned int numberOfSamplesPerGroup, const unsigned int dimensionality)
{
  std::default_random_engine generator(GetSeed(true));

  // Mean 0, variance 1
  Eigen::MatrixXd data1(dimensionality, numberOfSamplesPerGroup);
  std::normal_distribution<double> normalDistribution1(0,sqrt(1));
  for(unsigned int i = 0; i < numberOfSamplesPerGroup; i++)
  {
    Eigen::VectorXd v(1);
    v(0) = normalDistribution1(generator);
    data1.col(i) = v;
  }

  // Mean 5, variance 2
  Eigen::MatrixXd data2(dimensionality, numberOfSamplesPerGroup);
  std::normal_distribution<double> normalDistribution2(5,sqrt(2));
  for(unsigned int i = 0; i < numberOfSamplesPerGroup; i++)
  {
    Eigen::VectorXd v(1);
    v(0) = normalDistribution2(generator);
    data2.col(i) = v;
  }

  // Concatentate the matrices horiztonally
  Eigen::MatrixXd data(dimensionality, data1.cols() + data2.cols());
  data << data1, data2;

  return data;
}


long int GetSeed(const bool random)
{
  if(random)
  {
    auto now = std::chrono::system_clock::now();
    return now.time_since_epoch().count();
  }
  else
  {
    return 0;
  }
}
