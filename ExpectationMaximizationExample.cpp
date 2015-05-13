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
