#include "ExpectationMaximization.h"
#include "Gaussian1D.h"

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
    std::vector<Model*> models(3);

    bool random = true;

    std::default_random_engine generator(GetSeed(random));

    std::uniform_real_distribution<double> meanDistribution(-5.0,5.0);
    std::uniform_real_distribution<double> varianceDistribution(0,3.0);

    for(unsigned int i = 0; i < models.size(); i++)
    {
      Model* model = new Gaussian1D;
      model->SetDimensionality(dimensionality);

      dynamic_cast<Gaussian1D*>(model)->SetMean(meanDistribution(generator));

      dynamic_cast<Gaussian1D*>(model)->SetVariance(varianceDistribution(generator));
      model->SetMixingCoefficient(1./models.size());
      models[i] = model;
    }

    ExpectationMaximization expectationMaximization;
    expectationMaximization.SetData(data);
    expectationMaximization.SetModels(models);
    expectationMaximization.SetMinChange(1e-5);
    expectationMaximization.SetMaxIterations(100);

    expectationMaximization.Compute();

    std::cout << "Final models:" << std::endl;
    for(unsigned int i = 0; i < expectationMaximization.GetNumberOfModels(); ++i)
    {
      expectationMaximization.GetModel(i)->Print();
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

  // Mean 5, variance 1
  Eigen::MatrixXd data2(dimensionality, numberOfSamplesPerGroup);
  std::normal_distribution<double> normalDistribution2(5,sqrt(1));
  for(unsigned int i = 0; i < numberOfSamplesPerGroup; i++)
  {
    Eigen::VectorXd v(1);
    v(0) = normalDistribution2(generator);
    data2.col(i) = v;
  }

  // Mean -5, variance 1.5
  Eigen::MatrixXd data3(dimensionality, numberOfSamplesPerGroup);
  std::normal_distribution<double> normalDistribution3(-5,sqrt(1.5));
  for(unsigned int i = 0; i < numberOfSamplesPerGroup; i++)
  {
      Eigen::VectorXd v(1);
      v(0) = normalDistribution3(generator);
      data3.col(i) = v;
  }

  // Concatentate the matrices horiztonally
  Eigen::MatrixXd data(dimensionality, data1.cols() + data2.cols() + data3.cols());
  data << data1, data2, data3;

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
