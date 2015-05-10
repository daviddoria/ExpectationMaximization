#include "ExpectationMaximization.h"
#include "GaussianModel.h"

#include <iostream>
#include <random>

Eigen::MatrixXd GenerateData(const unsigned int numberOfSamplesPerGroup, const unsigned int dimensionality);

void TestNDEvaluation();

void TestND();

int main(int, char*[])
{
  TestNDEvaluation();

  TestND();

  return EXIT_SUCCESS;
}

Eigen::MatrixXd GenerateData(const unsigned int numberOfSamplesPerGroup, const unsigned int dimensionality)
{
  // Mean 0, Standard deviation 1
  Eigen::MatrixXd data1(dimensionality, numberOfSamplesPerGroup);
  for(unsigned int i = 0; i < numberOfSamplesPerGroup; i++)
  {
    Eigen::VectorXd v = Eigen::VectorXd::Random(dimensionality);
    data1.col(i) = v;
  }

  // Mean 5, Standard deviation 1
  Eigen::MatrixXd data2(dimensionality, numberOfSamplesPerGroup);
  for(unsigned int i = 0; i < numberOfSamplesPerGroup; i++)
  {
    Eigen::VectorXd v = Eigen::VectorXd::Random(dimensionality);
    v = v.array() + 5;
    data2.col(i) = v;
  }

  // Mean -5, Standard deviation 1.5
  Eigen::MatrixXd data3(dimensionality, numberOfSamplesPerGroup);
  for(unsigned int i = 0; i < numberOfSamplesPerGroup; i++)
  {
      Eigen::VectorXd v = 1.5 * Eigen::VectorXd::Random(dimensionality);
      v = v.array() - 5;
      data3.col(i) = v;
  }

  // Concatentate the matrices horiztonally
  Eigen::MatrixXd data(dimensionality, data1.cols() + data2.cols() + data3.cols());
  data << data1, data2, data3;

  return data;
}

void TestNDEvaluation()
{
  GaussianModel p(2);

  Eigen::VectorXd mean(2);
  mean(0) = 0;
  mean(1) = 5;
  p.SetMean(mean);

  Eigen::MatrixXd variance = Eigen::MatrixXd::Identity(2,2);
  variance(0,0) = 2.0;

  p.SetVariance(variance);

  Eigen::VectorXd x(2);
  x(0) = 0.3;
  x(0) = 0.4;
  std::cout << "ND: " << p.Evaluate(x) << std::endl;
}

void TestND()
{
  int dimensionality = 2;

  // Generate some data
  Eigen::MatrixXd data = GenerateData(40, dimensionality);

  // Initialize the model
  std::vector<Model*> models(3);

  for(unsigned int i = 0; i < models.size(); i++)
  {
    Model* model = new GaussianModel(dimensionality);
    models[i] = model;
  }

  ExpectationMaximization expectationMaximization;

  expectationMaximization.SetData(data);
  expectationMaximization.SetModels(models);
  expectationMaximization.SetMinChange(1e-5);
  //expectationMaximization->SetInitializationTechniqueToRandom();
  expectationMaximization.SetInitializationTechniqueToKMeans();

  Eigen::VectorXd a(2);
  a(0) = -5;
  a(1) = -5;
  std::cout << "a: " << a << " likelihood: " << expectationMaximization.GetModel(0)->WeightedEvaluate(a) << std::endl;
  std::cout << "a: " << a << " likelihood: " << expectationMaximization.WeightedEvaluate(a) << std::endl;

  std::cout << "Final models:" << std::endl;
  //OutputModelInfo(expectationMaximization.GetModels());
  //PlotModels2D(expectationMaximization->GetModels(), data);
}
