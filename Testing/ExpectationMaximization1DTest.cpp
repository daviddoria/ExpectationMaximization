#include "ExpectationMaximization.h"
#include "GaussianModel.h"

#include <iostream>
#include <random>

Eigen::MatrixXd GenerateData(const unsigned int numberOfSamplesPerGroup, const unsigned int dimensionality);

void Test1DEvaluation();

void Test1D();

int main(int, char*[])
{
  Test1DEvaluation();

  Test1D();

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


void Test1DEvaluation()
{
  GaussianModel p(1);
  Eigen::VectorXd mean = Eigen::VectorXd::Zero(1);
  p.SetMean(mean);
  Eigen::VectorXd variance = Eigen::MatrixXd::Identity(1,1);
  p.SetVariance(variance);
  Eigen::VectorXd v(1);
  v(0) = .3;
  std::cout << "1D: " << p.Evaluate(v) << std::endl;
}

void Test1D()
{
  int dimensionality = 1;

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

  std::cout << "Final models:" << std::endl;
  //OutputModelInfo(expectationMaximization->GetModels());
  //PlotModels1D(expectationMaximization->GetModels(), range);

  Eigen::VectorXd v(1);
  v(0) = 2;
  std::cout << "WeightedEvaluate(2) = " << expectationMaximization.WeightedEvaluate(v) << std::endl;
}
