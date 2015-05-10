#include "ExpectationMaximization.h"
#include "GaussianND.h"

#include <iostream>
#include <random>

Eigen::MatrixXd GenerateData(const unsigned int numberOfSamplesPerGroup, const unsigned int dimensionality);

void Test1DEvaluation();
void TestNDEvaluation();

void Test1D();
void TestND();

int main(int, char*[])
{
  //Test1DEvaluation();
  //TestNDEvaluation();

  //Test1D();
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


void Test1DEvaluation()
{
  Gaussian1D p;
  p.SetDimensionality(1);
  p.SetMean(0);
  p.SetVariance(2);
  std::cout << "1D: " << p.Evaluate(.3) << std::endl;
}

void TestNDEvaluation()
{
  GaussianND p;
  p.SetDimensionality(1);

  Eigen::VectorXd mean(1);
  mean(0) = 0;
  p.SetMean(mean);

  Eigen::MatrixXd variance(1,1);
  variance(0,0) = 2.0;

  p.SetVariance(variance);

  Eigen::VectorXd x(1);
  x(0) = 0.3;
  std::cout << "ND: " << p.Evaluate(x) << std::endl;
}

void Test1D()
{
  int dimensionality = 1;

  // Generate some data
  Eigen::MatrixXd data = GenerateData(40, dimensionality);

  // Initialize the model
  std::vector<Model*> models(3);

  std::default_random_engine generator;
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

  // Display initial model
  std::cout << "Randomly initialized model:" << std::endl;
  //OutputModelInfo(models);
  //double range[2] = {-15, 15};
  //PlotModels1D(models, range);

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

void TestND()
{
  int dimensionality = 2;

  // Generate some data
  Eigen::MatrixXd data = GenerateData(40, dimensionality);

  // Initialize the model
  std::vector<Model*> models(3);

  for(unsigned int i = 0; i < models.size(); i++)
  {
    Model* model = new GaussianND;
    model->SetDimensionality(dimensionality);
    model->Init();
    models[i] = model;
  }

  // Display initial model
  //std::cout << "Randomly initialized model:" << std::endl;
  //OutputModelInfo(models);
  //PlotModels2D(models, data);

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
