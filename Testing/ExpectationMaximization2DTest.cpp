#include "ExpectationMaximization.h"
#include "GaussianModel.h"

#include <iostream>
#include <random>

Eigen::MatrixXd GenerateData(const unsigned int numberOfSamplesPerGroup, const unsigned int dimensionality);

bool Test2DEvaluation();

bool Test2D_EM();

int main(int, char*[])
{
  bool passEval = Test2DEvaluation();
  std::cout << "passEval? " << passEval << std::endl;

  bool passEM = Test2D_EM();
  std::cout << "passEM? " << passEM << std::endl;

  return EXIT_SUCCESS;
}

bool Test2DEvaluation()
{
  unsigned int dimensionality = 2;
  GaussianModel p(dimensionality);

  Eigen::VectorXd mean(dimensionality);
  mean(0) = 0;
  mean(1) = 5;
  p.SetMean(mean);

  Eigen::MatrixXd variance = Eigen::MatrixXd::Identity(dimensionality, dimensionality);
  variance(0,0) = 2.0;

  p.SetVariance(variance);

  Eigen::VectorXd x(dimensionality);
  x(0) = 0.1;
  x(1) = 5.1;

  double eval = p.Evaluate(x);
  std::cout << "2D eval: " << eval << std::endl;

  double test = 0.111699;

  double epsilon = 1e-5;
  if(fabs(eval - test) > epsilon)
  {
      return false;
  }

  return true;
}

bool Test2D_EM()
{
  int dimensionality = 2;

  // Generate some data
  Eigen::MatrixXd data = GenerateData(40, dimensionality);

  // Initialize the model
  std::vector<Model*> models(2);

  for(unsigned int i = 0; i < models.size(); i++)
  {
    Model* model = new GaussianModel(dimensionality);
    models[i] = model;
  }

  ExpectationMaximization expectationMaximization;

  expectationMaximization.SetData(data);
  expectationMaximization.SetRandom(false);
  expectationMaximization.SetModels(models);
  expectationMaximization.SetMaxIterations(3);
  expectationMaximization.SetInitializationTechniqueToKMeans();
  expectationMaximization.Compute();

  // This is where we got the test output
  for(unsigned int i = 0; i < expectationMaximization.GetNumberOfModels(); ++i)
  {
    std::cout << "Model " << i << ":" << std::endl;
    expectationMaximization.GetModel(i)->Print();
  }

  Eigen::VectorXd mean0(2);
  mean0 << 0.0631675, -0.147669;

  Eigen::VectorXd mean1(2);
  mean1 << 10.23, 10.069;

  Eigen::MatrixXd var0(2,2);
  var0 << 0.818243, -0.027182,
          -0.027182,  0.836947;

  Eigen::MatrixXd var1(2,2);
  var1 << 2.49997, 0.10499,
          0.10499, 1.85563;

  double epsilon = 1e-4;

  // Check means
  if((expectationMaximization.GetModel(0)->GetMean() - mean0).norm() > epsilon)
  {
      return false;
  }
  if((expectationMaximization.GetModel(1)->GetMean() - mean1).norm() > epsilon)
  {
      return false;
  }

  // Check variances
  if((expectationMaximization.GetModel(0)->GetVariance() - var0).norm() > epsilon)
  {
      return false;
  }
  if((expectationMaximization.GetModel(1)->GetVariance() - var1).norm() > epsilon)
  {
      return false;
  }

  return true;
}


Eigen::MatrixXd GenerateData(const unsigned int numberOfSamplesPerGroup, const unsigned int dimensionality)
{
  std::default_random_engine generator(0); // Seed with 0, so results are repeatable

  // Mean 0, Standard deviation 1 for both dimensions
  Eigen::MatrixXd data1(dimensionality, numberOfSamplesPerGroup);
  std::normal_distribution<double> normalDistribution1(0,sqrt(1));
  for(unsigned int i = 0; i < numberOfSamplesPerGroup; i++)
  {
    Eigen::VectorXd v(dimensionality);
    v(0) = normalDistribution1(generator);
    v(1) = normalDistribution1(generator);
    data1.col(i) = v;
  }

  // Mean 10, Standard deviation 2 for both dimensions
  Eigen::MatrixXd data2(dimensionality, numberOfSamplesPerGroup);
  std::normal_distribution<double> normalDistribution2(10,sqrt(2));
  for(unsigned int i = 0; i < numberOfSamplesPerGroup; i++)
  {
    Eigen::VectorXd v(dimensionality);
    v(0) = normalDistribution2(generator);
    v(1) = normalDistribution2(generator);
    data2.col(i) = v;
  }

  // Concatentate the matrices horiztonally
  Eigen::MatrixXd data(dimensionality, data1.cols() + data2.cols());
  data << data1, data2;

  return data;
}
