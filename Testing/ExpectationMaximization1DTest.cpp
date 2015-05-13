#include "ExpectationMaximization.h"
#include "GaussianModel.h"

#include <iostream>
#include <random>
#include <chrono>

Eigen::MatrixXd GenerateData(const unsigned int numberOfSamplesPerGroup, const unsigned int dimensionality);

bool Test1DEvaluation();

bool Test1D_EM();

int main(int, char*[])
{
  bool passEval = Test1DEvaluation();
  std::cout << "passEval? " << passEval << std::endl;

  bool passEM = Test1D_EM();
  std::cout << "passEM? " << passEM << std::endl;

  return EXIT_SUCCESS;
}

bool Test1DEvaluation()
{
  GaussianModel p(1);
  Eigen::VectorXd mean = Eigen::VectorXd::Zero(1);
  p.SetMean(mean);
  Eigen::VectorXd variance = Eigen::MatrixXd::Identity(1,1);
  p.SetVariance(variance);
  Eigen::VectorXd v(1);
  v(0) = .3;
  double eval = p.Evaluate(v);
  //std::cout << "1D: " << eval << std::endl;
  double epsilon = 1e-6;

  if(fabs(eval -  0.381388) > epsilon)
  {
      return false;
  }

  return true;
}

bool Test1D_EM()
{
  int dimensionality = 1;

  // Generate some data
  Eigen::MatrixXd data = GenerateData(40, dimensionality);

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
  expectationMaximization.SetRandom(false);
  expectationMaximization.SetMixtureModel(mixtureModel);
  expectationMaximization.SetMaxIterations(3);
  expectationMaximization.Compute();

  // This is where we got the test output
  MixtureModel finalModel = expectationMaximization.GetMixtureModel();

  for(unsigned int i = 0; i < finalModel.GetNumberOfModels(); ++i)
  {
    std::cout << "Model " << i << ":" << std::endl;
    finalModel.GetModel(i)->Print();
  }

  std::vector<double> testMeans;
  testMeans.push_back(4.90469);
  testMeans.push_back(0.00576815);

  std::vector<double> testVariances;
  testVariances.push_back(0.830928);
  testVariances.push_back(0.862264);

  double epsilon = 1e-4;

  for(unsigned int i = 0; i < finalModel.GetNumberOfModels(); ++i)
  {
    if(fabs(testMeans[i] - finalModel.GetModel(i)->GetMean()(0)) > epsilon)
    {
        return false;
    }

    if(fabs(testVariances[i] - finalModel.GetModel(i)->GetVariance()(0,0)) > epsilon)
    {
        return false;
    }
  }

  return true;
}


Eigen::MatrixXd GenerateData(const unsigned int numberOfSamplesPerGroup, const unsigned int dimensionality)
{
  std::default_random_engine generator(0); // Seed with 0, so results are repeatable

  // Mean 0, Standard deviation 1
  Eigen::MatrixXd data1(dimensionality, numberOfSamplesPerGroup);
  std::normal_distribution<double> normalDistribution1(0,sqrt(1));
  for(unsigned int i = 0; i < numberOfSamplesPerGroup; i++)
  {
    Eigen::VectorXd v(dimensionality);
    v(0) = normalDistribution1(generator);
    data1.col(i) = v;
  }

  // Mean 5, Standard deviation 1
  Eigen::MatrixXd data2(dimensionality, numberOfSamplesPerGroup);
  std::normal_distribution<double> normalDistribution2(5,sqrt(1));
  for(unsigned int i = 0; i < numberOfSamplesPerGroup; i++)
  {
    Eigen::VectorXd v = Eigen::VectorXd::Random(dimensionality);
    v(0) = normalDistribution2(generator);
    data2.col(i) = v;
  }

  // Concatentate the matrices horiztonally
  Eigen::MatrixXd data(dimensionality, data1.cols() + data2.cols());
  data << data1, data2;

  return data;
}
