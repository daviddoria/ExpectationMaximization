// STL
#include <iostream>
#include <random>
#include <chrono>

// Custom
#include "ExpectationMaximization.h"
#include "GaussianModel.h"

// Submodules
#include "KMeansClustering/KMeansClustering.h"

// good slides: http://www.slideshare.net/petitegeek/expectation-maximization-and-gaussian-mixture-models

bool IsNaN(const double a);

ExpectationMaximization::ExpectationMaximization()
{
  // Defaults
  this->MaxIterations = 10;
  this->MinChange = 1e-3;
  this->InitializationTechnique = KMEANS;
}

void ExpectationMaximization::Compute()
{
  std::cout << "ExpectationMaximization::Compute" << std::endl;
  this->Responsibilities.resize(this->Data.size(),this->Model.GetNumberOfModels());

  if(this->InitializationTechnique == RANDOM)
  {
    RandomlyInitializeModels();
  }
  else if(this->InitializationTechnique == KMEANS)
  {
    KMeansInitializeModels();
  }

  bool done = false;
  int iter = 0;

  // An object used to save the old model
  MixtureModel previousModel;

  do
    {
    std::cout << "EM iteration " << iter << "..." << std::endl;
    // Save a snapshot of the current models
    previousModel.CopyFrom(this->Model);

    // E step - evaluate responsibilities
    for(unsigned int pointId = 0; pointId < this->NumberOfDataPoints(); pointId++)
    {
      double normalization = this->Model.WeightedEvaluate(this->Data.col(pointId));

      //std::cout << "normalization: " << normalization << std::endl;
      for(unsigned int modelId = 0; modelId < this->Model.GetNumberOfModels(); modelId++)
      {
        double resp = this->Model.GetModel(modelId)->WeightedEvaluate(this->Data.col(pointId));
        this->Responsibilities(pointId, modelId) = resp / normalization;
      }
    }

    // M step - estimate parameters based on the new responsibilities
    for(unsigned int modelId = 0; modelId < this->Model.GetNumberOfModels(); modelId++)
    {
      double sumResponsibilities = 0;
      for(unsigned int i = 0; i < this->NumberOfDataPoints(); i++)
      {
        sumResponsibilities += this->Responsibilities(i, modelId);
      }
      //std::cout << "sumResponsibilities : " << sumResponsibilities << std::endl;;

      Eigen::VectorXd newMean = Eigen::VectorXd::Zero(this->Model.GetDimensionality());
      for(unsigned int i = 0; i < this->NumberOfDataPoints(); i++)
      {
        newMean += this->Responsibilities(i, modelId) * this->Data.col(i);
      }
      newMean /= sumResponsibilities;
      this->Model.GetModel(modelId)->SetMean(newMean);

      Eigen::MatrixXd newVariance(this->Model.GetDimensionality(), this->Model.GetDimensionality());
      newVariance.setZero(this->Model.GetDimensionality(), this->Model.GetDimensionality());
      for(unsigned int i = 0; i < this->NumberOfDataPoints(); i++)
      {
        //SIGMA = (1/N) sum_{n=1}^N (x_n - u)(x_n - u)^T
        newVariance += this->Responsibilities(i,modelId) * (this->Data.col(i) - newMean) * (this->Data.col(i) - newMean).transpose();
      }
      newVariance /= sumResponsibilities;
      this->Model.GetModel(modelId)->SetVariance(newVariance);
    } // end M step loop

    // Update mixing coefficients
    for(unsigned int modelId = 0; modelId < this->Model.GetNumberOfModels(); modelId++)
    {
      double sumResponsibilities = 0;
      for(unsigned int i = 0; i < this->NumberOfDataPoints(); i++)
      {
        sumResponsibilities += this->Responsibilities(i, modelId);
      }
      this->Model.GetModel(modelId)->SetMixingCoefficient(sumResponsibilities/static_cast<double>(this->NumberOfDataPoints()));
    }

    // Check if we have converged
    done = true; // assume we have converged until we see otherwise
    for(unsigned int modelId = 0; modelId < this->Model.GetNumberOfModels(); modelId++)
    {
      double meanDiff = (this->Model.GetModel(modelId)->GetMean() - previousModel.GetModel(modelId)->GetMean()).norm();

      // If any of the models have not converged, we need to continue
      if((meanDiff > this->MinChange))
      {
        done = false;
      }
    }
  iter++;
  //OutputModelInfo(this->Models);
  }while(!done && (iter < this->MaxIterations) );

  std::cout << "EM performed " << iter << " iterations." << std::endl;
}

void ExpectationMaximization::RandomlyInitializeModels()
{
  unsigned int dim = this->Model.GetDimensionality();

  std::default_random_engine generator;

  for(unsigned int modelId = 0; modelId < this->Model.GetNumberOfModels(); modelId++)
  {
    // Mean
    Eigen::VectorXd mean(dim);
    std::uniform_real_distribution<double> meanDistribution(-5.0,5.0);
    for(unsigned int d = 0; d < dim; d++)
    {
      mean(d) = meanDistribution(generator);
    }
    this->Model.GetModel(modelId)->SetMean(mean);

    // Variance
    Eigen::VectorXd variance(dim);

    std::uniform_real_distribution<double> varianceDistribution(0.0,3.0);
    for(unsigned int d = 0; d < dim; d++)
    {
      variance(d) = varianceDistribution(generator);
    }

    Eigen::MatrixXd covariance = Eigen::MatrixXd::Identity(dim, dim);
    covariance.diagonal() = variance;

    this->Model.GetModel(modelId)->SetVariance(covariance);

    // Mixing coefficient
    this->Model.GetModel(modelId)->SetMixingCoefficient(1./this->Model.GetNumberOfModels());
  }
}

void ExpectationMaximization::KMeansInitializeModels()
{
  // Initialize with KMeans
  KMeansClustering kmeans;
  kmeans.SetK(this->Model.GetNumberOfModels());
  kmeans.SetPoints(this->Data);
  kmeans.Cluster();

  Eigen::MatrixXd clusterCenters = kmeans.GetClusterCenters();

  unsigned int dim = this->Model.GetDimensionality();

  for(unsigned int modelId = 0; modelId < this->Model.GetNumberOfModels(); modelId++)
  {
    // Set initial EM means from KMeans cluster centers
    this->Model.GetModel(modelId)->SetMean(clusterCenters.col(modelId));

    // Set initial EM variances from KMeans cluster bounds
    // See 'Partial Reductions': http://eigen.tuxfamily.org/dox/group__TutorialReductionsVisitorsBroadcasting.html
    Eigen::MatrixXd points = kmeans.GetPointsWithLabel(modelId);
    Eigen::VectorXd minValues = points.rowwise().minCoeff();
    Eigen::VectorXd maxValues = points.rowwise().maxCoeff();

    Eigen::VectorXd center(dim);
    for(unsigned int d = 0; d < dim; d++)
    {
      center(d) = (maxValues(d) - minValues(d))/2.0f;
    }

    Eigen::MatrixXd covariance = Eigen::MatrixXd::Identity(dim, dim);
    covariance.diagonal() = center;

    this->Model.GetModel(modelId)->SetVariance(covariance);

    // The ratio of points in the cluster to total points
    this->Model.GetModel(modelId)->SetMixingCoefficient(static_cast<float>(kmeans.GetIndicesWithLabel(modelId).size())/static_cast<float>(this->Data.cols()));
  }
}

bool IsNaN(const double a)
{
  if(a!=a)
  {
    return true;
  }
  return false;
}

void ExpectationMaximization::SetMinChange(const float minChange)
{
  this->MinChange = minChange;
}

void ExpectationMaximization::SetMaxIterations(const unsigned int maxIterations)
{
  this->MaxIterations = maxIterations;
}

void ExpectationMaximization::SetRandom(const bool r)
{
  this->Random = r;
}

long int ExpectationMaximization::GetSeed() const
{
  if(this->Random)
  {
    auto now = std::chrono::system_clock::now();
    return now.time_since_epoch().count();
  }

  return 0;
}
