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
  this->Responsibilities.resize(this->Data.size(),this->GetNumberOfModels());

  if(this->InitializationTechnique == RANDOM)
  {
    RandomlyInitializeModels();
  }
  else if(this->InitializationTechnique == KMEANS)
  {
    KMeansInitializeModels();
  }

  //OutputModelInfo(this->Models);

  bool done = false;
  int iter = 0;

  // We need to actually copy the data
  std::vector<Model*> lastModels(this->Models.size());
  for(unsigned int i = 0; i < this->Models.size(); i++)
  {
    Model* model = new GaussianModel(1); // This can be any derived type, we just need access to the member data
    lastModels[i] = model;
  }

  do
    {
    std::cout << "EM iteration " << iter << std::endl;
    // Save a snapshot of the current models
    //std::copy(this->Models.begin(), this->Models.end(), lastModels.begin());
    for(unsigned int i = 0; i < this->Models.size(); i++)
    {
      lastModels[i]->SetVariance(this->Models[i]->GetVariance());
      lastModels[i]->SetMean(this->Models[i]->GetMean());
    }

    // E step - evaluate responsibilities
    for(unsigned int point = 0; point < this->NumberOfDataPoints(); point++)
    {
      double normalization = 0.0;
      for(unsigned int model = 0; model < this->GetNumberOfModels(); model++)
      {
        double val = this->Models[model]->WeightedEvaluate(this->Data.col(point));
        if(IsNaN(val))
        {
          std::cout << "val of " << this->Data.col(point) << " is nan for model " << model << "!" << std::endl;
          //OutputModelInfo(this->Models);
          exit(-1);
        }
        normalization += val;
      }
      //std::cout << "normalization: " << normalization << std::endl;
      for(unsigned int model = 0; model < this->GetNumberOfModels(); model++)
      {
        double resp = this->Models[model]->WeightedEvaluate(this->Data.col(point));
        this->Responsibilities(point, model) = resp / normalization;
      }
    }

    // M step - estimate parameters based on the new responsibilities
    for(unsigned int model = 0; model < this->GetNumberOfModels(); model++)
    {
      double sumResponsibilities = 0;
      for(unsigned int i = 0; i < this->NumberOfDataPoints(); i++)
      {
        sumResponsibilities += this->Responsibilities(i, model);
      }
      //std::cout << "sumResponsibilities : " << sumResponsibilities << std::endl;;

      Eigen::VectorXd newMean = Eigen::VectorXd::Zero(this->Models[0]->GetDimensionality());
      for(unsigned int i = 0; i < this->NumberOfDataPoints(); i++)
      {
        newMean += this->Responsibilities(i, model) * this->Data.col(i);
      }
      newMean /= sumResponsibilities;
      this->Models[model]->SetMean(newMean);

      Eigen::MatrixXd newVariance(this->Models[0]->GetDimensionality(), this->Models[0]->GetDimensionality());
      newVariance.setZero(this->Models[0]->GetDimensionality(), this->Models[0]->GetDimensionality());
      for(unsigned int i = 0; i < this->NumberOfDataPoints(); i++)
      {
        //SIGMA = (1/N) sum_{n=1}^N (x_n - u)(x_n - u)^T
        newVariance += this->Responsibilities(i,model) * (this->Data.col(i) - newMean) * (this->Data.col(i) - newMean).transpose();
      }
      newVariance /= sumResponsibilities;
      this->Models[model]->SetVariance(newVariance);
    } // end M step loop

    // Update mixing coefficients
    for(unsigned int model = 0; model < this->GetNumberOfModels(); model++)
    {
      double sumResponsibilities = 0;
      for(unsigned int i = 0; i < this->NumberOfDataPoints(); i++)
      {
        sumResponsibilities += this->Responsibilities(i, model);
      }
      this->Models[model]->SetMixingCoefficient(sumResponsibilities/static_cast<double>(this->NumberOfDataPoints()));
    }

    // Check if we have converged
    done = true; // assume we have converged until we see otherwise
    for(unsigned int model = 0; model < this->GetNumberOfModels(); model++)
    {
      double meanDiff = (this->Models[model]->GetMean() - lastModels[model]->GetMean()).norm();

      // If any of the models have not converged, we need to continue
      std::cout << "model " << model << " meanDiff: " << meanDiff << std::endl;
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
  unsigned int dim = this->Models[0]->GetDimensionality();

  std::default_random_engine generator;

  for(unsigned int model = 0; model < this->Models.size(); model++)
  {
    // Mean
    Eigen::VectorXd mean(dim);
    std::uniform_real_distribution<double> meanDistribution(-5.0,5.0);
    for(unsigned int d = 0; d < dim; d++)
    {
      mean(d) = meanDistribution(generator);
    }
    this->Models[model]->SetMean(mean);

    // Variance
    Eigen::VectorXd variance(dim);

    std::uniform_real_distribution<double> varianceDistribution(0.0,3.0);
    for(unsigned int d = 0; d < dim; d++)
    {
      variance(d) = varianceDistribution(generator);
    }

    this->Models[model]->SetDiagonalCovariance(variance);

    // Mixing coefficient
    this->Models[model]->SetMixingCoefficient(1./this->GetNumberOfModels());
  }
}

void ExpectationMaximization::KMeansInitializeModels()
{
  // Initialize with KMeans
  KMeansClustering kmeans;
  kmeans.SetK(this->GetNumberOfModels());
  kmeans.SetPoints(this->Data);
  kmeans.Cluster();

  Eigen::MatrixXd clusterCenters = kmeans.GetClusterCenters();

  unsigned int dim = this->Models[0]->GetDimensionality();

  for(unsigned int modelId = 0; modelId < this->Models.size(); modelId++)
  {
    // Set initial EM means from KMeans cluster centers
    this->Models[modelId]->SetMean(clusterCenters.col(modelId));

    // Set initial EM variances from KMeans cluster bounds
    // See 'Partial Reductions': http://eigen.tuxfamily.org/dox/group__TutorialReductionsVisitorsBroadcasting.html
    Eigen::MatrixXd points = kmeans.GetPointsWithLabel(modelId);
    Eigen::VectorXd minValues = points.rowwise().minCoeff();
    Eigen::VectorXd maxValues = points.rowwise().maxCoeff();

    Eigen::VectorXd range(dim);
    for(unsigned int d = 0; d < dim; d++)
    {
      range(d) = maxValues(d) - minValues(d);
    }

    this->Models[modelId]->SetDiagonalCovariance(range);

    // The ratio of points in the cluster to total points
    this->Models[modelId]->SetMixingCoefficient(static_cast<float>(kmeans.GetIndicesWithLabel(modelId).size())/static_cast<float>(this->Data.cols()));
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

double ExpectationMaximization::WeightedEvaluate(const Eigen::VectorXd& x) const
{
  double sum = 0;
  for(unsigned int i = 0; i < this->Models.size(); i++)
  {
    sum += this->Models[i]->WeightedEvaluate(x);
  }
  return sum;
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
