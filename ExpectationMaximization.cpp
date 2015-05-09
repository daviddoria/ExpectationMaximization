// Custom
#include "ExpectationMaximization.h"

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
    Model* model = new Gaussian1D; // This can be any derived type, we just need access to the member data
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
        double val = this->Models[model]->WeightedEvaluate(this->Data[point]);
        if(IsNaN(val))
          {
          std::cout << "val of " << this->Data[point] << " is nan for model " << model << "!" << std::endl;
          //OutputModelInfo(this->Models);
          exit(-1);
          }
        normalization += val;
        }
      //std::cout << "normalization: " << normalization << std::endl;
      for(unsigned int model = 0; model < this->GetNumberOfModels(); model++)
        {
        double resp = this->Models[model]->WeightedEvaluate(this->Data[point]);
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

      Eigen::VectorXd newMean(this->Models[0]->GetDimensionality(), 0);
      for(unsigned int i = 0; i < this->NumberOfDataPoints(); i++)
        {
        newMean += this->Responsibilities(i, model) * this->Data[i];
        }
      newMean /= sumResponsibilities;
      this->Models[model]->SetMean(newMean);

      Eigen::MatrixXd newVariance(this->Models[0]->GetDimensionality(), this->Models[0]->GetDimensionality(), 0);
      for(unsigned int i = 0; i < this->NumberOfDataPoints(); i++)
        {
        //SIGMA = (1/N) sum_{n=1}^N (x_n - u)(x_n - u)^T
        newVariance += this->Responsibilities(i,model) * (this->Data[i] - newMean) * (this->Data[i] - newMean).transpose();
        }
      newVariance /= sumResponsibilities;
      this->Models[model]->SetVariance(newVariance);
      }

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

    //std::cout << "Models at iteration " << iter << std::endl;
    //OutputModelInfo(this->Models);
    //PlotModels(this->Models, range);

    // Check if we have converged
    done = true; // assume we have converged until we see otherwise
    for(unsigned int model = 0; model < this->GetNumberOfModels(); model++)
      {
      //double meanDiff = fabs(this->Models[model]->GetMean()(0) - lastModels[model]->GetMean()(0));
      double meanDiff = (this->Models[model]->GetMean() - lastModels[model]->GetMean()).norm();

      //double varianceDiff = fabs(this->Models[model]->GetVariance()(0,0) - lastModels[model]->GetVariance()(0,0));

      //if((meanDiff > this->MinChange) || (varianceDiff > this->MinChange) )
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

  for(unsigned int model = 0; model < this->Models.size(); model++)
    {
    // Mean
    Eigen::VectorXd mean(dim);
    for(unsigned int d = 0; d < dim; d++)
      {
      mean(d) = vtkMath::Random(-5, 5);
      }
    this->Models[model]->SetMean(mean);

    // Variance
    std::vector<double> variance(dim);
    for(unsigned int d = 0; d < dim; d++)
      {
      variance[d] = vtkMath::Random(0, 3);
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

  VectorOfPoints clusterCenters = kmeans.GetClusterCenters();

  unsigned int dim = this->Models[0]->GetDimensionality();

  for(unsigned int model = 0; model < this->Models.size(); model++)
    {
    // Set initial EM means from KMeans cluster centers
    this->Models[model]->SetMean(clusterCenters[model]);

    // Set initial EM variances from KMeans cluster bounds
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    kmeans.GetPointsWithLabel(model, points);
    double bounds[6];
    points->GetBounds(bounds);

    std::vector<double> variance(dim);
    for(unsigned int d = 0; d < dim; d++)
      {
      variance[d] = bounds[2*d + 1] - bounds[2*d]; // max - min
      }

    this->Models[model]->SetDiagonalCovariance(variance);

    // The ratio of points in the cluster to total points
    this->Models[model]->SetMixingCoefficient(static_cast<float>(kmeans->GetIndicesWithLabel(model).size())/static_cast<float>(this->Data.size()));
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

double ExpectationMaximization::WeightedEvaluate(Eigen::VectorXd x)
{
  double sum = 0;
  for(unsigned int i = 0; i < this->Models.size(); i++)
    {
    sum += this->Models[i]->WeightedEvaluate(x);
    }
  return sum;
}
