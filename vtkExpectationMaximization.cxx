// Custom
#include "vtkExpectationMaximization.h"
#include "Plotting.h"

// Submodules
#include "KMeansClustering/vtkKMeansClustering.h"

// VTK
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkMath.h"

// VNL
#include <vnl/vnl_vector.h>

vtkStandardNewMacro(vtkExpectationMaximization);

// good slides: http://www.slideshare.net/petitegeek/expectation-maximization-and-gaussian-mixture-models

bool IsNaN(const double a);

vtkExpectationMaximization::vtkExpectationMaximization()
{
  // Filter
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(0);

  // Defaults
  this->MaxIterations = 10;
  this->MinChange = 1e-3;
  this->InitializationTechnique = KMEANS;

  // Initializations
  this->Responsibilities = vtkSmartPointer<vtkDenseArray<double> >::New();
}

int vtkExpectationMaximization::RequestData(
                              vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector)
{
  std::cout << "vtkExpectationMaximization::RequestData" << std::endl;
  this->Responsibilities->Resize(this->Data.size(),this->GetNumberOfModels());

  if(this->InitializationTechnique == RANDOM)
    {
    RandomlyInitializeModels();
    }
  else if(this->InitializationTechnique == KMEANS)
    {
    KMeansInitializeModels();
    }

  OutputModelInfo(this->Models);

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
	  OutputModelInfo(this->Models);
	  exit(-1);
	  }
	normalization += val;
	}
      //std::cout << "normalization: " << normalization << std::endl;
      for(unsigned int model = 0; model < this->GetNumberOfModels(); model++)
	{
	double resp = this->Models[model]->WeightedEvaluate(this->Data[point]);
	this->Responsibilities->SetValue(point, model, resp / normalization);
	}
      }

    // M step - estimate parameters based on the new responsibilities
    for(unsigned int model = 0; model < this->GetNumberOfModels(); model++)
      {
      double sumResponsibilities = 0;
      for(unsigned int i = 0; i < this->NumberOfDataPoints(); i++)
        {
        sumResponsibilities += this->Responsibilities->GetValue(i, model);
        }
      //std::cout << "sumResponsibilities : " << sumResponsibilities << std::endl;;

      vnl_vector<double> newMean(this->Models[0]->GetDimensionality(), 0);
      for(unsigned int i = 0; i < this->NumberOfDataPoints(); i++)
        {
        newMean += this->Responsibilities->GetValue(i, model) * this->Data[i];
        }
      newMean /= sumResponsibilities;
      this->Models[model]->SetMean(newMean);

      vnl_matrix<double> newVariance(this->Models[0]->GetDimensionality(), this->Models[0]->GetDimensionality(), 0);
      for(unsigned int i = 0; i < this->NumberOfDataPoints(); i++)
        {
	//SIGMA = (1/N) sum_{n=1}^N (x_n - u)(x_n - u)^T
        newVariance += this->Responsibilities->GetValue(i,model) * outer_product(this->Data[i] - newMean, this->Data[i] - newMean);
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
        sumResponsibilities += this->Responsibilities->GetValue(i, model);
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
      double meanDiff = (this->Models[model]->GetMean() - lastModels[model]->GetMean()).magnitude();

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

  return 1;
}

void vtkExpectationMaximization::RandomlyInitializeModels()
{
  unsigned int dim = this->Models[0]->GetDimensionality();

  for(unsigned int model = 0; model < this->Models.size(); model++)
    {
    // Mean
    vnl_vector<double> mean(dim);
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

void vtkExpectationMaximization::KMeansInitializeModels()
{
  // Convert to a format suitable for kmeans
  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();
  for(unsigned int i = 0; i < this->Data.size(); i++)
    {
    double p[3] = {0,0,0};
    for(unsigned int d = 0; d < this->Data[i].size(); d++) // Assumes data is 3d or lower
      {
      p[d] = this->Data[i](d);
      }
    points->InsertNextPoint(p);
    }

  vtkSmartPointer<vtkPolyData> polydata =
    vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);

  // Initialize with KMeans
  vtkSmartPointer<vtkKMeansClustering> kmeans =
    vtkSmartPointer<vtkKMeansClustering>::New();
  kmeans->SetK(this->GetNumberOfModels());
  std::cout << "Set K to " << this->GetNumberOfModels() << std::endl;
  kmeans->SetInputData(polydata);
  kmeans->Update();

  vtkPolyData* clusterCenters = kmeans->GetOutput(1);

  unsigned int dim = this->Models[0]->GetDimensionality();

  for(unsigned int model = 0; model < this->Models.size(); model++)
    {
    // Set initial EM means from KMeans cluster centers
    vnl_vector<double> mean(dim);
    double p[3];
    clusterCenters->GetPoint(model, p);
    for(unsigned int d = 0; d < dim; d++)
      {
      mean(d) = p[d];
      }
    this->Models[model]->SetMean(mean);

    // Set initial EM variances from KMeans cluster bounds

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    kmeans->GetPointsWithLabel(model, points);
    double bounds[6];
    points->GetBounds(bounds);

    std::vector<double> variance(dim);
    for(unsigned int d = 0; d < dim; d++)
      {
      variance[d] = bounds[2*d + 1] - bounds[2*d];
      }

    this->Models[model]->SetDiagonalCovariance(variance);

    // The ratio of points in the cluster to total points
    this->Models[model]->SetMixingCoefficient(static_cast<float>(kmeans->GetIndicesWithLabel(model).size())/static_cast<float>(this->Data.size()));
    }
}

void vtkExpectationMaximization::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "NumberOfModels: " << this->GetNumberOfModels() << std::endl
                << "MaxIterations: " << this->MaxIterations << std::endl
                << "MinChange: " << this->MinChange << std::endl;

}

bool IsNaN(const double a)
{
  if(a!=a)
  {
    return true;
  }
  return false;
}

double vtkExpectationMaximization::WeightedEvaluate(vnl_vector<double> x)
{
  double sum = 0;
  for(unsigned int i = 0; i < this->Models.size(); i++)
    {
    sum += this->Models[i]->WeightedEvaluate(x);
    }
  return sum;
}