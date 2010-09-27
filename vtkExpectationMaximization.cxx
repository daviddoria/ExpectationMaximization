#include "vtkExpectationMaximization.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkMath.h"

vtkStandardNewMacro(vtkExpectationMaximization);

vtkExpectationMaximization::vtkExpectationMaximization()
{
  // Filter
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(0);

  // Defaults
  this->MaxIterations = 1000;
  this->MinChange = 1e-3;

  // Initializations
  this->Responsibilities = vtkSmartPointer<vtkDenseArray<double> >::New();
}

int vtkExpectationMaximization::RequestData(
                              vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector)
{
  this->Responsibilities->Resize(this->Data.size(),this->GetNumberOfModels());

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
        for(int model = 0; model < this->GetNumberOfModels(); model++)
          {
          normalization += this->Models[model]->WeightedEvaluate(this->Data[point]);
          }
       for(int model = 0; model < this->GetNumberOfModels(); model++)
          {
          double resp = this->Models[model]->WeightedEvaluate(this->Data[point]);
          this->Responsibilities->SetValue(point, model, resp / normalization);
          }
      }

    // M step - estimate parameters based on the new responsibilities
    for(int model = 0; model < this->GetNumberOfModels(); model++)
      {
      double sumResponsibilities = 0;
      for(int i = 0; i < this->NumberOfDataPoints(); i++)
        {
        sumResponsibilities += this->Responsibilities->GetValue(i, model);
        }

      vnl_vector<double> newMean(1, 0);
      for(int i = 0; i < this->NumberOfDataPoints(); i++)
        {
        newMean(0) += this->Responsibilities->GetValue(i, model) * this->Data[i](0);
        }
      newMean(0) /= sumResponsibilities;
      this->Models[model]->SetMean(newMean);

      vnl_matrix<double> newVariance(1,1,0);
      for(int i = 0; i < this->NumberOfDataPoints(); i++)
        {
        newVariance(0,0) += this->Responsibilities->GetValue(i,model) * pow(this->Data[i](0) - newMean(0), 2);
        }
      newVariance(0,0) /= sumResponsibilities;
      this->Models[model]->SetVariance(newVariance);
      }

    // Update mixing coefficients
    for(int model = 0; model < this->GetNumberOfModels(); model++)
      {
      double sumResponsibilities = 0;
      for(int i = 0; i < this->NumberOfDataPoints(); i++)
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
      double meanDiff = fabs(this->Models[model]->GetMean()(0) - lastModels[model]->GetMean()(0));
      //std::cout << "meanDiff: " << meanDiff << std::endl;
      double varianceDiff = fabs(this->Models[model]->GetVariance()(0,0) - lastModels[model]->GetVariance()(0,0));
      //std::cout << "varianceDiff: " << varianceDiff << std::endl;
      if((meanDiff > this->MinChange) || (varianceDiff > this->MinChange) )
        {
        done = false;
        }
      }
  iter++;
  }while(!done && (iter < this->MaxIterations) );

  std::cout << "Performed " << iter << " iterations." << std::endl;
  
  return 1;
}

void vtkExpectationMaximization::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "NumberOfModels: " << this->GetNumberOfModels() << std::endl
                << "MaxIterations: " << this->MaxIterations << std::endl
                << "MinChange: " << this->MinChange << std::endl;

}
