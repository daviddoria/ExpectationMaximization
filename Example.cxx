#include <vtkSmartPointer.h>
#include <vtkMath.h>

#include "vtkExpectationMaximization.h"
#include "Model.h"

std::vector<double> GenerateData(int numberOfSamplesPerGroup);

int main(int, char*[])
{
  // Generate some data
  std::vector<double> data = GenerateData(40);

  // Initialize the model
  std::vector<Model> models(3);

  for(int i = 0; i < models.size(); i++)
    {
    Model model;
    model.SetMean(vtkMath::Random(-5, 5));
    model.SetVariance(vtkMath::Random(0, 3));
    model.SetMixingCoefficient(1./models.size());
    models[i] = model;
    }

  // Display initial model
  std::cout << "Randomly initialized model:" << std::endl;
  OutputModelInfo(models);
  double range[2] = {-15, 15};
  PlotModels(models, range);

  vtkSmartPointer<vtkExpectationMaximization> expectationMaximization =
    vtkSmartPointer<vtkExpectationMaximization>::New();
  expectationMaximization->SetData(data);
  expectationMaximization->SetModels(models);
  expectationMaximization->SetMinChange(1e-2);

  expectationMaximization->Update();

  std::cout << "Final models:" << std::endl;
  OutputModelInfo(expectationMaximization->GetModels());
  PlotModels(expectationMaximization->GetModels(), range);

  return EXIT_SUCCESS;
}



std::vector<double> GenerateData(int numberOfSamplesPerGroup)
{
  // This demo creates 3 groups of 40 points each. It then combines them
  // without maintaining any information of which groups the points came from.

  vtkMath::RandomSeed(time(NULL));

  // Mean 0, Standard deviation 2
  std::vector<double> data1;
  for(unsigned int i = 0; i < numberOfSamplesPerGroup; i++)
    {
    double a = vtkMath::Gaussian(0.0,1.0);
    data1.push_back(a);
    }

  // Mean 3, Standard deviation 1
  std::vector<double> data2;
  for(unsigned int i = 0; i < numberOfSamplesPerGroup; i++)
    {
    double a = vtkMath::Gaussian(5.0,1.0);
    data2.push_back(a);
    }

  // Mean -4, Standard deviation 1.5
  std::vector<double> data3;
  for(unsigned int i = 0; i < numberOfSamplesPerGroup; i++)
    {
    double a = vtkMath::Gaussian(-5.0,1.0);
    data3.push_back(a);
    }

  std::vector<double> data;
  data.insert(data.end(), data1.begin(), data1.end());
  data.insert(data.end(), data2.begin(), data2.end());
  data.insert(data.end(), data3.begin(), data3.end());

  return data;
}