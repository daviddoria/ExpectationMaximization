#include <vtkSmartPointer.h>
#include <vtkMath.h>

#include "vtkExpectationMaximization.h"
#include "GaussianND.h"
#include "Plotting.h"

std::vector<vnl_vector<double> > GenerateData(int numberOfSamplesPerGroup, int dimensionality);

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



std::vector<vnl_vector<double> > GenerateData(int numberOfSamplesPerGroup, int dimensionality)
{
  // This demo creates 3 groups of 40 points each. It then combines them
  // without maintaining any information of which groups the points came from.

  vtkMath::RandomSeed(time(NULL));

  // Mean 0, Standard deviation 1
  std::vector<vnl_vector<double> > data1;
  for(unsigned int i = 0; i < numberOfSamplesPerGroup; i++)
    {
    vnl_vector<double> v(dimensionality);
    for(unsigned int d = 0; d < dimensionality; d++)
      {
      v(d) = vtkMath::Gaussian(0.0,1.0);
      }
    data1.push_back(v);
    }

  // Mean 5, Standard deviation 1
  std::vector<vnl_vector<double> > data2;
  for(unsigned int i = 0; i < numberOfSamplesPerGroup; i++)
    {
    vnl_vector<double> v(dimensionality);
    for(unsigned int d = 0; d < dimensionality; d++)
      {
      v(d) = vtkMath::Gaussian(5.0,1.0);
      }
    data2.push_back(v);
    }

  // Mean -5, Standard deviation 1.5
  std::vector<vnl_vector<double> > data3;
  for(unsigned int i = 0; i < numberOfSamplesPerGroup; i++)
    {
    vnl_vector<double> v(dimensionality);
    for(unsigned int d = 0; d < dimensionality; d++)
      {
      v(d) = vtkMath::Gaussian(-5.0,1.0);
      }
    data3.push_back(v);
    }

  std::vector<vnl_vector<double> > data;
  data.insert(data.end(), data1.begin(), data1.end());
  data.insert(data.end(), data2.begin(), data2.end());
  data.insert(data.end(), data3.begin(), data3.end());

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
  
  vnl_vector<double> mean(1);
  mean(0) = 0;
  p.SetMean(mean);
  
  vnl_matrix<double> variance(1,1);
  variance(0,0) = 2.0;
  
  p.SetVariance(variance);
  
  vnl_vector<double> x(1);
  x(0) = 0.3;
  std::cout << "ND: " << p.Evaluate(x) << std::endl;
}

void Test1D()
{
   
  int dimensionality = 1;
  
  // Generate some data
  std::vector<vnl_vector<double> > data = GenerateData(40, dimensionality);

  // Initialize the model
  std::vector<Model*> models(3);

  for(int i = 0; i < models.size(); i++)
    {
    Model* model = new Gaussian1D;
    model->SetDimensionality(dimensionality);
    
    dynamic_cast<Gaussian1D*>(model)->SetMean(vtkMath::Random(-5, 5));
    
    dynamic_cast<Gaussian1D*>(model)->SetVariance(vtkMath::Random(0, 3));
    model->SetMixingCoefficient(1./models.size());
    models[i] = model;
    }

  // Display initial model
  std::cout << "Randomly initialized model:" << std::endl;
  //OutputModelInfo(models);
  double range[2] = {-15, 15};
  PlotModels1D(models, range);

  vtkSmartPointer<vtkExpectationMaximization> expectationMaximization =
    vtkSmartPointer<vtkExpectationMaximization>::New();
  expectationMaximization->SetData(data);
  expectationMaximization->SetModels(models);
  expectationMaximization->SetMinChange(1e-5);

  expectationMaximization->Update();

  std::cout << "Final models:" << std::endl;
  //OutputModelInfo(expectationMaximization->GetModels());
  PlotModels1D(expectationMaximization->GetModels(), range);

}

void TestND()
{
  
  int dimensionality = 2;
  
  // Generate some data
  std::vector<vnl_vector<double> > data = GenerateData(40, dimensionality);

  // Initialize the model
  std::vector<Model*> models(3);

  for(int i = 0; i < models.size(); i++)
    {
    Model* model = new GaussianND;
    model->SetDimensionality(dimensionality);
    model->Init();
    vnl_vector<double> mean(2);
    mean(0) = vtkMath::Random(-5, 5);
    mean(1) = vtkMath::Random(-5, 5);
    model->SetMean(mean);
    vnl_matrix<double> variance(2,2);
    variance(0,0) = vtkMath::Random(0, 3);
    variance(1,1) = vtkMath::Random(0, 3);
    variance(1,0) = 0;
    variance(0,1) = 0;
    model->SetVariance(variance);
    model->SetMixingCoefficient(1./models.size());
    models[i] = model;
    }

  // Display initial model
  std::cout << "Randomly initialized model:" << std::endl;
  //OutputModelInfo(models);
  PlotModels2D(models, data);

  vtkSmartPointer<vtkExpectationMaximization> expectationMaximization =
    vtkSmartPointer<vtkExpectationMaximization>::New();
  expectationMaximization->SetData(data);
  expectationMaximization->SetModels(models);
  expectationMaximization->SetMinChange(1e-5);

  expectationMaximization->Update();

  std::cout << "Final models:" << std::endl;
  //OutputModelInfo(expectationMaximization->GetModels());
  PlotModels2D(expectationMaximization->GetModels(), data);
 
}