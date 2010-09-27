#ifndef __vtkExpectationMaximization_h
#define __vtkExpectationMaximization_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkSmartPointer.h"
#include "vtkDenseArray.h"
#include "vtkMath.h"

#include "Gaussian1D.h"

class vtkExpectationMaximization : public vtkPolyDataAlgorithm
{
public:
  static vtkExpectationMaximization *New();
  vtkTypeMacro(vtkExpectationMaximization,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  int GetNumberOfModels() {return this->Models.size();}

  vtkSetMacro(MaxIterations, int);

  vtkSetMacro(MinChange, double);

  void SetData(std::vector<vnl_vector<double> > data){this->Data = data;}
  void SetModels(std::vector<Model*> models) {this->Models = models;}
  void AddModel(Model* model) {this->Models.push_back(model);}
  int NumberOfDataPoints() {return this->Data.size();}

  std::vector<Model*> GetModels() {return this->Models;}
  Model* GetModel(int i) {return this->Models[i];}

protected:
  vtkExpectationMaximization();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  int MaxIterations;
  int MinChange;

  std::vector<Model*> Models;
  std::vector<vnl_vector<double> > Data;

  vtkSmartPointer<vtkDenseArray<double> > Responsibilities;

private:
  vtkExpectationMaximization(const vtkExpectationMaximization&);  // Not implemented.
  void operator=(const vtkExpectationMaximization&);  // Not implemented.
};

#endif
