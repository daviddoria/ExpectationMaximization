#ifndef __vtkExpectationMaximization_h
#define __vtkExpectationMaximization_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkSmartPointer.h"
#include "vtkDenseArray.h"
#include "vtkMath.h"

#include "Model.h"

class vtkExpectationMaximization : public vtkPolyDataAlgorithm
{
public:
  static vtkExpectationMaximization *New();
  vtkTypeMacro(vtkExpectationMaximization,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  int GetNumberOfModels() {return this->Models.size();}

  vtkSetMacro(MaxIterations, int);

  vtkSetMacro(MinChange, double);

  void SetData(std::vector<double> data){this->Data = data;}
  void SetModels(std::vector<Model> models) {this->Models = models;}
  int NumberOfDataPoints() {return this->Data.size();}

  std::vector<Model> GetModels() {return this->Models;}

protected:
  vtkExpectationMaximization();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  int MaxIterations;
  int MinChange;

  std::vector<Model> Models;
  std::vector<double> Data;

  vtkSmartPointer<vtkDenseArray<double> > Responsibilities;

private:
  vtkExpectationMaximization(const vtkExpectationMaximization&);  // Not implemented.
  void operator=(const vtkExpectationMaximization&);  // Not implemented.
};

#endif
