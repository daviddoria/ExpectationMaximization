#ifndef PLOTTING_H
#define PLOTTING_H

#include "Model.h"
#include <vector>
#include <vtkSmartPointer.h>
#include <vtkDenseArray.h>

class vtkPolyData;

void PlotModels1D(std::vector<Model*> models, double range[2]);
void PlotModels2D(std::vector<Model*> models, std::vector<vnl_vector<double> > data);
void PlotPoints(std::vector<double> points, vtkSmartPointer<vtkDenseArray<double> > responsibilities);
void CreateSurface(Model* model, double xrange[2], double yrange[2], int divisions, vtkPolyData* surface);

#endif