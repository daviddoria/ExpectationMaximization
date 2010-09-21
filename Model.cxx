#include "Model.h"

#include <vtkDenseArray.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkChartXY.h>
#include <vtkPlot.h>
#include <vtkTable.h>
#include <vtkFloatArray.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>

#include <sstream>

void OutputModelInfo(std::vector<Model> models)
{
  for(unsigned int i = 0; i < models.size(); i++)
    {
    std::cout << "Model " << i << " : Mean = " << models[i].GetMean()
              << " Variance = " << models[i].GetVariance()
              << " Mixing coefficient = " << models[i].GetMixingCoefficient() << std::endl;
    }
}

void PlotModels(std::vector<Model> models, double range[2])
{
  // Create a table with some points in it
  vtkSmartPointer<vtkTable> table =
    vtkSmartPointer<vtkTable>::New();

  vtkSmartPointer<vtkFloatArray> axis =
    vtkSmartPointer<vtkFloatArray>::New();
  axis->SetName("X Axis");
  table->AddColumn(axis);

  for(unsigned int i = 0; i < models.size(); i++)
    {
    std::stringstream ss;
    ss << "Model " << i;
    vtkSmartPointer<vtkFloatArray> modelArray =
      vtkSmartPointer<vtkFloatArray>::New();
    modelArray->SetName(ss.str().c_str());
    table->AddColumn(modelArray);
    }

  int numPoints = 100;

  float step = (range[1] - range[0]) / (numPoints-1);
  table->SetNumberOfRows(numPoints);

  // Setup axis
  for (int i = 0; i < numPoints; ++i)
    {
    double xval = range[0] + i * step;
    table->SetValue(i, 0, xval);
    }

  for(int model = 0; model < models.size(); model++)
    {
    for (int i = 0; i < numPoints; ++i)
      {
      double xval = range[0] + i * step;
      table->SetValue(i, model+1, models[model].Evaluate(xval));
      }
    }

  // Set up the view
  vtkSmartPointer<vtkContextView> view =
    vtkSmartPointer<vtkContextView>::New();
  view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);

  // Add multiple line plots, setting the colors etc
  vtkSmartPointer<vtkChartXY> chart =
    vtkSmartPointer<vtkChartXY>::New();
  view->GetScene()->AddItem(chart);
  vtkPlot *line = chart->AddPlot(vtkChart::LINE);
  line->SetInput(table, 0, 1);
  line->SetColor(0, 255, 0, 255);
  line->SetWidth(1.0);
  line = chart->AddPlot(vtkChart::LINE);
  line->SetInput(table, 0, 2);
  line->SetColor(255, 0, 0, 255);
  line->SetWidth(5.0);
  line = chart->AddPlot(vtkChart::LINE);
  line->SetInput(table, 0, 3);
  line->SetColor(0, 0, 255, 255);
  line->SetWidth(5.0);

  // Set up an interactor and start
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(view->GetRenderWindow());
  renderWindowInteractor->Initialize();
  renderWindowInteractor->Start();

}

void PlotPoints(std::vector<double> points, vtkSmartPointer<vtkDenseArray<double> > responsibilities)
{
  // Create a table with some points in it
  vtkSmartPointer<vtkTable> table =
    vtkSmartPointer<vtkTable>::New();

  vtkSmartPointer<vtkFloatArray> axis =
    vtkSmartPointer<vtkFloatArray>::New();
  axis->SetName("X Axis");
  table->AddColumn(axis);

  vtkSmartPointer<vtkFloatArray> zeros =
    vtkSmartPointer<vtkFloatArray>::New();
  zeros->SetName("zeros");
  table->AddColumn(zeros);

  table->SetNumberOfRows(points.size());

  for (int i = 0; i < points.size(); ++i)
    {
    table->SetValue(i, 0, points[i]);
    table->SetValue(i, 1, 0);
    }

  // Set up the view
  vtkSmartPointer<vtkContextView> view =
    vtkSmartPointer<vtkContextView>::New();
  view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);

  // Add multiple line plots, setting the colors etc
  vtkSmartPointer<vtkChartXY> chart =
    vtkSmartPointer<vtkChartXY>::New();
  view->GetScene()->AddItem(chart);
  vtkPlot *line = chart->AddPlot(vtkChart::LINE);
  line->SetInput(table, 0, 1);
  line->SetColor(0, 255, 0, 255);
  line->SetWidth(1.0);
  line = chart->AddPlot(vtkChart::LINE);
  line->SetInput(table, 0, 2);
  line->SetColor(255, 0, 0, 255);
  line->SetWidth(5.0);
  line = chart->AddPlot(vtkChart::LINE);
  line->SetInput(table, 0, 3);
  line->SetColor(0, 0, 255, 255);
  line->SetWidth(5.0);

  // Set up an interactor and start
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(view->GetRenderWindow());
  renderWindowInteractor->Initialize();
  renderWindowInteractor->Start();

}
