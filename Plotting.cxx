#include "Plotting.h"

#include <vtkDenseArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkDelaunay2D.h>
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

void PlotModels1D(std::vector<Model*> models, double range[2])
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
      vnl_vector<double> v(1);
      v(0) = range[0] + i * step;
      table->SetValue(i, model+1, models[model]->Evaluate(v));
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

void CreateSurface(Model* model, double xrange[2], double yrange[2], int divisions, vtkPolyData* surface)
{
  std::cout << "Creating surface on (" << xrange[0] << " , " << xrange[1] << ") , (" << yrange[0] << " , " << yrange[1] << " )" << std::endl;
  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();
 
  double xstep = (xrange[1] - xrange[0]) / static_cast<double>(divisions);
  double ystep = (yrange[1] - yrange[0]) / static_cast<double>(divisions);
  
  std::cout << "xstep: " << xstep << " ystep: " << ystep << std::endl;
  
  for(int x = 0; x < divisions; x++)
    {
    for(int y = 0; y < divisions; y++)
      {
      vnl_vector<double> p(2);
      p(0) = xrange[0] + x*xstep;
      p(1) = yrange[0] + y*ystep;
      points->InsertNextPoint(p(0), p(1), model->Evaluate(p));
      }
    }
 
  // Add the grid points to a polydata object
  vtkSmartPointer<vtkPolyData> polydata =
    vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);

  vtkSmartPointer<vtkDelaunay2D> delaunay =
    vtkSmartPointer<vtkDelaunay2D>::New();
  delaunay->SetInput(polydata);
  delaunay->Update();
 
  surface = delaunay->GetOutput();
  
}

void PlotModels2D(std::vector<Model*> models, std::vector<vnl_vector<double> > data)
{
  
  // Create a renderer, render window, and interactor
  vtkSmartPointer<vtkRenderer> renderer = 
    vtkSmartPointer<vtkRenderer>::New();
    
  vtkSmartPointer<vtkRenderWindow> renderWindow = 
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  
  // Set up an interactor and start
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  
  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();
  for(unsigned int i = 0; i < data.size(); i++)
    {
    points->InsertNextPoint(data[i](0), data[i](1), 0);
    }
    
  vtkSmartPointer<vtkPolyData> pointsPolyData = 
    vtkSmartPointer<vtkPolyData>::New();
  pointsPolyData->SetPoints(points);
  
  vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter = 
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  glyphFilter->SetInputConnection(pointsPolyData->GetProducerPort());
  glyphFilter->Update();
  
  // Create a mapper and actor
  vtkSmartPointer<vtkPolyDataMapper> pointsMapper = 
    vtkSmartPointer<vtkPolyDataMapper>::New();
  pointsMapper->SetInputConnection(glyphFilter->GetOutputPort());

  vtkSmartPointer<vtkActor> pointsActor = 
    vtkSmartPointer<vtkActor>::New();
  pointsActor->SetMapper(pointsMapper);
  renderer->AddActor(pointsActor);
  
  for(unsigned int i = 0; i < models.size(); i++)
    {
    vtkSmartPointer<vtkPolyData> surface = 
      vtkSmartPointer<vtkPolyData>::New();
      
    std::cout << "Model " << i << " dimensionality: " << models[i]->GetDimensionality() << std::endl;
    double xrange[2] = {models[i]->GetMean()(0) - 2.0, models[i]->GetMean()(0) + 2.0};
    double yrange[2] = {models[i]->GetMean()(1) - 2.0, models[i]->GetMean()(1) + 2.0};
    CreateSurface(models[i], xrange, yrange, 20, surface);

    // Create a mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> mapper = 
      vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(surface->GetProducerPort());
  
    vtkSmartPointer<vtkActor> actor = 
      vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    
    renderer->AddActor(actor);
    }

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
