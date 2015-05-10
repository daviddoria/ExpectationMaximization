#ifndef __ExpectationMaximization_h
#define __ExpectationMaximization_h

#include "Gaussian1D.h"

#include <Eigen/StdVector> // Required (http://eigen.tuxfamily.org/dox-devel/TopicStlContainers.html)

#include "KMeansClustering/KMeansClustering.h"

class ExpectationMaximization
{
public:
  /* Constructor. */
  ExpectationMaximization();

  unsigned int GetNumberOfModels() {return this->Models.size();}

  void SetData(Eigen::MatrixXd data){this->Data = data;}
  void SetModels(std::vector<Model*> models) {this->Models = models;}
  void AddModel(Model* model) {this->Models.push_back(model);}
  unsigned int NumberOfDataPoints() {return this->Data.cols();}

  std::vector<Model*> GetModels() {return this->Models;}
  Model* GetModel(int i) {return this->Models[i];}

  void SetInitializationTechniqueToRandom(){this->InitializationTechnique = RANDOM;}
  void SetInitializationTechniqueToKMeans(){this->InitializationTechnique = KMEANS;}
  
  double WeightedEvaluate(Eigen::VectorXd x);

  void Compute();

  void SetMinChange(const float minChange);

protected:

  int MaxIterations;
  float MinChange;

  std::vector<Model*> Models;
  Eigen::MatrixXd Data;

  Eigen::MatrixXd Responsibilities;
  
  void RandomlyInitializeModels();
  void KMeansInitializeModels();
  
  enum InitializationEnum {RANDOM, KMEANS};
  
  int InitializationTechnique;

};

#endif
