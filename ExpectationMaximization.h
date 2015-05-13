#ifndef __ExpectationMaximization_h
#define __ExpectationMaximization_h

#include "Model.h"
#include "MixtureModel.h"

#include <Eigen/StdVector> // Required (http://eigen.tuxfamily.org/dox-devel/TopicStlContainers.html)

#include "KMeansClustering/KMeansClustering.h"

class ExpectationMaximization
{
public:
  /* Constructor. */
  ExpectationMaximization();

  void SetData(const Eigen::MatrixXd& data){this->Data = data;}

  void SetMixtureModel(const MixtureModel& mixtureModel) {this->Model = mixtureModel;}

  unsigned int NumberOfDataPoints() const {return this->Data.cols();}

  MixtureModel GetMixtureModel() const {return this->Model;}

  void SetInitializationTechniqueToRandom(){this->InitializationTechnique = RANDOM;}
  void SetInitializationTechniqueToKMeans(){this->InitializationTechnique = KMEANS;}

  void Compute();

  /** Set a stopping criteria that determines how much the models have changed iteration to iteration. */
  void SetMinChange(const float minChange);

  /** Set a stopping criteria that tracks the number of iterations that have been run. */
  void SetMaxIterations(const unsigned int maxIterations);

  /**
   * Specify random or not for repeatability for testing.
   */
  void SetRandom(const bool r);

  /** Get a seed for a random number generator based on the Random flag. */
  long int GetSeed() const;

protected:

  int MaxIterations;
  float MinChange;

  MixtureModel Model;
  Eigen::MatrixXd Data;

  Eigen::MatrixXd Responsibilities;
  
  void RandomlyInitializeModels();
  void KMeansInitializeModels();
  
  enum InitializationEnum {RANDOM, KMEANS};
  
  int InitializationTechnique = RANDOM;

  bool Random = true;
};

#endif
