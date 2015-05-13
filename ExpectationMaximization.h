/*
Copyright (C) 2015 David Doria, daviddoria@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __ExpectationMaximization_h
#define __ExpectationMaximization_h

#include "Model.h"
#include "MixtureModel.h"

#include <Eigen/StdVector> // Required (http://eigen.tuxfamily.org/dox-devel/TopicStlContainers.html)

#include "KMeansClustering/KMeansClustering.h"

class ExpectationMaximization
{
public:
  /** Constructor. */
  ExpectationMaximization();

  /** The data to cluster. */
  void SetData(const Eigen::MatrixXd& data){this->Data = data;}

  /** The mixture model whose parameters we want to find. */
  void SetMixtureModel(const MixtureModel& mixtureModel) {this->Model = mixtureModel;}

  /** Get the number of data points. */
  unsigned int NumberOfDataPoints() const {return this->Data.cols();}

  /** Get the mixture model. */
  MixtureModel GetMixtureModel() const {return this->Model;}

  /** Set the initialization method to random. */
  void SetInitializationTechniqueToRandom(){this->InitializationTechnique = RANDOM;}

  /** Set the initialization method to use the KMeans to initialize the models. */
  void SetInitializationTechniqueToKMeans(){this->InitializationTechnique = KMEANS;}

  /** Perform the parameter optimization. */
  void Compute();

  /** Set a stopping criteria that determines how much the models have changed iteration to iteration. */
  void SetMinChange(const float minChange);

  /** Set a stopping criteria that tracks the number of iterations that have been run. */
  void SetMaxIterations(const unsigned int maxIterations);

  /** Specify if the "random" numbers should be random (set to 'false' for repeatability for testing). */
  void SetRandom(const bool r);

  /** Get a seed for a random number generator based on the Random flag. */
  long int GetSeed() const;

protected:

  /** A stopping condition (the number of iterations to perform). */
  int MaxIterations;

  /** A stopping condition (if the error changes more than this between iterations, continue running). */
  float MinChange;

  /** The mixture model of which we are estimating the parameters. */
  MixtureModel Model;

  /** The data to cluster. */
  Eigen::MatrixXd Data;

  /** The current weighted likelihoods. */
  Eigen::MatrixXd Responsibilities;
  
  /** Randomly initialize the models. */
  void RandomlyInitializeModels();

  /** Use KMeans to initialize the models. */
  void KMeansInitializeModels();
  
  /** The possible initialization methods. */
  enum InitializationEnum {RANDOM, KMEANS};
  
  /** The technique to use to initialize the models. */
  int InitializationTechnique = RANDOM;

  /** The flag that determine if we use a fixed or random seed. */
  bool Random = true;
};

#endif
