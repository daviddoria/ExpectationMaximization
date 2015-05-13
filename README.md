This is a c++ implementation of the Expectation-Maximization algorithm.


Example usage:

    // Initialize the models
    std::vector<Model*> models(2);

    for(unsigned int i = 0; i < models.size(); i++)
    {
      Model* model = new GaussianModel(dimensionality);
      models[i] = model;
    }

    MixtureModel mixtureModel;
    mixtureModel.SetModels(models);
    ExpectationMaximization expectationMaximization;
    expectationMaximization.SetData(data);
    expectationMaximization.SetModels(models);
    expectationMaximization.SetMinChange(1e-5);
    expectationMaximization.SetMaxIterations(100);
    expectationMaximization.Compute();
