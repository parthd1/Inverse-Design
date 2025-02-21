# Inverse-Design

Repo for Inverse Design Code.

## Table of Contents

- [Description](#description)
- [Installation](#installation)
- [Usage](#usage)
- [Functions](#functions)
- [Examples](#examples)
- [License](#license)

## Description

This repository contains MATLAB code for implementing inverse design algorithms, primarily focusing on optimizing geometric parameters to achieve a power splitter.
## Installation

1. Clone the repository:
   ```sh
   git clone https://github.com/parthd1/Inverse-Design.git
   ```
2. Navigate to the repository directory:
   ```sh
   cd Inverse-Design
   ```
## Functions

### powerSpliter_DBS.m

This MATLAB script is used for simulating a power splitter design. Below are the key components and their descriptions:

- **Dependencies**:
  - Adds paths to COMSOL libraries.
  - Imports necessary COMSOL classes: `com.comsol.model.*` and `com.comsol.model.util.*`.

- **COMSOL Server**:
  - Attempts to start a COMSOL server connection using `mphstart(2036);`.

- **Model Initialization**:
  - Removes any existing model named 'Model' using `ModelUtil.remove('Model');`.

- **Parameters**:
  - Defines geometry parameters (`geom_para`) and coupler parameters (`coupler_geom`).
  - Stores parameters in a structure `param`.

- **Split Ratio**:
  - Sets the target splitting ratio (`splitRatio`) and a hyperparameter (`hyperparam`).
  - Outputs the chosen parameters and target ratio.

- **Model Selection**:
  - Selects between different power splitter models based on the value of `modelUsed`.
  - Models:
    - Model 1: No PML boundary.
    - Model 2: PML with full cover.
    - Model 3: No PML boundary with coupler.
    - Model 4: PML with full cover with coupler.

- **Simulation**:
  - Calls the `DBS` function with the defined parameters: `[totalLoss, bestImages] = DBS(param, 0.001, splitRatio, hyperparam, modelUsed);`.
  - Saves the generated data to `./genData/data0.mat`.
  - Prints a completion message.

### random_iter.m

The `random_iter.m` file is a MATLAB script that performs an iterative optimization process on a binary image using COMSOL modeling. Below is an overview of its main components and functionality:

- **Function Definition**: The main function `random_iter` takes several parameters including `binaryImage`, `geom_params`, `target`, `hyperparam`, `modelUsed`, and `iterNum`.
- **Initialization**:
  - Imports necessary COMSOL model utilities.
  - Initializes storage for data such as `binaryImages`, `power1`, `power2`, `sMatrix`, and `losses`.
  - Creates a COMSOL model and updates it with initial parameters.
- **Simulation Loop**:
  - Iteratively updates the binary image by modifying pixels randomly.
  - Runs simulations to compute the loss and other metrics (`p1`, `p2`, `s11`).
  - Stores the results of each iteration.
- **Data Saving**:
  - Saves the best image and associated data to a `.mat` file at the end of the iterations.
- **Helper Function**:
  - `next_iter` function updates the binary image by flipping random pixels.

Overall, this script is used to optimize a binary image iteratively by running simulations and updating the image to minimize a loss function.
Sure, here is a formatted explanation of the `modelResult.m` file for your README:



### modelResult.m

The `modelResult.m` file in this repository defines a MATLAB function to evaluate and extract results from a given model. Below is a summary of its functionality:

#### Function Definition

```matlab
function [s11, p1, p2] = modelResult(model)
```

This function takes a `model` as input and returns three outputs: `s11`, `p1`, and `p2`.

#### Creating Evaluation Group

```matlab
eg1 = model.result.evaluationGroup.create('eg1', 'EvaluationGroup');
```

An evaluation group `eg1` is created to hold different evaluation features.

#### Adding Evaluation Features

Several evaluation features are created and configured, such as `s11`, `totalPower`, `outPower1`, and `outPower2`. For instance, the feature `s11` is set to evaluate `abs(solid.Smatrix11)`.

#### Running the Evaluation

```matlab
eg1.run;
data = mphtable(model, 'eg1');
s11 = data.data(2);
p1 = data.data(3);
p2 = data.data(4);
```

The evaluation group `eg1` is executed, and the results are extracted from the evaluation group and assigned to the output variables `s11`, `p1`, and `p2`.

#### Helper Function

```matlab
function pmlBlabels = getOutLabels(model)
```

This helper function `getOutLabels` calculates boundary labels for PML (Perfectly Matched Layer) regions based on model parameters and returns them.


### DBS.m
### Function Signature

```matlab
function [totalLoss, bestImages] = DBS(geom_params, sigma, target, hyperparam, modelUsed)
```
- **Inputs:**
  - `geom_params`: Geometrical parameters for the design.
  - `sigma`: Threshold value for the loss.
  - `target`: Target specifications for the design.
  - `hyperparam`: Hyperparameters for the model.
  - `modelUsed`: The model being used for the design process.
- **Outputs:**
  - `totalLoss`: Array of loss values from each iteration.
  - `bestImages`: Array of the best images from each iteration.

### Main Code

1. **Initialize random seed:**
   ```matlab
   rng("shuffle");
   ```

2. **Initialize binary image:**
   ```matlab
   binaryImage = get_init_image('init_best.mat');
   fprintf('The total holes for random initilization is : %d\n', sum(binaryImage, 'all'))
   ```

3. **Initialize variables:**
   ```matlab
   iterNum = 1;
   totalLoss = [];
   bestImages = [];
   ```

4. **First iteration:**
   ```matlab
   [loss, bestImage] = random_iter(binaryImage, geom_params, target, hyperparam, modelUsed, iterNum);
   totalLoss = [totalLoss, loss];
   bestImages = [bestImages, bestImage];
   ```

5. **Iterate until loss is below sigma:**
   ```matlab
   while loss > sigma
       iterNum = iterNum + 1;
       [loss, bestImage] = random_iter(bestImage, geom_params, target, hyperparam, modelUsed, iterNum);
       totalLoss = [totalLoss, loss];
       bestImages = cat(3, bestImages, bestImage);
   end
   ```

### Helper Function

- **Load initial image:**
  ```matlab
  function init_image = get_init_image(name)
      init_image = load(name).best;
  end
  ```

