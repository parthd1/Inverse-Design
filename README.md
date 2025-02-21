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

### `DBS.m`
- **Description**: Optimizes an initial binary image to minimize a loss function.
- **Usage**:
  ```matlab
  [totalLoss, bestImages] = DBS(geom_params, sigma, target, hyperparam, modelUsed);
  ```

### `modelResult.m`
- **Description**: Evaluates the result of the model and extracts specific performance metrics.
- **Usage**:
  ```matlab
  [s11, p1, p2] = modelResult(model);
  ```

## Examples

### Running `DBS.m`
```matlab
geom_params = [...]; % Define geometric parameters
sigma = 0.01;
target = [...]; % Define target values
hyperparam = [...]; % Define hyperparameters
modelUsed = ...; % Specify the model used

[totalLoss, bestImages] = DBS(geom_params, sigma, target, hyperparam, modelUsed);
```

### Running `modelResult.m`
```matlab
model = ...; % Specify the model
[s11, p1, p2] = modelResult(model);
```
