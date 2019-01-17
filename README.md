# dualnum

An example implementation of dual number-based tangent linear and adjoint models.
This repository contains the code for a simple nutrient-phytoplankton-zooplankton (NPZ) model, including its nonlinear code and dual number-based tangent linear and adjoint code. It also includes a code to test the consistency between models.

## Getting Started

All that is needed is a Python 3 installation with [numpy](https://www.numpy.org/) and (optionally) [matplotlib](http://www.matplotlib.org/).

## Running the tests

To see if it is all working, try running the NPZ nonlinear and tangent linear model (requires matplotlib):
```
python3 NPZmodel.py
```
or perform a series of dot-product tests:
```
python3 dotproduct_test.py
```

## Authors

* Jann Paul Mattern [jpmattern](https://github.com/jpmattern)


## License

This project is licensed under the Apache License - see the [LICENSE.md](LICENSE.md) file for details.
