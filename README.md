# Tiny_Project
# README: Tiny Project

## Overview
The Tiny Project is a C++-based implementation divided into two parts:

- **Part A**: A linear algebra library featuring `Vector`, `Matrix`, and `LinearSystem` classes. It supports vector and matrix operations, solves square linear systems (Ax = b) using Gaussian elimination, and handles symmetric positive definite systems with the conjugate gradient method. It also addresses under-determined and over-determined systems using the Moore-Penrose pseudo-inverse with Tikhonov regularization.
- **Part B**: A linear regression model predicting CPU performance (PRP) using the UCI Computer Hardware dataset (https://archive.ics.uci.edu/dataset/29/computer+hardware). The model trains on six features (MYCT, MMIN, MMAX, CACH, CHMIN, CHMAX) and evaluates performance with RMSE.

This project was completed as of Saturday, June 07, 2025, 11:06 PM +07.

## Setup Instructions
1. Ensure a C++ compiler (e.g., g++) is installed.
2. Clone or download the project repository/files.
3. Place the `machine.data` file (from the UCI dataset) in the project directory.
4. Compile the source files:
   - For `Vector.cpp`: `g++ -std=c++11 Vector.cpp -o vector_test`
   - For `Matrix.cpp`: `g++ -std=c++11 Matrix.cpp -o matrix_test`
   - For `LinearSystem.cpp`: `g++ -std=c++11 LinearSystem.cpp -o linear_system_test`
   - For `Source.cpp` (Part B): `g++ -std=c++11 Source.cpp -o linear_regression`
5. Run the compiled executables to test the implementations.

## Usage
- **Part A**:
  - Run `vector_test` to test vector operations (e.g., addition, norm).
  - Run `matrix_test` to test matrix operations (e.g., multiplication, inverse).
  - Run `linear_system_test` to solve linear systems using Gaussian elimination or conjugate gradient methods.
  - Input data interactively as prompted (e.g., vector or matrix elements).
- **Part B**:
  - Run `linear_regression` to train the model and output beta coefficients and RMSE.
  - Ensure `machine.data` is present; the program automatically loads and processes it.

## Notes
- The library assumes valid input sizes and may throw exceptions for singular matrices or incompatible dimensions.
- The linear regression model in Part B uses an 80/20 train-test split; results may vary due to random shuffling.
- For best results, consider normalizing features or adding regularization to improve RMSE.

## Acknowledgments
- Dataset source: UCI Machine Learning Repository (https://archive.ics.uci.edu/dataset/29/computer+hardware).
