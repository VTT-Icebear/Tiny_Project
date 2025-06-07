#ifndef MATRIX_H
#define MATRIX_H

#include "Vector.h"

class Matrix {
private:
    int mNumRows;     // Number of rows
    int mNumCols;     // Number of columns
    double** mData;   // Pointer to array of pointers (2D array)

public:
    // Constructors and destructor
    Matrix(int numRows, int numCols);     // Constructor with dimensions
    Matrix(const Matrix& other);          // Copy constructor
    ~Matrix();                           // Destructor

    // Assignment operator
    Matrix& operator=(const Matrix& other);

    // Access operators
    double& operator()(int i, int j);     // 1-based indexing
    const double& operator()(int i, int j) const;

    // Unary operators
    Matrix operator+() const;             // Unary plus
    Matrix operator-() const;             // Unary minus

    // Binary operators
    Matrix operator+(const Matrix& other) const;  // Matrix addition
    Matrix operator-(const Matrix& other) const;  // Matrix subtraction
    Matrix operator*(const Matrix& other) const;  // Matrix multiplication
    Matrix operator*(double scalar) const;        // Scalar multiplication
    Vector operator*(const Vector& vec) const;    // Matrix-vector multiplication
    friend Matrix operator*(double scalar, const Matrix& mat);  // Scalar multiplication (left)

    // Compound assignment operators
    Matrix& operator+=(const Matrix& other);
    Matrix& operator-=(const Matrix& other);
    Matrix& operator*=(double scalar);

    // Utility functions
    int GetNumRows() const { return mNumRows; }
    int GetNumCols() const { return mNumCols; }
    bool IsSquare() const { return mNumRows == mNumCols; }
    bool IsSymmetric() const;             // Check if matrix is symmetric
    double Determinant() const;           // Calculate determinant
    Matrix Inverse() const;               // Calculate inverse
    Matrix PseudoInverse() const;         // Calculate Moore-Penrose pseudo-inverse
    Matrix Transpose() const;             // Calculate transpose

    // Friend functions for I/O
    friend ostream& operator<<(ostream& os, const Matrix& mat);
    friend istream& operator>>(istream& is, Matrix& mat);
};

// Constructor
Matrix::Matrix(int numRows, int numCols) : mNumRows(numRows), mNumCols(numCols) {
    assert(numRows > 0 && numCols > 0);
    mData = new double*[numRows];
    for (int i = 0; i < numRows; i++) {
        mData[i] = new double[numCols];
        for (int j = 0; j < numCols; j++) {
            mData[i][j] = 0.0;
        }
    }
}

// Copy constructor
Matrix::Matrix(const Matrix& other) : mNumRows(other.mNumRows), mNumCols(other.mNumCols) {
    mData = new double*[mNumRows];
    for (int i = 0; i < mNumRows; i++) {
        mData[i] = new double[mNumCols];
        for (int j = 0; j < mNumCols; j++) {
            mData[i][j] = other.mData[i][j];
        }
    }
}

// Destructor
Matrix::~Matrix() {
    for (int i = 0; i < mNumRows; i++) {
        delete[] mData[i];
    }
    delete[] mData;
}

// Assignment operator
Matrix& Matrix::operator=(const Matrix& other) {
    if (this != &other) {
        // Delete old data
        for (int i = 0; i < mNumRows; i++) {
            delete[] mData[i];
        }
        delete[] mData;

        // Copy new data
        mNumRows = other.mNumRows;
        mNumCols = other.mNumCols;
        mData = new double*[mNumRows];
        for (int i = 0; i < mNumRows; i++) {
            mData[i] = new double[mNumCols];
            for (int j = 0; j < mNumCols; j++) {
                mData[i][j] = other.mData[i][j];
            }
        }
    }
    return *this;
}

// Access operator
double& Matrix::operator()(int i, int j) {
    assert(i >= 1 && i <= mNumRows && j >= 1 && j <= mNumCols);
    return mData[i-1][j-1];
}

const double& Matrix::operator()(int i, int j) const {
    assert(i >= 1 && i <= mNumRows && j >= 1 && j <= mNumCols);
    return mData[i-1][j-1];
}

// Unary operators
Matrix Matrix::operator+() const {
    return *this;
}

Matrix Matrix::operator-() const {
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            result.mData[i][j] = -mData[i][j];
        }
    }
    return result;
}

// Binary operators
Matrix Matrix::operator+(const Matrix& other) const {
    assert(mNumRows == other.mNumRows && mNumCols == other.mNumCols);
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            result.mData[i][j] = mData[i][j] + other.mData[i][j];
        }
    }
    return result;
}

Matrix Matrix::operator-(const Matrix& other) const {
    assert(mNumRows == other.mNumRows && mNumCols == other.mNumCols);
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            result.mData[i][j] = mData[i][j] - other.mData[i][j];
        }
    }
    return result;
}

Matrix Matrix::operator*(const Matrix& other) const {
    assert(mNumCols == other.mNumRows);
    Matrix result(mNumRows, other.mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < other.mNumCols; j++) {
            double sum = 0.0;
            for (int k = 0; k < mNumCols; k++) {
                sum += mData[i][k] * other.mData[k][j];
            }
            result.mData[i][j] = sum;
        }
    }
    return result;
}

Matrix Matrix::operator*(double scalar) const {
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            result.mData[i][j] = mData[i][j] * scalar;
        }
    }
    return result;
}

Vector Matrix::operator*(const Vector& vec) const {
    assert(mNumCols == vec.GetSize());
    Vector result(mNumRows);
    for (int i = 0; i < mNumRows; i++) {
        double sum = 0.0;
        for (int j = 0; j < mNumCols; j++) {
            sum += mData[i][j] * vec[j];
        }
        result[i] = sum;
    }
    return result;
}

Matrix operator*(double scalar, const Matrix& mat) {
    return mat * scalar;
}

// Compound assignment operators
Matrix& Matrix::operator+=(const Matrix& other) {
    assert(mNumRows == other.mNumRows && mNumCols == other.mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            mData[i][j] += other.mData[i][j];
        }
    }
    return *this;
}

Matrix& Matrix::operator-=(const Matrix& other) {
    assert(mNumRows == other.mNumRows && mNumCols == other.mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            mData[i][j] -= other.mData[i][j];
        }
    }
    return *this;
}

Matrix& Matrix::operator*=(double scalar) {
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            mData[i][j] *= scalar;
        }
    }
    return *this;
}

// Utility functions
bool Matrix::IsSymmetric() const {
    if (!IsSquare()) return false;
    for (int i = 0; i < mNumRows; i++) {
        for (int j = i + 1; j < mNumCols; j++) {
            if (mData[i][j] != mData[j][i]) return false;
        }
    }
    return true;
}

Matrix Matrix::Transpose() const {
    Matrix result(mNumCols, mNumRows);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            result.mData[j][i] = mData[i][j];
        }
    }
    return result;
}

double Matrix::Determinant() const {
    assert(IsSquare());
    // For simplicity, we'll implement a basic determinant calculation
    // This is not the most efficient method for large matrices
    if (mNumRows == 1) return mData[0][0];
    if (mNumRows == 2) {
        return mData[0][0] * mData[1][1] - mData[0][1] * mData[1][0];
    }

    double det = 0.0;
    for (int j = 0; j < mNumCols; j++) {
        // Create submatrix
        Matrix submatrix(mNumRows - 1, mNumCols - 1);
        for (int i = 1; i < mNumRows; i++) {
            for (int k = 0; k < mNumCols; k++) {
                if (k < j) submatrix.mData[i-1][k] = mData[i][k];
                else if (k > j) submatrix.mData[i-1][k-1] = mData[i][k];
            }
        }
        det += (j % 2 == 0 ? 1 : -1) * mData[0][j] * submatrix.Determinant();
    }
    return det;
}

Matrix Matrix::Inverse() const {
    assert(IsSquare());
    double det = Determinant();
    assert(abs(det) > 1e-10); // Check if matrix is singular

    Matrix result(mNumRows, mNumCols);
    // For simplicity, we'll implement a basic inverse calculation
    // This is not the most efficient method for large matrices
    if (mNumRows == 1) {
        result.mData[0][0] = 1.0 / mData[0][0];
        return result;
    }
    if (mNumRows == 2) {
        result.mData[0][0] = mData[1][1] / det;
        result.mData[0][1] = -mData[0][1] / det;
        result.mData[1][0] = -mData[1][0] / det;
        result.mData[1][1] = mData[0][0] / det;
        return result;
    }

    // For larger matrices, we'll use the adjugate method
    // This is not the most efficient method, but it's straightforward
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            // Create submatrix
            Matrix submatrix(mNumRows - 1, mNumCols - 1);
            for (int k = 0; k < mNumRows; k++) {
                for (int l = 0; l < mNumCols; l++) {
                    if (k < i && l < j) submatrix.mData[k][l] = mData[k][l];
                    else if (k < i && l > j) submatrix.mData[k][l-1] = mData[k][l];
                    else if (k > i && l < j) submatrix.mData[k-1][l] = mData[k][l];
                    else if (k > i && l > j) submatrix.mData[k-1][l-1] = mData[k][l];
                }
            }
            result.mData[j][i] = ((i + j) % 2 == 0 ? 1 : -1) * submatrix.Determinant() / det;
        }
    }
    return result;
}

Matrix Matrix::PseudoInverse() const {
    // Moore-Penrose pseudo-inverse using SVD
    // For simplicity, we'll use a basic implementation
    // In practice, you would want to use a more robust method
    Matrix A = *this;
    Matrix AT = A.Transpose();
    Matrix ATA = AT * A;
    Matrix AAT = A * AT;

    // Check if ATA is invertible
    if (abs(ATA.Determinant()) > 1e-10) {
        return AT * (A * AT).Inverse();
    }
    // Check if AAT is invertible
    else if (abs(AAT.Determinant()) > 1e-10) {
        return (AT * A).Inverse() * AT;
    }
    else {
        // If neither is invertible, use the formula: A+ = (A^T * A)^(-1) * A^T
        // with regularization
        double lambda = 1e-10; // Small regularization parameter
        Matrix I(mNumCols, mNumCols);
        for (int i = 0; i < mNumCols; i++) {
            I.mData[i][i] = 1.0;
        }
        return (AT * A + lambda * I).Inverse() * AT;
    }
}

// I/O operators
ostream& operator<<(ostream& os, const Matrix& mat) {
    for (int i = 0; i < mat.mNumRows; i++) {
        os << "[";
        for (int j = 0; j < mat.mNumCols; j++) {
            os << mat.mData[i][j];
            if (j < mat.mNumCols - 1) os << ", ";
        }
        os << "]" << (i < mat.mNumRows - 1 ? "\n" : "");
    }
    return os;
}

istream& operator>>(istream& is, Matrix& mat) {
    for (int i = 0; i < mat.mNumRows; i++) {
        for (int j = 0; j < mat.mNumCols; j++) {
            is >> mat.mData[i][j];
        }
    }
    return is;
}

#endif // MATRIX_H