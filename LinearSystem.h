#ifndef LINEAR_SYSTEM_H
#define LINEAR_SYSTEM_H

#include "Matrix.h"

class LinearSystem {
protected:
    int mSize;          // Size of the linear system
    Matrix* mpA;        // Pointer to the coefficient matrix
    Vector* mpb;        // Pointer to the right-hand side vector

public:
    // Constructor
    LinearSystem(const Matrix& A, const Vector& b);
    
    // Prevent default construction and copying
    LinearSystem() = delete;
    LinearSystem(const LinearSystem& other) = delete;
    LinearSystem& operator=(const LinearSystem& other) = delete;
    
    // Destructor
    virtual ~LinearSystem();

    // Virtual solve method to be overridden by derived classes
    virtual Vector Solve() const;

    // Utility functions
    int GetSize() const { return mSize; }
    const Matrix& GetMatrix() const { return *mpA; }
    const Vector& GetVector() const { return *mpb; }
};

// Derived class for positive definite symmetric linear systems
class PosSymLinSystem : public LinearSystem {
public:
    // Constructor
    PosSymLinSystem(const Matrix& A, const Vector& b);
    
    // Override solve method to use conjugate gradient method
    Vector Solve() const override;

private:
    // Helper functions for conjugate gradient method
    double ComputeAlpha(const Vector& r, const Vector& p, const Matrix& A) const;
    double ComputeBeta(const Vector& r, const Vector& rNext) const;
    bool CheckConvergence(const Vector& r, double tolerance = 1e-10) const;
};

// LinearSystem constructor
LinearSystem::LinearSystem(const Matrix& A, const Vector& b) {
    // Check if the system is valid
    if (!A.IsSquare()) {
        throw std::invalid_argument("Matrix A must be square");
    }
    if (A.GetNumRows() != b.GetSize()) {
        throw std::invalid_argument("Matrix A and vector b must have compatible dimensions");
    }

    mSize = A.GetNumRows();
    mpA = new Matrix(A);
    mpb = new Vector(b);
}

// LinearSystem destructor
LinearSystem::~LinearSystem() {
    delete mpA;
    delete mpb;
}

// LinearSystem Solve method using Gaussian elimination with partial pivoting
Vector LinearSystem::Solve() const {
    // Create copies of the matrix and vector to work with
    Matrix A = *mpA;
    Vector b = *mpb;
    Vector x(mSize);

    // Gaussian elimination with partial pivoting
    for (int k = 1; k <= mSize; k++) {
        // Find pivot
        int maxRow = k;
        double maxVal = std::abs(A(k, k));
        for (int i = k + 1; i <= mSize; i++) {
            double val = std::abs(A(i, k));
            if (val > maxVal) {
                maxVal = val;
                maxRow = i;
            }
        }

        // Check if matrix is singular
        if (maxVal < 1e-10) {
            throw std::runtime_error("Matrix is singular or nearly singular");
        }

        // Swap rows if necessary
        if (maxRow != k) {
            // Swap rows in matrix A
            for (int j = 1; j <= mSize; j++) {
                std::swap(A(k, j), A(maxRow, j));
            }
            // Swap elements in vector b (using 1-based indexing)
            std::swap(b(k), b(maxRow));
        }

        // Eliminate column k
        for (int i = k + 1; i <= mSize; i++) {
            double factor = A(i, k) / A(k, k);
            A(i, k) = 0.0;  // Explicitly set to zero for numerical stability
            for (int j = k + 1; j <= mSize; j++) {
                A(i, j) -= factor * A(k, j);
            }
            b(i) -= factor * b(k);
        }
    }

    // Back substitution
    for (int i = mSize; i >= 1; i--) {
        double sum = 0.0;
        for (int j = i + 1; j <= mSize; j++) {
            sum += A(i, j) * x(j);
        }
        x(i) = (b(i) - sum) / A(i, i);
    }

    return x;
}

// PosSymLinSystem constructor
PosSymLinSystem::PosSymLinSystem(const Matrix& A, const Vector& b) : LinearSystem(A, b) {
    // Check if the matrix is symmetric
    if (!A.IsSymmetric()) {
        throw invalid_argument("Matrix must be symmetric for PosSymLinSystem");
    }

    // Additional check for positive definiteness (optional)
    // This is a simple check that might not catch all cases
    for (int i = 1; i <= A.GetNumRows(); i++) {
        if (A(i, i) <= 0) {
            throw invalid_argument("Matrix must be positive definite");
        }
    }
}

// PosSymLinSystem Solve method using conjugate gradient method
Vector PosSymLinSystem::Solve() const {
    const int maxIterations = mSize;  // Maximum number of iterations
    const double tolerance = 1e-10;   // Convergence tolerance
    const double minResidual = 1e-15; // Minimum residual to prevent division by zero

    Vector x(mSize);                  // Initial guess (zero vector)
    Vector r = *mpb - (*mpA) * x;    // Initial residual
    Vector p = r;                     // Initial search direction
    Vector rNext(mSize);             // Next residual
    Vector Ap(mSize);                // A * p

    double initialResidual = r.Norm();
    if (initialResidual < minResidual) {
        return x;  // Initial guess is already solution
    }

    for (int iter = 0; iter < maxIterations; iter++) {
        Ap = (*mpA) * p;
        double alpha = ComputeAlpha(r, p, *mpA);
        
        // Update solution and residual
        x += alpha * p;
        rNext = r - alpha * Ap;

        // Check convergence
        if (CheckConvergence(rNext, tolerance * initialResidual)) {
            return x;
        }

        // Update search direction
        double beta = ComputeBeta(r, rNext);
        if (std::abs(beta) < minResidual) {
            // If beta is too small, restart the algorithm
            r = *mpb - (*mpA) * x;
            p = r;
        } else {
            p = rNext + beta * p;
            r = rNext;
        }
    }

    throw std::runtime_error("Conjugate gradient method did not converge within " + 
                           std::to_string(maxIterations) + " iterations");
}

// Helper function to compute alpha in conjugate gradient method
double PosSymLinSystem::ComputeAlpha(const Vector& r, const Vector& p, const Matrix& A) const {
    Vector Ap = A * p;
    double pAp = p.DotProduct(Ap);
    if (std::abs(pAp) < 1e-15) {
        throw std::runtime_error("Matrix is not positive definite");
    }
    return r.DotProduct(r) / pAp;
}

// Helper function to compute beta in conjugate gradient method
double PosSymLinSystem::ComputeBeta(const Vector& r, const Vector& rNext) const {
    double rDotR = r.DotProduct(r);
    if (std::abs(rDotR) < 1e-15) {
        return 0.0;  // Prevent division by zero
    }
    return rNext.DotProduct(rNext) / rDotR;
}

// Helper function to check convergence in conjugate gradient method
bool PosSymLinSystem::CheckConvergence(const Vector& r, double tolerance) const {
    return r.Norm() <= tolerance;
}

#endif // LINEAR_SYSTEM_H