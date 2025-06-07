#include "LinearSystem.h"

int main() {

    Matrix A(3, 3);
    Vector b(3);

    // Set up a simple 3x3 system
    // 2x + y - z = 8
    // -3x - y + 2z = -11
    // -2x + y + 2z = -3
    A(1,1) = 2.0;  A(1,2) = 1.0;  A(1,3) = -1.0;
    A(2,1) = -3.0; A(2,2) = -1.0; A(2,3) = 2.0;
    A(3,1) = -2.0; A(3,2) = 1.0;  A(3,3) = 2.0;

    b(1) = 8.0;
    b(2) = -11.0;
    b(3) = -3.0;

    // Solve using regular linear system
    LinearSystem system(A, b);
    Vector x = system.Solve();

    cout << "Solution using Gaussian elimination:" << endl;
    cout << "x = " << x << endl;

    // Create a symmetric positive definite system
    Matrix B(3, 3);
    Vector c(3);

    // Set up a symmetric positive definite matrix
    B(1,1) = 4.0;  B(1,2) = 1.0;  B(1,3) = 1.0;
    B(2,1) = 1.0;  B(2,2) = 5.0;  B(2,3) = 2.0;
    B(3,1) = 1.0;  B(3,2) = 2.0;  B(3,3) = 6.0;

    c(1) = 1.0;
    c(2) = 2.0;
    c(3) = 3.0;

    // Solve using conjugate gradient method
    PosSymLinSystem posSystem(B, c);
    Vector y = posSystem.Solve();

    cout << "\nSolution using conjugate gradient:" << endl;
    cout << "y = " << y << endl;

    // Verify solutions
    Vector residual1 = A * x - b;
    Vector residual2 = B * y - c;

    cout << "\nResidual norms:" << endl;
    cout << "Gaussian elimination residual: " << residual1.Norm() << endl;
    cout << "Conjugate gradient residual: " << residual2.Norm() << endl;
    return 0;
}