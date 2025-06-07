#include "Matrix.h"

int main(){
    Matrix A(3, 3);
    Matrix B(3, 3);

    cin >> A;
    cout << "--------------------------------" << endl;
    cin >> B;
    cout << "--------------------------------" << endl;
    cout << "A + B: " << endl;
    cout << A + B << endl;
    cout << "--------------------------------" << endl;
    cout << "A * B: " << endl;
    cout << A * B << endl;
    cout << "--------------------------------" << endl;
    cout << "A * 2: " << endl;
    cout << A * 2 << endl;
    cout << "--------------------------------" << endl;
    cout << "A += B: " << endl;
    cout << (A += B) << endl;
    cout << "--------------------------------" << endl;
    cout << "A -= B: " << endl;
    cout << (A -= B) << endl;
    cout << "--------------------------------" << endl;
    cout << "A *= 2: " << endl;
    cout << (A *= 2) << endl;
    cout << "--------------------------------" << endl;
    cout << "Determinant of A: " << endl;
    cout << A.Determinant() << endl;
    cout << "--------------------------------" << endl;
    cout << "Inverse of A: " << endl;
    cout << A.Inverse() << endl;
    cout << "--------------------------------" << endl;
    cout << "Determinant of B: " << endl;
    cout << B.Determinant() << endl;
    cout << "--------------------------------" << endl;
    cout << "Inverse of B: " << endl;
    cout << B.Inverse() << endl;
    cout << "--------------------------------" << endl;
    return 0;
}