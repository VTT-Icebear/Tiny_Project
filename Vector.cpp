#include "Vector.h"


int main(){
    Vector A(4);
    Vector B(4);
    cin >> A;
    cout << "--------------------------------" << endl;
    cin >> B;
    cout << "--------------------------------" << endl;
    cout << "A + B: " << endl;
    cout << A + B << endl;
    cout << "--------------------------------" << endl;
    cout << "A - B: " << endl;
    cout << A - B << endl;
    cout << "--------------------------------" << endl;
    cout << "A * 2: " << endl;
    cout << A * 2 << endl;
    cout << "--------------------------------" << endl;
    cout << "A *= 2: " << endl;
    cout << (A *= 2) << endl;
    cout << "--------------------------------" << endl;
    cout << "Norm of A: " << endl;
    cout << A.Norm() << endl;
    cout << "--------------------------------" << endl;
    cout << "Dot product of A and B: " << endl;
    cout << A.DotProduct(B) << endl;
    cout << "--------------------------------" << endl;
    return 0;
}