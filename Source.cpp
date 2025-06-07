#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

struct DataRow {
    double MYCT, MMIN, MMAX, CACH, CHMIN, CHMAX, PRP;
};

// ===========================
// MA TRẬN UTILITIES
// ===========================
typedef vector<vector<double>> Matrix;
typedef vector<double> Vector;

// Nhân ma trận A x B
Matrix multiply(const Matrix& A, const Matrix& B) {
    int n = A.size(), m = B[0].size(), p = B.size();
    Matrix result(n, vector<double>(m, 0));
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < p; ++k)
            for (int j = 0; j < m; ++j)
                result[i][j] += A[i][k] * B[k][j];
    return result;
}

// Chuyển vị ma trận
Matrix transpose(const Matrix& A) {
    int n = A.size(), m = A[0].size();
    Matrix result(m, vector<double>(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            result[j][i] = A[i][j];
    return result;
}

// Nghịch đảo ma trận bằng Gauss-Jordan
Matrix inverse(Matrix A) {
    int n = A.size();
    Matrix I(n, vector<double>(n, 0));
    for (int i = 0; i < n; ++i) I[i][i] = 1;

    for (int i = 0; i < n; ++i) {
        // Tìm hàng lớn nhất để pivot
        double maxVal = abs(A[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (abs(A[k][i]) > maxVal) {
                maxVal = abs(A[k][i]);
                maxRow = k;
            }
        }
        swap(A[i], A[maxRow]);
        swap(I[i], I[maxRow]);

        // Chia cho pivot
        double pivot = A[i][i];
        for (int j = 0; j < n; ++j) {
            A[i][j] /= pivot;
            I[i][j] /= pivot;
        }

        // Khử các hàng khác
        for (int k = 0; k < n; ++k) {
            if (k == i) continue;
            double factor = A[k][i];
            for (int j = 0; j < n; ++j) {
                A[k][j] -= factor * A[i][j];
                I[k][j] -= factor * I[i][j];
            }
        }
    }

    return I;
}

// ===========================
// HỒI QUY TUYẾN TÍNH
// ===========================

Vector linearRegression(const Matrix& X, const Vector& Y) {
    Matrix Xt = transpose(X);
    Matrix XtX = multiply(Xt, X);
    Matrix XtX_inv = inverse(XtX);

    // Chuyển Y thành ma trận cột
    Matrix Y_col(Y.size(), vector<double>(1));
    for (int i = 0; i < Y.size(); ++i) Y_col[i][0] = Y[i];

    Matrix XtY = multiply(Xt, Y_col);
    Matrix B = multiply(XtX_inv, XtY);

    // Trả về vector β
    Vector beta(B.size());
    for (int i = 0; i < B.size(); ++i) beta[i] = B[i][0];
    return beta;
}

// ===========================
// ĐỌC DATA
// ===========================

vector<DataRow> loadData(string filename) {
    ifstream file(filename);
    string line;
    vector<DataRow> data;
    while (getline(file, line)) {
        DataRow row;
        stringstream ss(line);
        string temp;
        int skip = 2;
        while (skip--) getline(ss, temp, ','); // skip vendor & model
        ss >> row.MYCT; ss.ignore();
        ss >> row.MMIN; ss.ignore();
        ss >> row.MMAX; ss.ignore();
        ss >> row.CACH; ss.ignore();
        ss >> row.CHMIN; ss.ignore();
        ss >> row.CHMAX; ss.ignore();
        ss >> row.PRP; ss.ignore();
        getline(ss, temp); // skip ERP
        data.push_back(row);
    }
    return data;
}

// ===========================
// RMSE
// ===========================

double rmse(const Matrix& X, const Vector& Y, const Vector& beta) {
    double sum = 0;
    for (int i = 0; i < Y.size(); ++i) {
        double y_pred = 0;
        for (int j = 0; j < beta.size(); ++j)
            y_pred += X[i][j] * beta[j];
        sum += pow(y_pred - Y[i], 2);
    }
    return sqrt(sum / Y.size());
}

// ===========================
// MAIN
// ===========================

int main() {
    auto data = loadData("machine.data");

    // Shuffle & split
    random_shuffle(data.begin(), data.end());
    int n = data.size();
    int train_size = n * 0.8;

    Matrix X_train, X_test;
    Vector Y_train, Y_test;

    auto fillXY = [](const DataRow& row, Matrix& X, Vector& Y) {
        X.push_back({ row.MYCT, row.MMIN, row.MMAX, row.CACH, row.CHMIN, row.CHMAX });
        Y.push_back(row.PRP);
        };

    for (int i = 0; i < train_size; ++i) fillXY(data[i], X_train, Y_train);
    for (int i = train_size; i < n; ++i) fillXY(data[i], X_test, Y_test);

    Vector beta = linearRegression(X_train, Y_train);

    cout << "Beta coefficients:\n";
    for (double b : beta) cout << b << " ";
    cout << "\n";

    double err = rmse(X_test, Y_test, beta);
    cout << "RMSE on test set: " << err << endl;
    return 0;
}
