#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
using namespace std;
#include <cassert>
#include <cmath>

class Vector {
private:
    int mSize;      // Size of the vector
    double* mData;  // Pointer to the data array

public:
    // Constructors and destructor
    Vector(int size);                    // Constructor with size
    Vector(const Vector& other);         // Copy constructor
    ~Vector();                          // Destructor

    // Assignment operator
    Vector& operator=(const Vector& other);

    // Access operators
    double& operator[](int i);           // 0-based indexing with bounds checking
    const double& operator[](int i) const;
    double& operator()(int i);           // 1-based indexing
    const double& operator()(int i) const;

    // Unary operators
    Vector operator+() const;            // Unary plus
    Vector operator-() const;            // Unary minus

    // Binary operators
    Vector operator+(const Vector& other) const;  // Vector addition
    Vector operator-(const Vector& other) const;  // Vector subtraction
    Vector operator*(double scalar) const;        // Scalar multiplication
    friend Vector operator*(double scalar, const Vector& vec);  // Scalar multiplication (left)

    // Compound assignment operators
    Vector& operator+=(const Vector& other);
    Vector& operator-=(const Vector& other);
    Vector& operator*=(double scalar);

    // Utility functions
    int GetSize() const { return mSize; }
    double Norm() const;                 // Euclidean norm
    double DotProduct(const Vector& other) const;  // Dot product

    // Friend functions for I/O
    friend ostream& operator<<(ostream& os, const Vector& vec);
    friend istream& operator>>(istream& is, Vector& vec);
};

// Constructor
Vector::Vector(int size) : mSize(size) {
    assert(size > 0);
    mData = new double[size];
    for (int i = 0; i < size; i++) {
        mData[i] = 0.0;
    }
}

// Copy constructor
Vector::Vector(const Vector& other) : mSize(other.mSize) {
    mData = new double[mSize];
    for (int i = 0; i < mSize; i++) {
        mData[i] = other.mData[i];
    }
}

// Destructor
Vector::~Vector() {
    delete[] mData;
}

// Assignment operator
Vector& Vector::operator=(const Vector& other) {
    if (this != &other) {
        delete[] mData;
        mSize = other.mSize;
        mData = new double[mSize];
        for (int i = 0; i < mSize; i++) {
            mData[i] = other.mData[i];
        }
    }
    return *this;
}

// Access operators
double& Vector::operator[](int i) {
    assert(i >= 0 && i < mSize);
    return mData[i];
}

const double& Vector::operator[](int i) const {
    assert(i >= 0 && i < mSize);
    return mData[i];
}

double& Vector::operator()(int i) {
    assert(i >= 1 && i <= mSize);
    return mData[i-1];
}

const double& Vector::operator()(int i) const {
    assert(i >= 1 && i <= mSize);
    return mData[i-1];
}

// Unary operators
Vector Vector::operator+() const {
    return *this;
}

Vector Vector::operator-() const {
    Vector result(mSize);
    for (int i = 0; i < mSize; i++) {
        result.mData[i] = -mData[i];
    }
    return result;
}

// Binary operators
Vector Vector::operator+(const Vector& other) const {
    assert(mSize == other.mSize);
    Vector result(mSize);
    for (int i = 0; i < mSize; i++) {
        result.mData[i] = mData[i] + other.mData[i];
    }
    return result;
}

Vector Vector::operator-(const Vector& other) const {
    assert(mSize == other.mSize);
    Vector result(mSize);
    for (int i = 0; i < mSize; i++) {
        result.mData[i] = mData[i] - other.mData[i];
    }
    return result;
}

Vector Vector::operator*(double scalar) const {
    Vector result(mSize);
    for (int i = 0; i < mSize; i++) {
        result.mData[i] = mData[i] * scalar;
    }
    return result;
}

Vector operator*(double scalar, const Vector& vec) {
    return vec * scalar;
}

// Compound assignment operators
Vector& Vector::operator+=(const Vector& other) {
    assert(mSize == other.mSize);
    for (int i = 0; i < mSize; i++) {
        mData[i] += other.mData[i];
    }
    return *this;
}

Vector& Vector::operator-=(const Vector& other) {
    assert(mSize == other.mSize);
    for (int i = 0; i < mSize; i++) {
        mData[i] -= other.mData[i];
    }
    return *this;
}

Vector& Vector::operator*=(double scalar) {
    for (int i = 0; i < mSize; i++) {
        mData[i] *= scalar;
    }
    return *this;
}

// Utility functions
double Vector::Norm() const {
    double sum = 0.0;
    for (int i = 0; i < mSize; i++) {
        sum += mData[i] * mData[i];
    }
    return sqrt(sum);
}

double Vector::DotProduct(const Vector& other) const {
    assert(mSize == other.mSize);
    double sum = 0.0;
    for (int i = 0; i < mSize; i++) {
        sum += mData[i] * other.mData[i];
    }
    return sum;
}

// I/O operators
ostream& operator<<(ostream& os, const Vector& vec) {
    os << "[";
    for (int i = 0; i < vec.mSize; i++) {
        os << vec.mData[i];
        if (i < vec.mSize - 1) os << ", ";
    }
    os << "]";
    return os;
}

istream& operator>>(istream& is, Vector& vec) {
    for (int i = 0; i < vec.mSize; i++) {
        is >> vec.mData[i];
    }
    return is;
}

#endif // VECTOR_H 