#pragma once
/*******************************************************************************
*   Author:  Pihai Sun                                                         *
*   Country: China                                                             *
*   University: Qingdao University                                             *
*   College: School of Computer Science and Technology                         *
*   Course teacher: Guodong Wang                                               *
*   Description:                                                               *
*   I will declare matrix classes                                              *
*   and functions related to matrix operations in matrix.h                     *
*******************************************************************************/
#ifndef _Matrix_h_
#define _Matrix_h_
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <algorithm>
using namespace std;
const int AddSubError = 1;
const int SolError = 2;
const int MulError = 3;
const int InvError = 4;
const int SinError = 5;
const int IniError = 6;
const int ConError = 7;
const int DetError = 8;
const int OpError = 9;
// How many decimal places after the decimal output
const int PRECISION = 5;
class Vector;
// TODO:Show Class
class Matrix
{
private:
    double* data;
    int rowNum;
    int columnNum;

public:
    /*******************************
    *        Constructors          *
    *******************************/
    Matrix();
    // Initialize with row and columb
    Matrix(int row, int column);
    // Initialize with array,row and column
    Matrix(double* da, int row, int column);
    // Copy Constructor
    Matrix(const Matrix& m);
    // Create an identity matrix
    Matrix(int dimension);

    /**************************
    *       Destructor        *
    **************************/
    ~Matrix();

    /*********************************
    *        Basic operations        *
    *********************************/
    // Allocation of space for matrix classes
    void GetDataSpace();
    // Copy the data from the array into the matrix
    void CopyData(const double* da);
    // Output Matrix Class
    void print();
    // Get the element according row and column
    double& operator()(int, int);
    // Get the element according the Serial num in data
    double& operator()(int);
    // Creat a new matrix
    void CreatNewMatrix();
    // Set matrix to identity matrix
    void SetIndentity(int dimension);

    /**************************************
    *     Basic arithmetic operations     *
    **************************************/
    void operator=(const Matrix& m);
    void operator+=(const Matrix& m);
    void operator-=(const Matrix& m);
    void operator*=(const Matrix& m);
    void operator*=(const double constant);

    /**********************************************
    *              Other operation                *
    * Transpose, Inverse, Determinant, Rank, Ax=b *
    **********************************************/

    /* ************* Matrix transpose in situ ************* */
    // Transpose call function
    void transpose();
    // Get the predecessor of the transpose cycle
    int GetPrePos(int cur_pos);
    // Find the successor of the transpose cycle
    int GetNexPos(int cur_pos);
    // Move the elements in the transpose loop
    void MoveDataInCir(int cur_pos);
    // Check whether the transpose loop has been traversed
    bool is_Traversed(int cur_pos);
    /* ******* End of functions related to transpose ******* */

    /* ***** Gaussian elimination inverse matrix ***** */
    bool GaussianInverse(Matrix& inv);
    //friend bool Gaussian(Matrix &temp, Matrix &ans);
    /* ***** End of function related to inverse ***** */

    /* ***** LUP decomposition method to solve linear equations ***** */
    bool LUP_solve(Vector& B, Vector& Ans);
    // LUP decomposition to find L, U, P
    bool LUP_decomposition(Matrix& temp_matrix, Vector& P);
    // Initialize L, U, P
    void Init_P(Vector& P);
    void Init_L(Matrix& L);
    void Init_U(Matrix& U);
    void GetL(Matrix& temp_matrix, Matrix& L);
    void GetU(Matrix& temp_matrix, Matrix& U);
    // Forward substitution to solve Ly=Pb
    void ForwardsSubstitution(Matrix& L, Vector& Y, Vector& P, Vector& B);
    // Backward replacement to solve Ux=y
    void ReverseSubstitution(Matrix& U, Vector& Ans, Vector& Y);
    bool canSolve(Vector& B);
    /* ***** End of LUP_Solve linear equations ***** */

    /* ***** Get the determinant of the matrix ***** */
    // Some functions in LUP decomposition are used to implement this function
    double GetDeterminant();
    // Calculate the number of reverse pairs.
    // Determine whether the value of the determinant needs to be negative
    bool GetReversePairInP(Vector& P);
    bool canDet();
    /* ***** End of Get determinant ***** */

    /* ***** Rank ***** */
    int Rank();
    void standard_echelon(Matrix& temp_matrix);
    /* *****End of rank ***** */

    /**************************************
    *     Friend function declaration     *
    **************************************/
    friend bool isSameType(const Matrix& m, const Matrix& n);
    friend Matrix operator+(const Matrix& m, const Matrix& n);
    friend Matrix operator-(const Matrix& m, const Matrix& n);
    friend Matrix operator*(const Matrix& m, const Matrix& n);
    friend int GetSerialNum(int row, int column, int colNum);
    friend bool Gaussian(Matrix& temp, Matrix& ans);
    friend bool canAddSub(const Matrix& m, const Matrix& n);
    friend bool canMul(const Matrix& m, const Matrix& n);
};

// Check whether the two matrices with the same shape
bool isSameType(const Matrix& m, const Matrix& n);
// void CopyShape(Matrix &);
// Friend overloading of the operator
Matrix operator+(const Matrix& m, const Matrix& n);
Matrix operator-(const Matrix& m, const Matrix& n);
Matrix operator*(const Matrix& m, const Matrix& n);
bool canAddSub(const Matrix& m, const Matrix& n);
bool canMul(const Matrix& m, const Matrix& n);
// Get the position of the element of the specified row and column in a one-dimensional array
int GetSerialNum(int row, int column, int coluNum);
// Gaussian elimination method operation processa
bool Gaussian(Matrix& temp, Matrix& ans);
// Print various error messages
void ReportError(int ErrorType);

class Vector
{
private:
    int length;
    double* data;

public:
    Vector();
    Vector(int l);
    Vector(const Vector& v);
    Vector(double* d, int l);
    void CreatNewVector();
    ~Vector();
    double& operator()(int serial);
    void operator=(const Vector& v);
    void print();
    friend Matrix;
};

#endif