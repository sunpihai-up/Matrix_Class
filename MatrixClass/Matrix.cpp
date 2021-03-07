/************************************************************************
*   Author:  Pihai Sun                                                  *
*   Country: China                                                      *
*   University: Qingdao University                                      *
*   College: School of Computer Science and Technology                  *
*   Course teacher: Guodong Wang                                        *
*   Description:                                                        *
*   I will give the definition of the function declared                 *
*   in Matrix.h in this file                                            *
************************************************************************/
#include "Matrix.h"

/* ************ Start of Member Function ************ */

Matrix::Matrix()
{
    // No-argument constructor
    data = NULL;
    rowNum = 0;
    columnNum = 0;
}

Matrix::Matrix(int row, int column)
{
    rowNum = row;
    columnNum = column;

    data = NULL;
    GetDataSpace();
}

Matrix::Matrix(double* da, int row, int column)
{
    rowNum = row;
    columnNum = column;

    data = NULL;
    GetDataSpace();

    CopyData(da);
}

Matrix::Matrix(const Matrix& m)
{
    // Copy constructor
    *this = m;
    //cout << this << "))))" << endl;
}

Matrix::Matrix(int dimension)
{
    rowNum = dimension;
    columnNum = dimension;

    data = NULL;
    GetDataSpace();

    for (int i = 0; i < rowNum; i++)
    {
        for (int j = 0; j < columnNum; j++)
        {
            int Serial = GetSerialNum(i, j, rowNum);
            data[Serial] = (i == j) ? 1 : 0;
        }
    }
}

Matrix::~Matrix()
{
    delete[] data;
    data = NULL;

    rowNum = 0;
    columnNum = 0;
}

void Matrix::CreatNewMatrix()
{
    int row = 0, col = 0;
    while (row < 1 || col < 1)
    {
        cout << "Rows = ";
        cin >> row;
        cout << "Columns = ";
        cin >> col;
        if (row < 1 || col < 1)
        {
            ReportError(IniError);
        }
    }

    rowNum = row;
    columnNum = col;
    GetDataSpace();

    cout << "Please enter the matrix elements in order from top left to bottom right\n";
    for (int i = 0; i < rowNum * columnNum; i++)
    {
        if (i % columnNum == 0)
        {
            cout << "Please enter the "
                << i / columnNum
                << "-th row element of the matrix :";
        }
        cin >> data[i];
    }

    cout << "Now that you have this matrix: \n";
    print();
    cout << endl;
}

double& Matrix::operator()(int row, int column)
{
    int Serial_num = GetSerialNum(row, column, columnNum);
    return data[Serial_num];
}

double& Matrix::operator()(int Serial_num)
{
    return data[Serial_num];
}

void Matrix::operator*=(const double constant)
{
    for (int i = 0; i < rowNum * columnNum; i++)
    {
        data[i] *= constant;
    }
}

void Matrix::GetDataSpace()
{

    // Assign a rowNum*columnNum space to the matrix
    //cout << data << endl;
    delete[] data;

    data = new double[rowNum * columnNum];
}

void Matrix::CopyData(const double* da)
{
    // Copy the data from the array into the matrix
    for (int i = 0; i < rowNum * columnNum; i++)
    {
        data[i] = da[i];
    }
}

void Matrix::print()
{
    cout << endl;
    for (int i = 0; i < rowNum; i++)
    {
        for (int j = 0; j < columnNum; j++)
        {
            cout << setiosflags(ios::fixed) << setprecision(PRECISION)
                << setiosflags(ios::left);
            cout << setw(15) << data[GetSerialNum(i, j, columnNum)];
        }
        cout << endl;
    }
}

void Matrix::operator=(const Matrix& m)
{
    // Overload the matrix equal sign operator
    rowNum = m.rowNum;
    columnNum = m.columnNum;
    // cout << &m << "))))(" << endl;
    //cout << "&&" << this << endl;
    //cout << this->data << "***" << endl;
    GetDataSpace();

    CopyData(m.data);

    return;
}

void Matrix::operator+=(const Matrix& m)
{
    if (!isSameType(*this, m))
    {
        ReportError(AddSubError);
        return;
    }

    for (int i = 0; i < rowNum * columnNum; i++)
    {
        data[i] += m.data[i];
    }
    return;
}

void Matrix::operator-=(const Matrix& m)
{
    if (!isSameType(*this, m))
    {
        ReportError(AddSubError);
        return;
    }

    for (int i = 0; i < rowNum * columnNum; i++)
    {
        data[i] -= m.data[i];
    }
    return;
}

void Matrix::operator*=(const Matrix& m)
{
    if (columnNum != m.rowNum)
    {
        ReportError(MulError);
        return;
    }

    Matrix tem_matrix(rowNum, m.columnNum);
    for (int i = 0; i < tem_matrix.rowNum; i++)
    {
        for (int j = 0; j < tem_matrix.columnNum; j++)
        {
            tem_matrix(i, j) = 0;
            for (int k = 0; k < tem_matrix.rowNum; k++)
            {
                int pos_m = GetSerialNum(k, j, m.columnNum);
                tem_matrix(i, j) += (*this)(i, k) * m.data[pos_m];
            }
        }
    }
    *this = tem_matrix;
    return;
}

void Matrix::SetIndentity(int dimension)
{
    rowNum = columnNum = dimension;
    GetDataSpace();

    for (int i = 0; i < rowNum; i++)
    {
        for (int j = 0; j < columnNum; j++)
        {
            if (i == j)
            {
                (*this)(i, j) = 1;
            }
            else
            {
                (*this)(i, j) = 0;
            }
        }
    }
}

// TODO: Explain Transpose
void Matrix::transpose()
{
    for (int cur_pos = 0; cur_pos < rowNum * columnNum; cur_pos++)
    {
        if (is_Traversed(cur_pos))
        {
            continue;
        }

        MoveDataInCir(cur_pos);
    }

    swap(rowNum, columnNum);
}

int Matrix::GetNexPos(int cur_pos)
{
    return (cur_pos % columnNum) * rowNum + cur_pos / columnNum;
}

int Matrix::GetPrePos(int cur_pos)
{
    return (cur_pos % rowNum) * columnNum + cur_pos / rowNum;
}

void Matrix::MoveDataInCir(int cur_pos)
{
    double temp = data[cur_pos];
    int cur = cur_pos;
    int pre = GetPrePos(cur);
    while (pre != cur_pos)
    {
        data[cur] = data[pre];
        cur = pre;
        pre = GetPrePos(cur);
    }
    data[cur] = temp;
}

bool Matrix::is_Traversed(int cur_pos)
{
    int next = GetNexPos(cur_pos);
    while (next > cur_pos)
    {
        next = GetNexPos(next);
    }

    if (next == cur_pos)
    {
        return false;
    }
    return true;
}

bool Matrix::GaussianInverse(Matrix& inv)
{
    if (rowNum != columnNum)
    {
        ReportError(InvError);
        return false;
    }

    Matrix temp_matrix(*this);
    inv.SetIndentity(rowNum);
    if (!Gaussian(temp_matrix, inv))
    {
        ReportError(InvError);
        return false;
    }

    return true;
}

// TODO: Explain LUP_SOLVE
bool Matrix::LUP_solve(Vector& B, Vector& Ans)
{
    int dimension = rowNum;
    Vector P(B.length);
    Matrix temp_matrix(*this);
    Init_P(P);
    if (!LUP_decomposition(temp_matrix, P))
    {
        ReportError(SinError);
        return false;
    }

    Matrix L(dimension);
    Matrix U(dimension);
    Init_L(L);
    Init_U(U);

    GetL(temp_matrix, L);
    GetU(temp_matrix, U);

    Vector temp_Ans(dimension);
    Vector Y(dimension);

    ForwardsSubstitution(L, Y, P, B);
    ReverseSubstitution(U, temp_Ans, Y);

    Ans = temp_Ans;
    return true;
}

bool Matrix::LUP_decomposition(Matrix& temp_matrix, Vector& P)
{
    // The vector P is already initialized
    int dimension = rowNum;

    for (int k = 0; k < dimension; k++)
    {
        double pivot = 0;
        int kk = k;
        for (int i = k; i < dimension; i++)
        {
            if (abs(temp_matrix(i, k)) > pivot)
            {
                pivot = abs(temp_matrix(i, k));
                kk = i;
            }
        }

        if (pivot == 0)
        {
            return false;
        }

        swap(P(k), P(kk));

        for (int i = 0; i < dimension; i++)
        {
            swap(temp_matrix(k, i), temp_matrix(kk, i));
        }

        for (int i = k + 1; i < dimension; i++) //row
        {
            temp_matrix(i, k) /= temp_matrix(k, k);
            for (int j = k + 1; j < dimension; j++) //col
            {
                temp_matrix(i, j) -= temp_matrix(i, k) * temp_matrix(k, j);
            }
        }
    }
    return true;
}

void Matrix::Init_P(Vector& P)
{
    for (int i = 0; i < P.length; i++)
    {
        P(i) = i;
    }
}

void Matrix::Init_L(Matrix& L)
{
    memset(L.data, 0, sizeof(L.data));
    for (int i = 0; i < L.rowNum; i++)
    {
        L(i, i) = 1;
    }
}

void Matrix::Init_U(Matrix& U)
{
    memset(U.data, 0, sizeof(U.data));
}

// TODO: Explain Forward
void Matrix::ForwardsSubstitution(Matrix& L, Vector& Y, Vector& P, Vector& B)
{
    int dimension = L.rowNum;
    for (int i = 0; i < dimension; i++)
    {
        double sub = 0;
        for (int j = 0; j < i; j++)
        {
            sub += L(i, j) * Y(j);
        }
        Y(i) = B((int)P(i)) - sub;
    }
}

// TODO: Explain Reverse
void Matrix::ReverseSubstitution(Matrix& U, Vector& Ans, Vector& Y)
{
    int dimension = U.rowNum;
    int i, j;
    for (i = dimension - 1; i >= 0; i--)
    {
        double sub = 0;
        for (j = i + 1; j < dimension; j++)
        {
            sub += U(i, j) * Ans(j);
        }

        Ans(i) = (Y(i) - sub) / U(i, i);
    }
}

void Matrix::GetL(Matrix& temp_matrix, Matrix& L)
{
    int dimension = L.rowNum;
    for (int i = 1; i < dimension; i++)
    {
        for (int j = 0; j < i; j++)
        {
            L(i, j) = temp_matrix(i, j);
        }
    }
}

void Matrix::GetU(Matrix& temp_matrix, Matrix& U)
{
    int dimension = U.rowNum;
    for (int i = 0; i < dimension; i++)
    {
        for (int j = i; j < dimension; j++)
        {
            U(i, j) = temp_matrix(i, j);
        }
    }
}

bool Matrix::canSolve(Vector& B)
{
    if (rowNum != B.length)
    {
        ReportError(SolError);
        return false;
    }

    return true;
}

double Matrix::GetDeterminant()
{
    int dimension = rowNum;
    Vector P(dimension);
    Matrix temp_matrix(*this);
    Init_P(P);

    // Singular matrix
    if (!LUP_decomposition(temp_matrix, P))
    {
        return 0;
    }

    // Get the lower triangular matrix U
    Matrix U(dimension);
    Init_U(U);
    GetU(temp_matrix, U);

    // Calculate the value of the determinant
    double det = 1;
    for (int i = 0; i < dimension; i++)
    {
        det *= U(i, i);
    }

    if (GetReversePairInP(P))
    {
        det = -det;
    }

    return det;
}

bool Matrix::GetReversePairInP(Vector& P)
{
    int RePairNum = 0;
    for (int i = 0; i < P.length; i++)
    {
        for (int j = i + 1; j < P.length; j++)
        {
            if (P(i) > P(j))
            {
                RePairNum++;
            }
        }
    }

    if (RePairNum % 2 == 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool Matrix::canDet()
{
    if (rowNum != columnNum)
    {
        ReportError(DetError);
        return false;
    }
    return true;
}

int Matrix::Rank()
{
    int none_zero = 0;
    int result = 0;
    Matrix temp_matrix(*this);
    standard_echelon(temp_matrix);

    for (int i = 0; i < rowNum; i++)
    {
        for (int j = 0; j < columnNum; j++)
        {
            if (temp_matrix(i, j) != 0)
            {
                none_zero = 1;
                break;
            }
        }

        if (none_zero == 1)
        {
            result++;
        }
        none_zero = 0;
    }
    return result;
}

void Matrix::standard_echelon(Matrix& temp_matrix)
{
    int i, j, k, l;
    double times;
    int total[1000];
    memset(total, 0, sizeof(total));

    for (i = 0; i < temp_matrix.columnNum - 1; i++)
    {
        for (k = i + 1; k < temp_matrix.columnNum; k++)
        {
            j = 0;
            while (temp_matrix(i, j) == 0)
            {
                j++;
            }

            if (temp_matrix(i, j) != 0)
            {
                times = temp_matrix(k, j) / temp_matrix(i, j);
                for (j = 0; j < temp_matrix.columnNum; j++)
                {
                    temp_matrix(k, j) -= temp_matrix(i, j) * times;
                }
            }
        }
    }

    for (i = 0; i < temp_matrix.rowNum; i++)
    {
        j = 0;
        while (temp_matrix(i, j) == 0)
        {
            j++;
        }

        if (temp_matrix(i, j) != 0)
        {
            times = temp_matrix(i, j);
            for (j = 0; j < temp_matrix.columnNum; j++)
            {
                temp_matrix(i, j) /= times;
            }
        }
    }

    for (i = 0; i < temp_matrix.rowNum; i++)
    {
        for (j = 0; j < temp_matrix.columnNum; j++)
        {
            if (temp_matrix(i, j) == 0)
            {
                total[i]++;
            }
            else
            {
                break;
            }
        }
    }

    for (l = temp_matrix.rowNum - 1; l > 0; l--)
    {
        for (i = 0; i < l; i++)
        {
            if (total[l] < total[i])
            {
                for (j = 0; j < temp_matrix.columnNum; j++)
                {
                    swap(temp_matrix(l, j), temp_matrix(i, j));
                }
            }
        }
    }

    for (i = 0; i < temp_matrix.rowNum; i++)
    {
        j = 0;
        while (temp_matrix(i, j) == 0)
        {
            j++;
        }

        if (temp_matrix(i, j) != 0)
        {
            for (k = 0; k < i; k++)
            {
                times = temp_matrix(k, j) / temp_matrix(i, j);
                for (l = 0; l < temp_matrix.columnNum; l++)
                {
                    temp_matrix(k, l) -= times * temp_matrix(i, l);
                }
            }
        }
    }
}

/* ************ End Of Member Function ************* */

/* ************ Start Of Friend Function ************* */

bool isSameType(const Matrix& m, const Matrix& n)
{
    // Determine if the two matrices have the same shape
    if (m.rowNum == n.rowNum && m.columnNum == n.columnNum)
    {
        return true;
    }
    else
    {
        return false;
    }
}

Matrix operator+(const Matrix& m, const Matrix& n)
{
    Matrix SumOfMatrix(m.rowNum, m.columnNum);

    for (int i = 0; i < m.columnNum * m.rowNum; i++)
    {
        SumOfMatrix.data[i] = m.data[i] + n.data[i];
    }

    return SumOfMatrix;
}

Matrix operator-(const Matrix& m, const Matrix& n)
{
    Matrix SubOfMatrix(m.rowNum, m.columnNum);

    for (int i = 0; i < m.rowNum * m.columnNum; i++)
    {
        SubOfMatrix.data[i] = m.data[i] - m.data[i];
    }

    return SubOfMatrix;
}

Matrix operator*(const Matrix& m, const Matrix& n)
{
    Matrix tem_matriax(m.rowNum, n.columnNum);

    for (int i = 0; i < m.rowNum; i++)
    {
        for (int j = 0; j < n.columnNum; j++)
        {
            int pos = GetSerialNum(i, j, tem_matriax.columnNum);
            for (int k = 0; k < m.columnNum; k++)
            {
                int pos_m = GetSerialNum(i, k, m.columnNum);
                int pos_n = GetSerialNum(k, j, n.columnNum);
                tem_matriax.data[pos] = m.data[pos_m] * n.data[pos_n];
            }
        }
    }
    cout << n.data << "((" << endl;
    return tem_matriax;
}

bool canAddSub(const Matrix& m, const Matrix& n)
{
    if (!isSameType(m, n))
    {
        ReportError(AddSubError);
        return false;
    }
    return true;
}

bool canMul(const Matrix& m, const Matrix& n)
{
    if (m.columnNum != n.rowNum)
    {
        ReportError(MulError);
        return false;
    }
    return true;
}

int GetSerialNum(int row, int column, int colNum)
{
    return row * colNum + column;
}

bool Gaussian(Matrix& temp, Matrix& ans)
{
    int dimension = temp.rowNum;
    double MAX;
    int k;
    int i, j;

    for (i = 0; i < dimension; i++)
    {
        // Find the pivot
        MAX = temp(i, i);
        k = i;
        for (j = i + 1; j < dimension; j++)
        {
            if (fabs(temp(j, i)) > fabs(MAX))
            {
                MAX = temp(j, i);
                k = j;
            }
        }

        // If the row of the pivot is not the i-th row,
        // perform row swap
        if (k != i)
        {
            for (j = 0; j < dimension; j++)
            {
                swap(temp(i, j), temp(k, j));
                swap(ans(i, j), ans(k, j));
            }
        }

        // Determine whether the pivot is 0,
        // if so, then matrix is not a full-rank matrix,
        // and there is no inverse matrix
        if (temp(i, i) == 0)
        {
            return false;
        }

        // Eliminate the i-th column of A
        // and remove the elements of each row except row i
        double t = temp(i, i);
        for (j = 0; j < dimension; j++)
        {
            // The element on the main diagonal becomes 1
            temp(i, j) = temp(i, j) / t;
            ans(i, j) = ans(i, j) / t;
        }

        for (j = 0; j < dimension; j++)
        {
            if (j != i)
            {
                t = temp(j, i);
                for (k = 0; k < dimension; k++)
                {
                    temp(j, k) = temp(j, k) - temp(i, k) * t;
                    ans(j, k) = ans(j, k) - ans(i, k) * t;
                }
            }
        }
    }

    return true;
}

/* ************ End Of Friend Function ************ */

/* ************ Start Of Ordinary Function ************ */

void ReportError(int ErrorType)
{
    switch (ErrorType)
    {
    case 1:
        cout << "Error: Matrices that are not the same shape cannot be added or subtracted!\n";
        break;
    case 2:
        cout << "Error: The length of vector B must be equal to the number of rows of matrix A\n";
    case 3:
        cout << "Error: These two matrices can not perform multiplication operations!\n";
        break;
    case 4:
        cout << "Error: Singular matrix has no inverse matrix\n";
        break;
    case 5:
        cout << "Error: Singular Error!\n";
        break;
    case 6:
        cout << "Error: Min size = 1*1\n";
        break;
    case 7:
        cout << "Error: Plesse input yes or no!\n";
        break;
    case 8:
        cout << "Error: A determinant cannot be calculated for a non-square matrix!\n";
        break;
    case 9:
        cout << "Error: Invalid operation code!\n";
        break;
    default:
        break;
    }
}

/* ************ End Of Ordinary Function ************ */

/* ************ Vector ************ */

Vector::Vector()
{
    data = 0;
    length = 0;
}

Vector::Vector(int l)
{
    length = l;
    data = new double[l];
}

void Vector::CreatNewVector()
{
    cout << "The length of this vector: ";
    cin >> length;

    data = new double[length];
    cout << "The element of this vector:\n";
    for (int i = 0; i < length; i++)
    {
        cin >> data[i];
    }
}

Vector::Vector(const Vector& v)
{
    length = v.length;
    data = new double[length];
    for (int i = 0; i < length; i++)
    {
        data[i] = v.data[i];
    }
}

Vector::~Vector()
{
    delete[] data;
    length = 0;
    data = NULL;
}

void Vector::print()
{
    for (int i = 0; i < length; i++)
    {
        cout << data[i] << " ";
    }
    cout << endl;
}

void Vector::operator=(const Vector& v)
{
    delete data;
    length = v.length;
    data = new double[length];
    for (int i = 0; i < length; i++)
    {
        data[i] = v.data[i];
    }
}

double& Vector::operator()(int serial)
{
    return data[serial];
}