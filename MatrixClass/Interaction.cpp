/*************************************************************************
*   Author:  Pihai Sun                                                   *
*   Country: China                                                       *
*   University: Qingdao University                                       *
*   College: School of Computer Science and Technology                   *
*   Course teacher: Guodong Wang                                         *
*   Description:                                                         *
*   I will give the definition of the function declared                  *
*   by Interactive.h in this file                                        *
*************************************************************************/
#include "Interaction.h"
#include "Matrix.h"

void Start()
{
    cout << "******Welcome to use the matrix calculator written by Pihai Sun******\n";
    cout << "This matrix calculator will provide you with functions "
        << "such as basic matrix operations, matrix inversion, \n"
        << "matrix rank, matrix transposition, matrix determinant solving, "
        << "and linear equations solving.\n";
    cout << "The creative background of this matrix calculator is coursework.\n";
    cout << "So it is inevitable to have various shortcomings.\n";
    cout << "If you have good suggestions for improvement, you can contact me by email sunpihai@qq.com\n";
    cout << "******Wish you a happy life******\n\n";

    cout << "Matrix calculator operation code: \n";
    cout << "Enter 1 for matrix multiplication.\n";
    cout << "Enter 2 for matrix addition.\n";
    cout << "Enter 3 for matrix subtraction.\n";
    cout << "Enter 4 to find the matrix determinant.\n";
    cout << "Enter 5 to transpose the matrix.\n";
    cout << "Enter 6 to find the inverse of your matrix.\n";
    cout << "Enter 7 to find the matrix rank.\n";
    cout << "Enter 8 to solve linear equations.\n\n";

    cout << "Start by making your first matrix!\n";
}

bool AskContinue()
{
    string choice;
    while (1)
    {
        cout << "Continueï¼?ï¼ˆEnter yes or no) ";
        cin >> choice;
        if (choice == "no")
        {
            return false;
        }
        else if (choice == "yes")
        {
            return true;
        }
        else
        {
            ReportError(ConError);
        }
    }
}

bool AskCreat(Matrix& cur_matrix)
{
    string choice;
    cout << "Now enter yes to creat a new matrix or no to use current matrix.\n";
    cin >> choice;
    while (1)
    {
        if (choice == "yes")
        {
            cur_matrix.CreatNewMatrix();
            return true;
        }
        else if (choice == "no")
        {
            return true;
        }
        else
        {
            ReportError(ConError);
        }
    }
}

void AskOp(Matrix& cur_matrix, Matrix& op_matrix)
{
    string choice;
    while (1)
    {
        cout << "Please enter your operation choice: ";
        cin >> choice;
        if (choice == "1")
        {
            cout << "You need to create a new matrix for this operation\n";
            op_matrix.CreatNewMatrix();
            if (canMul(cur_matrix, op_matrix))
            {
                Matrix ans(cur_matrix);
                ans.print();
                ans *= op_matrix;
                ans.print();
            }
        }
        else if (choice == "2")
        {
            cout << "You need to create a new matrix for this operation\n";
            op_matrix.CreatNewMatrix();
            if (canAddSub(cur_matrix, op_matrix))
            {
                op_matrix += cur_matrix;
                op_matrix.print();
            }
        }
        else if (choice == "3")
        {
            cout << "You need to create a new matrix for this operation\n";
            op_matrix.CreatNewMatrix();
            if (canAddSub(cur_matrix, op_matrix))
            {
                op_matrix -= cur_matrix;
                op_matrix.print();
            }
        }
        else if (choice == "4")
        {
            if (cur_matrix.canDet())
            {
                double result = cur_matrix.GetDeterminant();
                cout << "The det of this matrix: " << result << endl;
            }
        }
        else if (choice == "5")
        {
            op_matrix = cur_matrix;
            op_matrix.transpose();
            op_matrix.print();
        }
        else if (choice == "6")
        {
            Matrix inv;
            bool result = cur_matrix.GaussianInverse(inv);
            if (result)
            {
                inv.print();
            }
        }
        else if (choice == "7")
        {
            int rank = cur_matrix.Rank();
            cout << "The rank of this matrix: " << rank << endl;
        }
        else if (choice == "8")
        {
            cout << "You need to enter a vector b to do this operation.\n";
            Vector B;
            B.CreatNewVector();

            if (cur_matrix.canSolve(B))
            {
                Vector Ans;
                bool result = cur_matrix.LUP_solve(B, Ans);
                if (result)
                {
                    Ans.print();
                }
            }
        }
        else
        {
            ReportError(OpError);
            continue;
        }
        break;
    }
    return;
}