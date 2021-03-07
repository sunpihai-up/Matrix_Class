#include "Matrix.h"
#include "Interaction.h"

int main()
{
    Matrix cur_matrix, op_matrix;
    Start();
    cur_matrix.CreatNewMatrix();

    do
    {
        AskOp(cur_matrix, op_matrix);
    } while (AskContinue() && AskCreat(cur_matrix));

    system("pause");
    return 0;
}