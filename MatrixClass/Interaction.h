#pragma once
/**************************************************************************
*   Author:  Pihai Sun                                                    *
*   Country: China                                                        *
*   University: Qingdao University                                        *
*   College: School of Computer Science and Technology                    *
*   Course teacher: Guodong Wang                                          *
*   Description:                                                          *
*   I will declare functions related to user interaction in this file     *
**************************************************************************/
#include "Matrix.h"
#ifndef __Interaction_h__
#define __Interaction_h__
// Print welcome message and program introduction
void Start();
bool AskContinue();
bool AskCreat(Matrix& cur_matrix);
void AskOp(Matrix& cur_matrix, Matrix& op_matrix);

#endif