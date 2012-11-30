/**
 *
 * @file psbarrier.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * Barrier for algorithm mixing computation on tile/panel.
 *
 * @version 2.4.2
 * @author Mathieu Faverge
 * @date 2009-11-15
 *
 * @generated s Thu Sep 15 12:09:24 2011
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, float, m, n)

/***************************************************************************//**
 *  Barrier from tiles to panels
 **/
void plasma_psbarrier_tl2pnl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request)
{
    int m, n;
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
 
    for (n = 0; n < A.nt; n++)
    {
        /* Protection from previous GATHERV */
        QUARK_Insert_Task(plasma->quark, CORE_foo_quark, &task_flags,
                          sizeof(float)*A.mb*A.nb, A(0, n), INOUT,
                          0);

        for (m = 0; m < A.mt; m++)
        {
            QUARK_Insert_Task(plasma->quark, CORE_foo2_quark, &task_flags,
                              sizeof(float)*A.mb*A.nb, A(0, n), INOUT | GATHERV,
                              sizeof(float)*A.mb*A.nb, A(m, n), INOUT,
                              0);
        }

        /* Protection to next GATHERV */
        QUARK_Insert_Task(plasma->quark, CORE_foo_quark, &task_flags,
                          sizeof(float)*A.mb*A.nb, A(0, n), INOUT,
                          0);
    }
}

/***************************************************************************//**
 *  Barrier from panels to tiles
 **/
void plasma_psbarrier_pnl2tl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request)
{
    int m, n;
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
 
    for (n = 0; n < A.nt; n++)
    {
        /* Protection from previous GATHERV */
        QUARK_Insert_Task(plasma->quark, CORE_foo_quark, &task_flags,
                          sizeof(float)*A.mb*A.nb, A(0, n), INOUT,
                          0);

        for (m = 0; m < A.mt; m++)
        {
            QUARK_Insert_Task(plasma->quark, CORE_foo2_quark, &task_flags,
                              sizeof(float)*A.mb*A.nb, A(0, n), INPUT,
                              sizeof(float)*A.mb*A.nb, A(m, n), INOUT,
                              0);
        }
    }
}

