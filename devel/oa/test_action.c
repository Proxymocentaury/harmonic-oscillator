/*******************************************************************************
*
* File test_action.c
*
* Check for the action of OA
*
*******************************************************************************/
#define MAIN_PROGRAM
#define TEST_ACTION_C

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "global.h"
#include "action.h"

int main(void)
{   
    double action_value = 0.0;
    double action_true = 0.0; 
    int i = 0; 
    double dim = (double)N;


    /*TEST1 valori del vettore tutti pari a 1*/

    while (i < N){
        xx[i] = 1;
        i++;
    }

    action_true= N*M*W*W/2;
    action_value = soa();


    printf("valore teorico dell'azione: %.3f\n", action_true);
    printf("valore calcolato: %.3f\n", action_value);
    printf("differenza: %.3f\n", action_value - action_true);


    /*TEST2 valori del vettore pari a i*/
    
    i=0;

    while (i < N){
        xx[i] = i;
        i++;
    }
    
    action_true = 0.5*M*((dim-1) + (dim-1)*(dim-1)) 
            + 0.5*M*W*W*(dim-1)*dim*(2*dim-1)/6.0;
    
            action_value = soa(); 

    
    printf("valore teorico dell'azione: %.3f\n", action_true);
    printf("valore calcolato: %.3f\n", action_value);
    printf("differenza: %.3f\n", action_value - action_true);

    return 0;
}


