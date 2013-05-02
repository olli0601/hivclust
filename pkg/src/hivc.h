/** \file nabc_fun.h
    \brief Main file that provides all functions that are callable from R.
*/

#ifndef NABC_FUN_H_
#define NABC_FUN_H_

#include <R.h>
#include <Rinternals.h>


extern "C" {

void hivc_printdna(unsigned char *x, int *n);
void hivc_dist_ambiguous_dna(unsigned char *x1, unsigned char *x2, int *n, double *ans);

}


#endif /*NABC_FUN_H_*/
