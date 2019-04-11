#include <cstring>

#ifndef LIBRARY_H
#define LIBRARY_H

double read_timer();
void fill( double *p, int n );
void printMatrix( double *A, int n );

int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
long long read_long_long( int argc, char **argv, const char *option, long long default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif 
