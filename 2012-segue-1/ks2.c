/*

Muac - A fast algorithm for the 2d KS test (v1.1, 13 June 1999).

Copyright (C) 1998 Andrew Cooke  (Jara Software)

Comments to jara@andrewcooke.free-online.co.uk

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA

These conditions also apply to any other files distributed with this
program.

*/

#include "stdlib.h"
#include "stdio.h"

/* --- global structs and routines -------------------- */

#include "muac.h"


/* --- local structs and routines --------------------- */

void readFile( float **x, float **y, int *n, char *name ) ;


/* --- test muac routine ------------------------------ */

void main( int argc, char **argv ) {

  float *x1 ; /* pointer to data for first set */
  float *y1 ;
  float *x2 ; /* pointer to data for second set */
  float *y2 ;
  int n1 ;    /* number of points in each set */
  int n2 ; 
  
  if ( argc != 3 ) {
    fprintf( stderr, "\nUsage: ks2 file1 file2\n" ) ;
    fprintf( stderr, "Diagnostic messages to stderr; result to stdout.\n" ) ;
    fprintf( stderr, "Result is maximum fractional difference, NOT a transformed statistic.\n" ) ;
    fprintf( stderr, "2D KS test using fast muac algorithm (v1.1) (c) Andrew Cooke 1999\n" ) ;
    fprintf( stderr, "\nThis program is free software; you can redistribute it and/or modify\nit under the terms of the GNU General Public License as published by\nthe Free Software Foundation; either version 2 of the License, or\n(at your option) any later version.\n\nThis program is distributed in the hope that it will be useful,\nbut WITHOUT ANY WARRANTY; without even the implied warranty of\nMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\nGNU General Public License for more details.\n\nYou should have received a copy of the GNU General Public License\nalong with this program; if not, write to the Free Software\nFoundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.\n" ) ;
  } else {
    readFile( &x1, &y1, &n1, argv[1] ) ;
    readFile( &x2, &y2, &n2, argv[2] ) ;
    fprintf( stdout, "%f\n", unsortedMuac( x1, y1, n1, x2, y2, n2 )) ;
  }

}


/* --- read from a file ------------------------------- */

/* expect the format
   1.0  2.0
   3.0  4.0
   ...
   etc. */

/* memory is allocated in inc chunks - this shoud be set to a little
   larger than the number of points typically read for efficient
   running */

void readFile( float **x, float **y, int *n, char *name ) {

  FILE *in ;
  int nRead ;
  float xx ;
  float yy ;
  int size ;
  int inc ;

  *x = NULL ;
  *y = NULL ;
  *n = 0 ;
  inc = 100 ;
  size = 0 ;

#ifdef DEBUG
  fprintf( stderr, "reading data from %s\n", name ) ;
#endif
  in = fopen( name, "r" ) ;
  for (;;) {
    nRead = fscanf( in, " %f %f", &xx, &yy ) ;
    if ( nRead != 2 ) break ;
    *n = *n + 1 ;
    if ( *n > size ) {
      size += inc ;
      *x = (float*)realloc( *x, size * sizeof( float )) ;
      *y = (float*)realloc( *y, size * sizeof( float )) ;
    }
    (*x)[*n-1] = xx ;
    (*y)[*n-1] = yy ;
#ifdef DEBUG
    fprintf( stderr, "%d %f %f\n", *n, xx, yy ) ;
#endif
  }
}

