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


typedef struct {
  float x ;
  float y ;
  char type ;
} MuacRecord ;

/* this routine expects an array of data which has already been sorted
   by y */
float muac( MuacRecord* in, int n1, int n2 ) ;

/* this routine expects unsorted data - the sorting is randomized, so
   using with sorted data will not cause problems */
float unsortedMuac( float* x1, float* y1, int n1, 
                    float* x2, float* y2, int n2 ) ;

