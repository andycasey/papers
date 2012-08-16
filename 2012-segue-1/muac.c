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

/* --- the main routines for a fast 2d ks test -------- */

#include "stdlib.h"
#include "stdio.h"

/* --- global structs and routines -------------------- */

#include "muac.h"
#undef DEBUG

/* --- local structs and routines --------------------- */

typedef struct node {
  struct node *parent ;
  struct node *left ;
  struct node *right ;
  MuacRecord *data ;
  float score ;
  float min ;
  float max ;
  float delta ;
} Node ;

float sweep( MuacRecord *in, int nPoints, int reverse, float *score ) ;
void reCalc( Node *node ) ;
Node *addNode( MuacRecord *data, Node *root, Node *treeData, int *treeIndex,
               float *score ) ;
Node *newNode( Node *parent, Node **side, MuacRecord *data, Node *treeData, 
               int *treeIndex, float *score ) ;
Node *allocTree( int nPoints, int *treeIndex ) ;
Node *allocNode( Node *treeData, int *treeIndex ) ;
void freeTree( Node *treeData, int *treeIndex ) ;
void quickSort( MuacRecord *in, int nPoints ) ;
void ranSort( MuacRecord *in, int start, int end ) ;
int ranPartition( MuacRecord *in, int start, int end ) ;
void exchange( MuacRecord *in, int a, int b ) ;
int iRandom( int left, int right ) ;
void printTree( Node *root, FILE *out ) ;
void printNode( Node* node, FILE *out, int depth ) ;
char* padding( int depth ) ;

#define MAX( a, b ) (a > b ? a : b)
#define MIN( a, b ) (a < b ? a : b)


/* --- main routine ----------------------------------- */

/* the main 2-dimensional ks routine, this expects data to be sorted
   by y.  the routine unsortedMuac (below) provides a less demanding
   interface. */

float muac( MuacRecord *in, int n1, int n2 ) {
  
  float maxDelta ;
  float *score ;
  int nPoints ;

  score = malloc( 2 * sizeof( float )) ;
  score[0] = 1.0 / n1 ;
  score[1] = -1.0 / n2 ;
  nPoints = n1 + n2 ;
  maxDelta = sweep( in, nPoints, 0, score ) ;
  maxDelta = MAX( maxDelta, sweep( in, nPoints, 1, score )) ;
  free( score ) ; /* fix for v1.1 */

  return maxDelta ;
}
    

/* --- sort before processing ------------------------- */

/* for most users, this routine is the simplest interface to use,
   whether data are sorted or not.

   a randomized quicksort is used to sort the data by y - it will be
   slower than calling muac directly if you have already sorted data,
   but it is not a significant delay in most cases (in contrast, an
   ordinary quicksort is notoriously slow with sorted data).

   if you are obsessed with speed: call muac directly with sorted
   data.  if you are sure that your data are not pre/almost sorted,
   check out qsort(3) (man 3 qsort) and sort data yourself before
   calling muac. */

float unsortedMuac( float *x1, float *y1, int n1, 
                    float *x2, float *y2, int n2 ) {

  int n ;
  int i ;
  MuacRecord *in ;
  float maxDelta ;

  n = n1 + n2 ;
  in = (MuacRecord*) malloc( n * sizeof( MuacRecord )) ;
  
  for ( i = 0 ; i < n1 ; ++i ) {
    in[i].x = x1[i] ;
    in[i].y = y1[i] ;
    in[i].type = 0 ;
  }
  for ( i = 0 ; i < n2 ; ++i ) {
    in[n1 + i].x = x2[i] ;
    in[n1 + i].y = y2[i] ;
    in[n1 + i].type = 1 ;
  }
  
  quickSort( in, n ) ;
  maxDelta = muac( in, n1, n2 ) ;

  free( in ) ;

  return maxDelta ;
}


/* --- sweep through the data ------------------------- */

/* this is a single pass through the data - either top to bottom or
   bottom to top - calculating the statistic (max cumulative
   difference) for the top or bottom two quadrants, respectively. */

/* re-using the total delta (as rightOrigin) allows us to calculate
   the right hand quadrant from the left (it's a simple offset that is
   non-zero if there are more of one type than the other) - this
   optimization is not in the corresponding python code. */

float sweep( MuacRecord *in, int nPoints, int reverse, float *score ) {

  int i ;
  int iStart ;
  int iEnd ;
  int iStep ;
  float maxDelta = 0.0 ;
  float rightOrigin ;
  Node *root ;
  Node *nextNode ;
  Node *treeData ;
  int treeIndex ;

  treeData = allocTree( nPoints, &treeIndex ) ;

  if ( reverse ) {
    iStart = nPoints - 1 ;
    iEnd = -1 ;
    iStep = -1 ;
  } else {
    iStart = 0 ;
    iEnd = nPoints ;
    iStep = 1 ;
  }

  /* second arg is same as return value (cheat for orphan node) */
  newNode( 0, &root, in + iStart, treeData, &treeIndex, score ) ;
  iStart += iStep ;

  for ( i = iStart ; i != iEnd ; i += iStep ) {
    nextNode = addNode( in + i, root, treeData, &treeIndex, score ) ;
    reCalc( nextNode ) ;

#ifdef DEBUG
    printTree( root, stdout ) ;
#endif

    rightOrigin = root->delta ;
    maxDelta = MAX( maxDelta, MAX( MAX( root->max, root->max - rightOrigin ),
                                   -1 * MIN( root->min, 
                                             root->min - rightOrigin ))) ;
  }

  freeTree( treeData, &treeIndex ) ;

  return maxDelta ;
}


/* --- calculate tree data ---------------------------- */

/* each time we add a node we need to work back up the tree,
   recalculating the values for each sub node.  */

void reCalc( Node *node ) {

  for ( node = node->parent ; node ; node = node->parent ) {
    if ( node->left && node->right ) {
      node->delta = node->left->delta + node->score + node->right->delta ;
      node->min = MIN( MIN( node->left->min, node->left->delta + node->score ),
                       node->left->delta + node->score + node->right->min ) ;
      node->max = MAX( MAX( node->left->max, node->left->delta + node->score ),
                       node->left->delta + node->score + node->right->max ) ;
    } else if ( node->left ) {
      node->delta = node->left->delta + node->score ;
      node->min = MIN( node->left->min, node->left->delta + node->score ) ;
      node->max = MAX( node->left->max, node->left->delta + node->score ) ;
    } else {
      node->delta = node->score + node->right->delta ;
      node->min = MIN( MIN( 0, node->score ),
                       node->score + node->right->min ) ;
      node->max = MAX( MAX( 0, node->score ),
                       node->score + node->right->max ) ;
    }
  }
}


/* --- descend tree, adding new node at bottom -------- */

Node *addNode( MuacRecord *data, Node *root, Node *treeData, int *treeIndex,
               float *score ) {

  Node *parent ;
  Node *node = root ;
  Node **side ;

  while ( node ) {
    parent = node ;
    if ( data->x < node->data->x ) {
      side = &(node->left) ;
      node = node->left ;
    } else {
      side = &(node->right) ;
      node = node->right ;
    }
  }

  return newNode( parent, side, data, treeData, treeIndex, score ) ;
}


Node *newNode( Node *parent, Node **side, MuacRecord *data, Node *treeData, 
               int *treeIndex, float *score ) {

  Node *node ;

  node = allocNode( treeData, treeIndex ) ;
  *side = node ;
  node->parent = parent ;
  node->data = data ;
  node->score = score[data->type] ;
  node->delta = node->score ;
  node->max = MAX( 0, node->delta ) ;
  node->min = MIN( 0, node->delta ) ;
  node->left = NULL ;
  node->right = NULL ;

  return node ;
}


/* --- manage our own chunk of memory ----------------- */

Node *allocTree( int nPoints, int *treeIndex ) {
  *treeIndex = 0 ;
  return (Node*) calloc( nPoints, sizeof( Node ) ) ;
}


Node *allocNode( Node *treeData, int *treeIndex ) {
  return &(treeData[(*treeIndex)++]) ;
}


void freeTree( Node *treeData, int *treeIndex ) {
  *treeIndex = 0 ;
  free( treeData ) ;
}


/* --- randomized quicksort --------------------------- */

/* randomized because i suspect people will use this routine without
   thinking, giving data already sorted, and then be surprised at the
   slow response. */

/* probably not the fastest implementation in the world (especially
   with randomisation), but (a) i don't like pointer arithmetic and
   (b) the ks calculation is (by a constant factor) slower than this
   anyway.  translated to c from cormen et al. */

void quickSort( MuacRecord *in, int nPoints ) {

#ifdef DEBUG
  int i ;
#endif

  ranSort( in, 0, nPoints - 1 ) ;

#ifdef DEBUG
  for ( i = 0 ; i < nPoints ; ++i ) {
    fprintf( stderr, "%d %f\n", i, in[i].y ) ;
  }
#endif
}


/* tail recursion applied - would gcc have optimised this if the while
   was an if and the last statement called ranSort with the right hand
   half of the data? */

void ranSort( MuacRecord *in, int start, int end ) {

  int pivot ;

  while ( start < end ) {
    pivot = ranPartition( in, start, end ) ;
    ranSort( in, start, pivot ) ;
    start = pivot + 1 ;
  }
}


/* this could be improved by using median of 3 (see cormen et al.) */

int ranPartition( MuacRecord *in, int start, int end ) {

  float y ;
  int left, right ;
  int r ;

  r = iRandom( start, end ) ;
  y = in[r].y ;

  left = start - 1 ;
  right = end + 1 ;

  for ( ; ; ) {
    do { --right ; } while( in[right].y > y ) ;
    do { ++left ; } while( in[left].y < y ) ;

    if ( left < right ) {
      exchange( in, left, right ) ;
    } else {
      return right ;
    }
  }
}


void exchange( MuacRecord *in, int a, int b ) {

  float xy ;
  char type ;

  xy = in[a].x ;
  in[a].x = in[b].x ;
  in[b].x = xy ;
  xy = in[a].y ;
  in[a].y = in[b].y ;
  in[b].y = xy ;
  type = in[a].type ;
  in[a].type = in[b].type ;
  in[b].type = type ;
}


int iRandom( int left, int right ) {
  return left + (int) ( ( 1.0 + right - left ) * rand( ) /
                        ( RAND_MAX + 1.0 ) ) ;
}

  
/* --- print tree for debugging ----------------------- */

void printTree( Node *root, FILE *out ) {
  fprintf( out, "\n" ) ;
  printNode( root, out, 0 ) ;
}

void printNode( Node* node, FILE *out, int depth ) {
  fprintf( out, "%sscr: %f min %f, max %f, dlt %f\n", padding(depth), 
           node->score, node->min, node->max, node->delta ) ;
  if ( node->left ) { printNode( node->left, out, depth + 1 ) ; }
  if ( node->right ) { printNode( node->right, out, depth + 1 ) ; }
  fflush( out ) ;
}

char* space = "                                                       " ;

char* padding( int depth ) {
  return ( space + strlen( space ) - depth ) ;
}
