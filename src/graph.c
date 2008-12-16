/*
 * All malloc/free calls have been changed to the R_alloc
 * construct from the R-system and all calloc/free calls
 * have been changed to the Calloc/Free construct also from
 * the R-system such that R takes over the management of the
 * memory allocation. Robert Castelo (17/11/2008)
*/

/*
 * This file contains the graph handling routines.
 *
 * Copyright (C) 2002 Sampo Niskanen, Patric Östergård.
 * Licensed under the GNU GPL, read the file LICENSE for details.
 */

/*
 * the functions to read and write graph DIMACS files have been removed
 */


#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include <Rinternals.h>
#include <R_ext/RS.h>

#include "graph.h"


/*
 * graph_new()
 *
 * Returns a newly allocated graph with n vertices all with weight 1,
 * and no edges.
 */
graph_t *graph_new(int n) {
	graph_t *g;
	int i;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(n>0);

	/* g=malloc(sizeof(graph_t)); */
	g=(graph_t *) R_alloc(1, sizeof(graph_t));
	g->n=n;
	/* g->edges=malloc(g->n * sizeof(set_t)); */
	g->edges=(set_t *) R_alloc(g->n, sizeof(set_t));
	/* g->weights=malloc(g->n * sizeof(int)); */
	g->weights=(int *) R_alloc(g->n, sizeof(int));
	for (i=0; i < g->n; i++) {
		g->edges[i]=set_new(n);
		g->weights[i]=1;
	}
	return g;
}

/*
 * graph_free()
 *
 * Frees the memory associated with the graph g.
 */
void graph_free(graph_t *g) {
	int i;

	ASSERT((sizeof(setelement)*8)==ELEMENTSIZE);
	ASSERT(g!=NULL);
	ASSERT(g->n > 0);

	for (i=0; i < g->n; i++) {
		set_free(g->edges[i]);
	}
/*
	free(g->weights);
	free(g->edges);
	free(g);
*/
	return;
}


/*
 * graph_resize()
 *
 * Resizes graph g to given size.  If size > g->n, the new vertices are
 * not connected to any others and their weights are set to 1.
 * If size < g->n, the last g->n - size vertices are removed.
 */
void graph_resize(graph_t *g, int size) {
	int i;

	ASSERT(g!=NULL);
	ASSERT(g->n > 0);
	ASSERT(size > 0);

	if (g->n == size)
		return;

	/* Free/alloc extra edge-sets */
	for (i=size; i < g->n; i++)
		set_free(g->edges[i]);
	/* g->edges=realloc(g->edges, size * sizeof(set_t)); */
	g->edges=Realloc(g->edges, size, set_t);
	for (i=g->n; i < size; i++)
		g->edges[i]=set_new(size);

	/* Resize original sets */
	for (i=0; i < MIN(g->n,size); i++) {
		g->edges[i]=set_resize(g->edges[i],size);
	}

	/* Weights */
	/* g->weights=realloc(g->weights,size * sizeof(int)); */
	g->weights=Realloc(g->weights, size, int);
	for (i=g->n; i<size; i++)
		g->weights[i]=1;
	
	g->n=size;
	return;
}

/*
 * graph_crop()
 *
 * Resizes the graph so as to remove all highest-valued isolated vertices.
 */
void graph_crop(graph_t *g) {
	int i;
	
	for (i=g->n-1; i>=1; i--)
		if (set_size(g->edges[i])>0)
			break;
	graph_resize(g,i+1);
	return;
}


/*
 * graph_weighted()
 *
 * Returns TRUE if all vertex weights of graph g are all the same.
 *
 * Note: Does NOT require weights to be 1.
 */
boolean graph_weighted(graph_t *g) {
	int i,w;

	w=g->weights[0];
	for (i=1; i < g->n; i++)
		if (g->weights[i] != w)
			return TRUE;
	return FALSE;
}

/*
 * graph_edge_count()
 *
 * Returns the number of edges in graph g.
 */
int graph_edge_count(graph_t *g) {
	int i;
	int count=0;

	for (i=0; i < g->n; i++) {
		count += set_size(g->edges[i]);
	}
	return count/2;
}
