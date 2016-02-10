/* ============================================================================
 * Copyright (c) 2013 Charles A. Bouman (Purdue University)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice, this
 * list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * Neither the name of Charles A. Bouman, Purdue
 * University, nor the names of its contributors may be used
 * to endorse or promote products derived from this software without specific
 * prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


#include "allocate.h"

void *get_spc(size_t num, size_t size)
{
	void *pt;

/*	fprintf(stdout, "==> calloc(), num=%zu, size=%zu\n", num, size);*/
	if( (pt=calloc((size_t)num,size)) == NULL ) {
		fprintf(stderr, "==> calloc() error, num=%zu, size=%zu\n", num, size);
		exit(-1);
		}
	return(pt);
}

void *mget_spc(size_t num,size_t size)
{
	void *pt;

	if( (pt=malloc((size_t)(num*size))) == NULL ) {
		fprintf(stderr, "==> malloc() error, num=%zu, size=%zu\n", num, size);
		exit(-1);
		}
	return(pt);
}

void **get_img(size_t wd,size_t ht,size_t size)
{
	size_t i;
	void  **ppt;
	char   *pt;

	ppt = (void **)mget_spc(ht,sizeof(void *));
	pt = (char *)mget_spc(wd*ht,size);

	for(i=0; i<ht; i++) ppt[i] = pt + i*wd*size;

	return(ppt);
}

void free_img(void **pt)
{
	free( (void *)pt[0]);
	free( (void *)pt);
}




/* modified from dynamem.c on 4/29/91 C. Bouman                           */
/* Converted to ANSI on 7/13/93 C. Bouman         	                  */
/* Modified for 1-D case on 6/29/95 C. Bouman         	                  */
/* multialloc( s, d,  d1, d2 ....) allocates a d dimensional array, whose */
/* dimensions are stored in a list starting at d1. Each array element is  */
/* of size s.                                                             */


void *multialloc(size_t s, size_t d, ...)
{
        va_list ap;             /* varargs list traverser */
        size_t max,                /* size of array to be declared */
        *q;                     /* pointer to dimension list */
        char **r,               /* pointer to beginning of the array of the
                                 * pointers for a dimension */
        **s1, *t, *tree;        /* base pointer to beginning of first array */
        size_t i, j;               /* loop counters */
        size_t *d1;                /* dimension list */

        va_start(ap,d);
        d1 = (size_t *) mget_spc(d,sizeof(size_t));

        for(i=0;i<d;i++)
          d1[i] = va_arg(ap,size_t);


        /* Take care of 1-D case separately (6/29/95) */
        if( d==1 ) {
          tree = (char *)mget_spc(d1[0],s*sizeof(char));
          free((void *)d1);
          return((void *)tree);              /* return base pointer */
        }


        r = &tree;
        q = d1;                /* first dimension */
        max = 1;
        for (i = 0; i < d - 1; i++, q++) {      /* for each of the dimensions
                                                 * but the last */
          max *= (*q);
          r[0]=(char *)mget_spc(max,sizeof(char **));
          r = (char **) r[0];     /* step through to beginning of next
                                   * dimension array */
        }
        max *= (size_t)s * (*q);        /* grab actual array memory */
        r[0] = (char *)mget_spc(max,sizeof(char));

        /*
         * r is now set to posize_t to the beginning of each array so that we can
         * use it to scan down each array rather than having to go across and
         * then down 
         */
        r = (char **) tree;     /* back to the beginning of list of arrays */
        q = d1;                 /* back to the first dimension */
        max = 1;
        for (i = 0; i < d - 2; i++, q++) {      /* we deal with the last
                                                 * array of pointers later on */
          max *= (*q);    /* number of elements in this dimension */
          for (j=1, s1=r+1, t=r[0]; j<max; j++) { /* scans down array for
                                                   * first and subsequent
                                                   * elements */

          /*  modify each of the pointers so that it points to
           * the correct position (sub-array) of the next
           * dimension array. s1 is the current position in the
           * current array. t is the current position in the
           * next array. t is incremented before s1 is, but it
           * starts off one behind. *(q+1) is the dimension of
           * the next array. */

            *s1 = (t += sizeof (char **) * *(q + 1));
            s1++;
          }
          r = (char **) r[0];     /* step through to begining of next
                                   * dimension array */
        }
        max *= (*q);              /* max is total number of elements in the
                                   * last pointer array */

        /* same as previous loop, but different size factor */
        for (j = 1, s1 = r + 1, t = r[0]; j < max; j++) 
          *s1++ = (t += s * *(q + 1));

        va_end(ap);
        free((void *)d1);
        return((void *)tree);              /* return base pointer */
}



/*
 * multifree releases all memory that we have already declared analogous to
 * free() when using malloc() 
 */
void multifree(void *r,size_t d)
{
        void **p;
        void *next;
        size_t i;

        for (p = (void **)r, i = 0; i < d; p = (void **) next,i++)
          if (p != NULL) {
            next = *p;
            free((void *)p);
            }
}



