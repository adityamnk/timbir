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

#ifndef _TIFF_H_
#define _TIFF_H_


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include "typeutil.h"

struct TIFF_img {
  int32_t			height;
  int32_t			width;
  char			TIFF_type;  	/* 'g' = grayscale;               */
					/* 'p' = palette-color;           */
					/* 'c' = RGB full color           */

  uint8_t		**mono;		/* monochrome data, or indices    */
					/* into color-map; indexed as     */
					/* mono[row][col]                 */

  uint8_t		***color;	/* full-color RGB data; indexed   */
					/* as color[plane][row][col],     */
					/* with planes 0, 1, 2 being red, */
					/* green, and blue, respectively  */

  char 			compress_type;	/* 'u' = uncompressed             */

  uint8_t		**cmap;		/* for palette-color images;      */
					/* for writing, this array MUST   */
					/* have been allocated with       */
					/* height=256 and width=3         */
};


/* For the following routines: 1 = error reading file; 0 = success  */

/* This routine allocates space and reads TIFF image */
int32_t read_TIFF ( FILE *fp, struct TIFF_img *img ); 

/* This routine writes out a valid TIFF image */
int32_t write_TIFF ( FILE *fp, struct TIFF_img *img );

/* This routine allocates a TIFF image.     */
/* height, width, TIFF_type must be defined */
int32_t get_TIFF ( struct TIFF_img *img, int32_t height, 
               int32_t width, char TIFF_type );

/* This routine frees memory allocated for TIFF image */
void free_TIFF ( struct TIFF_img *img );

#endif /* _TIFF_H_ */

