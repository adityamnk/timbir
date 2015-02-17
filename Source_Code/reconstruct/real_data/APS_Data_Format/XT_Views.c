#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include "XT_Constants.h"
#include <stdlib.h>

int32_t gen_interlaced_views_0_to_Inf (float* views, float* times, int32_t K, int32_t N_theta, int32_t proj_start, int32_t N_p, float min_acq_time, float rot_speed, FILE* debug_file_ptr)
{
	int32_t L = N_theta/K;
        int32_t maxbits = floor(log(K)/log(2.0) + 0.5);/*Implements round off function.*/
        float delta_theta = M_PI/N_theta;
	float min_dTheta = rot_speed*min_acq_time;
	int32_t buf1 = 0, buf2 = 0, i, j, i_clip, temp;

	fprintf(debug_file_ptr, "The number of bits required to represent K = %d is %d.\n", K, maxbits);
	float* proj_angles = (float*)calloc (N_p + proj_start, sizeof(float));
	float* proj_times = (float*)calloc (N_p + proj_start, sizeof(float));
	for (i_clip = 0, i = 0; i_clip < proj_start + N_p; i++)
	{
        	buf1 = i*K;
        	temp = i/L % K;
		buf2 = 0;
		for (j = 0; j < maxbits; j++)
    		{
        		if((temp & (1 << j)))
       			   buf2 |= 1 << ((maxbits - 1) - j);  
   		}

        	proj_angles[i_clip] = (buf1 + buf2);
        	proj_times[i_clip] = proj_angles[i_clip]/K;
    	    	proj_angles[i_clip] = proj_angles[i_clip]*delta_theta;
		
		if (i == 0 || proj_angles[i_clip] - proj_angles[i_clip-1] > min_dTheta)
			i_clip++;
	}

	for (i = 0; i < N_p; i++)
	{
		views[i] = proj_angles[proj_start+i];
		times[i] = proj_times[proj_start+i];
	}

	free(proj_angles);
	free(proj_times);
	return (i_clip);
}
