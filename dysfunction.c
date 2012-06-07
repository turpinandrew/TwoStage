/*
** Program to generate input files for model that
** describe various types of dysfuntion in the ganglion cells.
**
** Author: Andrew Turpin and Stuart Gardiner
** Date: Tue Nov  7 2006
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "model.h"
#include "dysfunction.h"
#include "gauss.h"

/*
** Output x y f 0.0 for percent of the cells
**
** percent: floating point string in range 0..100
*/
void
death(char *percent) {
    float p;
    sscanf(percent, "%f", &p);
    printf("# Dysfunction death of %0.4f%%)\n",p);
    printf("# Note indices start at 0, and dimensions of gcell array are %d * %d\n",X,Y);
    p /= 100.0;
    for(int x = 0 ; x < X ; x++)
        for(int y = 0 ; y < Y ; y++)
            if ((double)rand()/(double)RAND_MAX < p)
                printf("%3d %3d f 0.0\n",x,y);
}//death()

/*
** Draw randomly from gaussian distribution for every cell
*/
void
gaussian(float mean, float stdev, char type) {
    printf("# Dysfunction gaussian(%0.4f, %0.4f) of ",mean, stdev);
    switch (type) {
        case DYSFUNCTION_HEIGHT : printf("height\n"); break;
        case DYSFUNCTION_WIDTH  : printf("spatial\n"); break;
        case DYSFUNCTION_FIRE   : printf("fire\n"); break;
    }
    printf("# Note indices start at 0, and dimensions of gcell array are %d * %d\n",X,Y);
    
    for(int x = 0 ; x < X ; x++)
        for(int y = 0 ; y < Y ; y++) {
            float value;
            if (stdev == 0)
                value = mean;
            else {
                    // XXX simple linear search with function calls. Can be sped up heaps!
                value = (double)rand()/(double)RAND_MAX;
                double p;
                for(p = -mean ; p < 1 - mean; p += 0.00001) 
                    if (gsl_cdf_gaussian_P(p, stdev) >= value) 
                        break;
                value = mean + p;
            }
            switch (type) {
                case DYSFUNCTION_HEIGHT : printf("%3d %3d h %f\n",x,y,value); break;
                case DYSFUNCTION_WIDTH  : printf("%3d %3d w %f\n",x,y,value); break;
                case DYSFUNCTION_FIRE   : printf("%3d %3d f %f\n",x,y,value); break;
            }
        }
}//gaussian()

void
usage(char *name) {
    fprintf(stderr,"Usage: %s type params\n",name);
    fprintf(stderr,"   where\n");
    fprintf(stderr,"      DEATH   type = d params = %%-of-cells-affected-[0..100]\n");
    fprintf(stderr,"      HEIGHT  type = h params = mean-gaussian stdev-gaussian\n");
    fprintf(stderr,"      SPATIAL type = s params = mean-gaussian stdev-gaussian\n");
    fprintf(stderr,"      FIRE    type = f params = mean-gaussian stdev-gaussian\n");
}

int
main(int argc, char*argv[]) {
    srand(time(NULL));

    if ((argc == 3)  && (argv[1][0] == 'd')) 
        death(argv[2]);
    else if (argc == 4) {
        float mean, stdev;
        sscanf(argv[2], "%f", &mean);
        sscanf(argv[3], "%f", &stdev);
        if (argv[1][0] == 'h') 
            gaussian(mean, stdev, DYSFUNCTION_HEIGHT);
        else if (argv[1][0] == 's') 
            gaussian(mean, stdev, DYSFUNCTION_WIDTH);
        else if (argv[1][0] == 'f') 
            gaussian(mean, stdev, DYSFUNCTION_FIRE);
        else {
            usage(argv[0]);
            return -1;
        }
    } else {
        usage(argv[0]);
        return -1;
    }

    return 0;
}//main()
