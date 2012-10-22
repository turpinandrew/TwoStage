
/*
** Attempt to model cellular responses to circular stimuli allowing 
** for temporal and spatial characteristics.
** Each ganglion cell is a DoG fed by a grid of receptors.
** Cortical filter is fed by hex array of ganglion cells.
**
** Version 1.0
** In this version, GCells drive everything, pushing there firing to the
** CCells at every time step.
** Version 1.1 Mon Nov  6 06:22:15 EST 2006
**   Has relative refractory period, fixed memory leaks, binary search for cortical firing thresh.
** Version 1.2 Tue Nov  7 19:45:16 EST 2006
**   Read in dysfunction from text files
** Version 1.3 Fri Nov 10 15:34:52 EST 2006
**   Replace setCorticalTriggerPoints() binary search with "average case"
**
** Author: Andrew Turpin and Stuart Gardiner
** Date: Oct 27 2006
*/
#include <float.h>
#define MAXDOUBLE __DBL_MAX__

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "dist.h"
#include "model.h"
#include "gauss.h"

const DogParams DOG_PCELL = {  77.8, 0.07, 0.6, 0.54 };
const DogParams DOG_MCELL = { 115.0, 0.18, 2.0, 1.19 };

//const DogParams DOG_CCELL = {  1.0, 0.25*0.83252/1.0, -0.25, 0.5*0.83252/1.0 };  // spatial freq 1.0
const DogParams DOG_CCELL = {  2.0, 0.25*0.83252/0.5/1.414, 0.5, 0.5*0.83252/0.5/1.414 };  // spatial freq 0.5


double **receptorDistances;  // distance of each receptor from centre of GCell
GCell  **gCellArray;         // Hex array of ganglion cells (independant of stimulus)
double **spikesPerSecond;    // rate of firing for each G cell given a stimulus location
CCell  **cCellArray;         // Hex array of cortical cells
PoissonTable *poissonTable;  // lookup table for poisson noise in CCells
double responseWave[(int)CORTICAL_EPSP_PERIOD]; // lookup table for EPSP response to an incoming spike
double ***EPSPs;             // time series of EPSPs for each cortical cell to allow dynamic threshold setting
double *uzzellTable;         // pre-computed Uzzell weights indexed by (time since last fired)*TIME_UNITS_PER_SECOND
double maxUzzellTime;        // max number of seconds that is contained in uzzellTable
double *cdfSpikesNoStim;      // cdf of trial outputs without a stimulus
double *cdfSpikes;              // cdf of trial outputs with a stimulus
double fosStimuli[NO_FOS_POINTS];  // Allow 1dB steps - increase this number if you ever allow smaller steps!
double fosProbs[NO_FOS_POINTS];
double fosFP;                      // probability of false detection of 100dB stimulus
double threshold;                  // output from FOS curve fitting
double slope;                      // output from FOS curve fitting
double fp;                      // output from FOS curve fitting
double fn;                      // output from FOS curve fitting

/*
** Allocate memory
**
** Return: FALSE if any malloc fails, TRUE otherwise
*/
char
allocateMem() {
        // allocate memory for spikesPerSecond[][]
    if ((spikesPerSecond = (double **)malloc(sizeof(double *) * X)) == NULL)
        return FALSE;
    for(int x = 0 ; x < X ; x++)
        if((spikesPerSecond[x] = (double *)malloc(sizeof(double)*Y)) == NULL)
            return FALSE;

        // allocate memory for receptorDistances[][]
    if((receptorDistances = (double **)malloc(sizeof(double *) * RECEPTOR_X)) == NULL)
            return FALSE;
    for(int i = 0 ; i < RECEPTOR_X ; i++) 
        if((receptorDistances[i] = (double *)malloc(sizeof(double)*RECEPTOR_Y)) == NULL)
            return FALSE;

        // allocate memory for gCellArray[][]
    if((gCellArray = (GCell **)malloc(sizeof(GCell *) * X)) == NULL)
            return FALSE;
    for(int x = 0 ; x < X ; x++)
        if((gCellArray[x] = (GCell *)malloc(sizeof(GCell)*Y)) == NULL)
            return FALSE;

        // allocate memory for cCellArray[][]
    if((cCellArray = (CCell **)malloc(sizeof(CCell *) * X)) == NULL)
            return FALSE;
    for(int x = 0 ; x < X ; x++) 
        if((cCellArray[x] = (CCell *)malloc(sizeof(CCell)*Y)) == NULL)
            return FALSE;

        // allocate memory for ESPSs[][][]
    if((EPSPs = (double ***)malloc(sizeof(double **) * X)) == NULL)
            return FALSE;
    for(int x = 0 ; x < X ; x++){
        if((EPSPs[x] = (double **)malloc(sizeof(double *) * Y)) == NULL)
            return FALSE;
        for(int y = 0 ; y < Y ; y++) 
            if((EPSPs[x][y] = (double *)malloc(sizeof(double) * (int)TRIGGER_END_TIME)) == NULL)
                return FALSE;
    }

    return TRUE;
}//allocateMem()

/*
** Create lookup table for EPSP response to an incoming cortical spike
** Uses t as time in TIME_UNITs
*/
void
createResponseLookup() {
double tFactor;
for(int t = 0 ; t < CORTICAL_EPSP_PERIOD ; t++) {
    tFactor = 2 * M_PI * CORTICAL_CORNER_FREQ * t / TIME_UNITS_PER_SECOND;
    responseWave[t] = tFactor * exp(1-tFactor);
//    printf("t=%2i, Response Lookup %4.2f\n",t,responseWave[t]);
}
return;
}

/*
** Allocate memory and populate table that is indexed by seconds*TIME_UNITS_PER_SECOND
** Also set maxUzzellTime
** uzzellTable[t] = (t/TIME_UNITS_PER_SECOND)^4 / ((t/TIME_UNITS_PER_SECOND)^4 + CORTICAL_REL_REF_PER^4)
*/
void
makeUzzellTable() {
    int n = 32;   // arbirary, can be as small as 1
    uzzellTable = (double *)malloc(sizeof(double) * n);
    int i = 0;

    fprintf(stderr,"Initialising uzzellTable ");
    const double t_rel_to_4 = CORTICAL_REL_REF_PER*CORTICAL_REL_REF_PER*CORTICAL_REL_REF_PER*CORTICAL_REL_REF_PER;
    double time   = 0.0;
    double u      = 0.0;
    do {
            // check that uzzellTable[] is big enough
        if (i == n) {
            n *= 2;
            if ((uzzellTable = (double *)realloc(uzzellTable, n*sizeof(double))) == NULL) {
                fprintf(stderr,"Out of memory for Uzzell table!!!\n");
                exit(-1);
            }
            fprintf(stderr,"dbl ");fflush(stderr);
        }
        uzzellTable[i++] = u;
        time += 1.0/TIME_UNITS_PER_SECOND;
        double temp = pow(time, 4.0);
        u = temp / (temp + t_rel_to_4);
    } while (u < 1.0 - 0.000000001);

    maxUzzellTime = time - 1.0/TIME_UNITS_PER_SECOND;
    fprintf(stderr,"DONE\n");fflush(stderr);
}//makeUzzellTable()

/*
** Read in list of x y type value quadrulples in ascii format
** Skip any lines beginning with #
*/
void
getDysfunction(FILE *f) {
    fprintf(stderr,"Getting dysfunction file...");fflush(stderr);
        // default is no dysfunction
    for(int x = 0 ; x < X ; x++) 
        for(int y = 0 ; y < Y ; y++) {
            gCellArray[x][y].heightScale  = 1.0;
            gCellArray[x][y].spatialScale  = 1.0;
            gCellArray[x][y].fireReduction = 1.0;
        }
        // now read in lines from dysfunction file
    int x, y;
    char type;
    double value;

    #define BUFF_LEN 1024
    char buff[BUFF_LEN];
    while (fgets(buff, BUFF_LEN, f) != NULL) {
        if (buff[0] == '#')
            continue;

        sscanf(buff, "%d %d %c %lf\n",&x,&y,&type,&value);
        if ((x >= 0) && (x < X) && (y >= 0) && (y < Y)) {
            if (type == 'h')
                gCellArray[x][y].heightScale = value;
            else if (type == 's')
                gCellArray[x][y].spatialScale = value;
            else if (type == 'f')
                gCellArray[x][y].fireReduction = value;
            else
                fprintf(stderr,"Unknown dysfunction type %c: ignoring this line\n",type);
        }
        else
                fprintf(stderr,"Dysfunction file: index out of bounds (%d,%d)\n",x,y);
    }

        // Print into output file the mean heightScale, spatialScale, 
        // fireReduction, and %dead with this file.
        // Note: no newline character, so that FOS curve fit results 
        // appear on same line in output file
    double hSum = 0;
    double sSum = 0;
    double fSum = 0;
    double dead = 0;
    for (x = 0 ; x < X ; x++)
       for (y = 0 ; y < Y ; y++) {
           hSum += gCellArray[x][y].heightScale;
           sSum += gCellArray[x][y].spatialScale;
           fSum += gCellArray[x][y].fireReduction;
           if (gCellArray[x][y].fireReduction == 0.0)  
               dead++;
       }

    double divisor = X * Y;
    printf("Dysfunction: h= %4.3f s= %4.3f f= %4.3f d= %4.2f\n",hSum/divisor,sSum/divisor,
                                           fSum/divisor,dead/divisor);

    fprintf(stderr,"DONE\n");fflush(stderr);
}//getDysfunction()

/*
** Initialise distance of each receptor from centre of G cell.
** Initialise each GCell.
** Read in dysfunction from file f.
** Initialise DoG values for each receptor array for each Gcell.
** Initialise every CCell, putting in links from Gcell to Ccell.
*/
void
initialiseCells(FILE *f) {
        // initialise receptor distances from gcell centre array
        // assuming a square grid for now...
    double xSpacing = RECEPTOR_FIELD_EXTENT / (RECEPTOR_X-1);
    double ySpacing = RECEPTOR_FIELD_EXTENT / (RECEPTOR_Y-1);
    double yReceptor, 
           xReceptor = - RECEPTOR_FIELD_EXTENT / 2;
    for(int i = 0 ; i < RECEPTOR_X ; i++) {
        yReceptor = -RECEPTOR_FIELD_EXTENT / 2;
        for(int j = 0 ; j < RECEPTOR_Y ; j++) {
            receptorDistances[i][j] = sqrt(xReceptor*xReceptor + yReceptor*yReceptor);
            yReceptor += ySpacing;
        }
        xReceptor += xSpacing;
    }

        //
        // initialise Gcell array
        //
    getDysfunction(f);
    fprintf(stderr,"Setting up gCells...");fflush(stderr);
    for(int x = 0 ; x < X ; x++) 
        for(int y = 0 ; y < Y ; y++) {
            gCellArray[x][y].refractoryPeriod = RGC_REFRACTORY_PERIOD;
            gCellArray[x][y].numberOfCCells   = 0; 
            double receptorArea = RECEPTOR_FIELD_EXTENT * RECEPTOR_FIELD_EXTENT
                                / (double)(RECEPTOR_X - 1.0) / (double)(RECEPTOR_Y - 1.0);
            for(int i  = 0 ; i < RECEPTOR_X ; i++)
                for(int j  = 0 ; j < RECEPTOR_Y ; j++) {
                    gCellArray[x][y].receptorWeights[i][j] = 
                        gCellArray[x][y].heightScale 
                      * dog(DOG_MCELL, receptorDistances[i][j], gCellArray[x][y].spatialScale) * receptorArea;
                }
        }

        // set GCell hex grid spacing (x,y) + little jitter
    double gX, gY, r;
    double yCell = - floor((double)Y/2.0) * sqrt(3)/2.0 ;
    for(int y = 0 ; y < Y ; y++) {
        double xCell;
        if (y % 2)
            xCell = -floor((double)X/2.0) ;
        else
            xCell = -floor((double)X/2.0) + 0.5;
        for(int x = 0 ; x < X ; x++) {
            for(gX=-2, r=rand()/(double)RAND_MAX ; r > gsl_cdf_gaussian_P (gX, RGC_SPACING_JITTER_X_SD) ; gX += 0.005);
            for(gY=-2, r=rand()/(double)RAND_MAX ; r > gsl_cdf_gaussian_P (gY, RGC_SPACING_JITTER_Y_SD) ; gY += 0.005);
            gCellArray[x][y].x = xCell * GCELL_SPACING + gX;
            gCellArray[x][y].y = yCell * GCELL_SPACING + gY;
            xCell += 1.0;
        }
        yCell += sqrt(3) / 2.0;
    }
    fprintf(stderr,"DONE\n");fflush(stderr);

        //
        // Initialise CCell array and 
        // link each gcell to its ccells by looping over ccells.
        //
    fprintf(stderr,"Setting up cCells...");fflush(stderr);
    for(int x = 0 ; x < X ; x++) {
        for(int y = 0 ; y < Y ; y++) {
            double xCell = gCellArray[x][y].x;  // assume CCell centred on GCell
            double yCell = gCellArray[x][y].y;  // assume CCell centred on GCell

            cCellArray[x][y].timeLastFired = -CORTICAL_ABS_REF_PER; // give em a chance!
            for(int t = 0 ; t < END_TIME ; t++) cCellArray[x][y].sumOfEPSP[t] = 0.0;
                cCellArray[x][y].RGCOut             = 0.0;
            cCellArray[x][y].meanEPSP      = 0.0;  // cause gets accumulated later
        for(int i = 0 ; i < END_TIME ; i++) 
            cCellArray[x][y].spikesList[i] = 0;

                // for each Gcell
            for(int i  = 0 ; i < X ; i++)
                for(int j  = 0 ; j < Y ; j++) {
                    double xDelta = xCell - gCellArray[i][j].x;
                    double yDelta = yCell - gCellArray[i][j].y;
                    double weight = dog(DOG_CCELL, sqrt(xDelta*xDelta + yDelta*yDelta), 1.0);
                    if (fabs(weight) > 0.000001) {
                        int n = gCellArray[i][j].numberOfCCells++;
                        gCellArray[i][j].cCells[n]  = &cCellArray[x][y];
                        gCellArray[i][j].weights[n] = weight; 
                    }
                }
        }
    }
    fprintf(stderr,"DONE\n");fflush(stderr);
}//initialiseCells()

/*
** Go through all ganglion cells and fire them with some prob
** based on how many receptors on under cell.
** If a gcell fires, propergate the spike to the ccells
**
** currentTime: time in seconds
** stimulusOn: boolean saying if stimulus is on or off
*/
int 
checkAllCells(double currentTime, char stimulusOn) {
    int count = 0;
     for(int x = 0 ; x < X ; x++) {
        for(int y = 0 ; y < Y ; y++) {
            double ips = RGC_RESPONSE(0.0) * gCellArray[x][y].fireReduction;
            GCell *g = &gCellArray[x][y];
            if (currentTime - g->timeLastFired >= g->refractoryPeriod) {
                if (stimulusOn)
                    ips = spikesPerSecond[x][y] * UZZELL(currentTime - g->timeLastFired);
                                                        // 13/11: sg added Uzzell for RGCs
                double r=(double)rand()/(double)RAND_MAX;
                if (r < ips/TIME_UNITS_PER_SECOND/(1 - ips * g->refractoryPeriod)) { // sg's magic formula! 
                    gCellArray[x][y].timeLastFired = currentTime;
                    count++;
                    for(int i = 0 ; i < gCellArray[x][y].numberOfCCells ; i++) {
                        CCell *c = gCellArray[x][y].cCells[i];
                                c->RGCOut += THALMIC_CELLS_PER_RGC * gCellArray[x][y].weights[i]; 
                    }
                }
            }
        }
    }
    return count;
} // checkAllCells()

/*
** Add poisson noise to the cortical cells and check for firing.
** Reset sumOfEPSP after checking for next time around.
**
** currentTime: in seconds
** returns: number of spikes
*/
int
cortical(double currentTime) {
    int numSpikes   = 0;
    double spikesIn = 0.0;
    int timePeriod  = (int) round(currentTime * TIME_UNITS_PER_SECOND);
    int stopTime    = (int)(CORTICAL_EPSP_PERIOD);
    if(stopTime + timePeriod > END_TIME) 
        stopTime = (int)END_TIME - timePeriod;

     for(int x = 0 ; x < X ; x++) {
          for(int y = 0 ; y < Y ; y++) {
                //spikesIn = cCellArray[x][y].RGCOut + (double)poissonRand(poissonTable);
                double prand;
                POISSON_RAND(poissonTable, prand);
                spikesIn = cCellArray[x][y].RGCOut + prand;
                cCellArray[x][y].RGCOut = 0.0;
                if(spikesIn > 0) 
                    for (int t = 0 ; t < stopTime ; t++)
                        cCellArray[x][y].sumOfEPSP[timePeriod + t] += spikesIn * responseWave[t];
                if (cCellArray[x][y].sumOfEPSP[timePeriod] * UZZELL(currentTime-cCellArray[x][y].timeLastFired) >= cCellArray[x][y].firingThreshold) {
                    numSpikes++;
                    cCellArray[x][y].fireCount++;
                    cCellArray[x][y].timeLastFired = currentTime;
                    cCellArray[x][y].spikesList[timePeriod] = 1;
              } else 
                    cCellArray[x][y].spikesList[timePeriod] = 0;
//          if((x==9)&&(y==9)) printf("%4i, %5.2f, %3i\n",timePeriod,cCellArray[x][y].sumOfEPSP[timePeriod],cCellArray[x][y].spikesList[timePeriod]);
             cCellArray[x][y].sumOfEPSP[timePeriod] = 0.0;
          }
     }

     return numSpikes;
}//cortical()

/*
** Ignoring any stimulus (ie random firing of gcells), what sumOfEPSP 
** do you get for cortical cells so that there is a base firing rate of 
** CORTICAL_BASE_FIRE_RATE?
** Use this number to set firingThreshold for each cortical cell,
*/
void
setCorticalTriggerPoints() {
      double spikesIn = 0.0;
      int stopTime = (int)(CORTICAL_EPSP_PERIOD);
      
      fprintf(stderr,"Establishing firingThreshold for cCells...");fflush(stderr);
          // First determine the maximum EPSP for each cCell by firing every gCell
     for(int x = 0 ; x < X ; x++) 
          for(int y = 0 ; y < Y ; y++) 
                for(int i = 0 ; i < gCellArray[x][y].numberOfCCells ; i++) {
                     CCell *c = gCellArray[x][y].cCells[i];
                     c->maxSumOfEPSP += THALMIC_CELLS_PER_RGC * gCellArray[x][y].weights[i];
                     c->maxSumOfEPSP += poissonMax(poissonTable);
                }

          // Loop not actually firing any cortical cells, just storing
          // EPSP for determination of firingThreshold.
          // And adding input from other cortical cells (poisson..)
     for(int time = 0 ; time < TRIGGER_END_TIME ; time++) {
          if (time % (int)TIME_UNITS_PER_SECOND == 0)
                fprintf(stderr,"%d ",time/(int)TIME_UNITS_PER_SECOND);fflush(stderr);

             if(stopTime + time > TRIGGER_END_TIME) stopTime = (int)TRIGGER_END_TIME - time;
             else stopTime = (int)(CORTICAL_EPSP_PERIOD);
          (void)checkAllCells((double)time/TIME_UNITS_PER_SECOND, 0);
          for(int x = 0 ; x < X ; x++)
                for(int y = 0 ; y < Y ; y++) {
                    //spikesIn = cCellArray[x][y].RGCOut + (double)poissonRand(poissonTable);
                    double prand;
                    POISSON_RAND(poissonTable, prand);
                    spikesIn = cCellArray[x][y].RGCOut + prand;
                    cCellArray[x][y].RGCOut = 0.0;
                              if(spikesIn > 0) for (int t = 0 ; t < stopTime ; t++)
                                     EPSPs[x][y][time + t] += spikesIn * responseWave[t];
//                     if((x==9)&&(y==9)) printf("SpikesIn=%5.3f; Search: EPSP=%5.3f\n",spikesIn,EPSPs[x][y][time]);
                }
     }

     fprintf(stderr,"bs ");fflush(stderr);
          // Now binary search for a good firingThreshold so that
          // each cell has CORTICAL_BASE_FIRE_RATE on average 
     for(int x = 0 ; x < X ; x++) {
          for(int y = 0 ; y < Y ; y++) {
                double lo = 0;
                double hi = cCellArray[x][y].maxSumOfEPSP; 
                double baseFireRate;
                do {
                     double mid = (hi + lo) / 2.0;
                     baseFireRate = 0.0;
                     cCellArray[x][y].timeLastFired = -CORTICAL_ABS_REF_PER;
                     for(int i = 0 ; i < TRIGGER_END_TIME ; i++) {
                          double currentTime = (double)i/TIME_UNITS_PER_SECOND;
                          if (EPSPs[x][y][i] * UZZELL(currentTime-cCellArray[x][y].timeLastFired) >= mid) {
                                cCellArray[x][y].timeLastFired = currentTime;
                                if (i >= TRIGGER_IGNORE_TIME)
                                     baseFireRate++;
                          }
                     }
                     baseFireRate /= (TRIGGER_END_TIME - TRIGGER_IGNORE_TIME)/TIME_UNITS_PER_SECOND;
                     if (baseFireRate > CORTICAL_BASE_FIRE_RATE) 
                          lo = mid;
                     else
                          hi = mid;
                } while (fabs(hi-lo) > 0.000001);
                cCellArray[x][y].firingThreshold = (hi+lo)/2.0;
             }
     }
          // now reset cortical cells ready for a run!
    //printf("\nFiring thresholds CC\n");
    for(int x = 0 ; x < X ; x++) {
        for(int y = 0 ; y < Y ; y++) {
            cCellArray[x][y].RGCOut = 0.0;
            cCellArray[x][y].timeLastFired = -CORTICAL_ABS_REF_PER;
            //printf("%4.2f ",cCellArray[x][y].firingThreshold);
        }
        //printf("\n");
    }
    fprintf(stderr,"DONE\n");fflush(stderr);
}//setCorticalTriggerPoints()

    // sort highest to lowest (decreasing)
int cmp(const void *a, const void *b) {
    int *aa = (int *)a;
    int *bb = (int *)b;

    return *bb - *aa;
}//cmp()

/*
** Decide if something was seen by finding biggest respose of 
** any cell within a window of width ATTENTION_WINDOW
int
decision(int numSpikes[]) {
    int noWindows = (int) END_TIME - (int) ATTENTION_WINDOW;
    int maxResponse[noWindows];
    for(int i = 0 ; i < noWindows ; i++) 
        maxResponse[i] = 0;
    for(int x = 0 ; x < X ; x++) 
        for(int y = 0 ; y < Y ; y++) {
        int sum = 0;
        for(int i = 0 ; i < ATTENTION_WINDOW ; i++)
            sum += cCellArray[x][y].spikesList[i];
        for(int i = 0 ; i < noWindows ; i++) {
            sum -= cCellArray[x][y].spikesList[i];
            sum += cCellArray[x][y].spikesList[i + (int)ATTENTION_WINDOW];
            if (sum > maxResponse[i]) 
                maxResponse[i] = sum;
        }
    }
    //qsort(maxResponse, noWindows, sizeof(int), cmp);
    //return maxResponse[3];

    int overallMax = 0;
    for(int i = 0 ; i < END_TIME - ATTENTION_WINDOW ; i++) 
        if(maxResponse[i] > overallMax) 
            overallMax = maxResponse[i];
    return overallMax;

}//decision()
*/

/*
** Get the DECISION_RANK'th highest element from a[0..n-1]
** Uses insertion sort but only maintaining the top DECISION_RANK elements
** Returns -1 if r > n.
*/
inline int 
getValueAtRank(int *a, int n) { 

    int sorted[DECISION_RANK+1];
    for(int i = 0 ; i < DECISION_RANK  ; i++) 
        sorted[i] = -1;

    for(int i = 0 ; i < n ; i++) {
        int j;
        for(j = DECISION_RANK; j > 0 && a[i] > sorted[j-1]; j--)
            sorted[j] = sorted[j-1];
        sorted[j] = a[i];
    }

    return sorted[DECISION_RANK-1];
} //getValueAtRank()

/*
** Return the decision variable "Maximum number of spikes by the fourth most
** active cell in any window of width ATTENTION_WINDOW"
**
** Inefficient version Tue 16 Oct 2012 10:49:22 EST
*/
int
decision(int numSpikes[]) {
    int *count = (int *)malloc(sizeof(int) * X * Y);
    if (count == NULL) {
        fprintf(stderr,"Out of memory for counts in decision()\n");
        return 0;
    }
        // get counts for first window
    int index = 0;
    for(int x = 0 ; x < X ; x++) 
        for(int y = 0 ; y < Y ; y++) {
            count[index] = 0;
            for(int j = 0 ; j < (int)ATTENTION_WINDOW ; j++) 
                count[index] += cCellArray[x][y].spikesList[j];
            index++;
        }
    int overallMax = getValueAtRank(count, X*Y);

        // now slide the windows along
    int noWindows = (int) END_TIME - (int) ATTENTION_WINDOW;
    for(int i = 1 ; i < noWindows ; i++) {
        int index = 0;
        for(int x = 0 ; x < X ; x++) 
            for(int y = 0 ; y < Y ; y++) {
                count[index] -= cCellArray[x][y].spikesList[i-1];
                count[index] += cCellArray[x][y].spikesList[i+(int)ATTENTION_WINDOW-1];
                index++;
            }
        int t = getValueAtRank(count, X*Y);
        if (t > overallMax)
            overallMax = t;
    }

    free(count);
    return overallMax;
}//decision()

/*
** Stimulus of intensity centered at (xOffset, yOffset)
** First check which cell's receptors are under the stimulus,
** and store the cell response
**
** Return seen or not seen.
*/
char
doStimulus(double intensity, double xOffset, double yOffset) {
          // For each GCell, add up the receptorWeights that fall under
          // the stimulus and multiply by the stim intensity
          // Store the results in spikesPerSecond[][] for ready access for this stim.
    for(int x = 0 ; x < X ; x++) {
        for(int y = 0 ; y < Y ; y++) {
                // distance of centre of cell from centre of stimulus
            double cellDist =  (gCellArray[x][y].x - xOffset) * (gCellArray[x][y].x - xOffset);
                   cellDist += (gCellArray[x][y].y - yOffset) * (gCellArray[x][y].y - yOffset);
                   cellDist = sqrt(cellDist);
            
            double sum = 0.0;
            for(int i  = 0 ; i < RECEPTOR_X ; i++) {
                for(int j  = 0 ; j < RECEPTOR_Y ; j++) {
                    if (cellDist + receptorDistances[i][j] < STIMULUS_RADIUS) {
                        sum += gCellArray[x][y].receptorWeights[i][j];
                    }
                }
            }
            spikesPerSecond[x][y] = RGC_RESPONSE(intensity * sum) * gCellArray[x][y].fireReduction;
            if (spikesPerSecond[x][y] < 0) 
                spikesPerSecond[x][y] = 0;
            
            gCellArray[x][y].timeLastFired = -RGC_REFRACTORY_PERIOD; // give em a chance!
            
            cCellArray[x][y].timeLastFired = -CORTICAL_ABS_REF_PER;  // give em a chance!
            cCellArray[x][y].fireCount     = 0;                      // ensure CCells all 0
        }
    }

    if (!STOP_AT_RGC) {
            // Now loop for all time steps, generating spikes
        int stimulusOn = 0;
        int numSpikes[(int)END_TIME];
        for(int time = 0 ; time < END_TIME ; time++) {
            if (time == STIM_START)
                stimulusOn = TRUE;
            if (time == STIM_START + STIM_DURATION)
                stimulusOn = FALSE;

            checkAllCells((double)time/TIME_UNITS_PER_SECOND, stimulusOn);

            numSpikes[time] = cortical((double)time/TIME_UNITS_PER_SECOND);
        }

            // return the max number of spikes per window - this is the decision variable
        return decision(numSpikes);
    }
    else
        return 0;
}//doStimulus()

/*
** Get cdf of max number of spikes for a given stimulus intensity
*/
int generateCdf(double intensity, int trials) {
    int *spikes = (int *)malloc(sizeof(int) * trials);
    int maxSpikes = 0;
    for(int i = 0 ; i < trials ; i++) {
        if (i % 10 == 0)
            fprintf(stderr, "%d ",i);
        double x, y, r;
            // jitter x and y for fixation movement
        for(x = -2, r = rand()/(double)RAND_MAX ; r > gsl_cdf_gaussian_P (x, FIXATION_LOSS_X_SD) ; x += 0.005)
            ;
        for(y = -2, r = rand()/(double)RAND_MAX ; r > gsl_cdf_gaussian_P (y, FIXATION_LOSS_Y_SD) ; y += 0.005)
            ;
        spikes[i] = doStimulus(intensity, x, y);
        if(spikes[i] > maxSpikes) 
            maxSpikes = spikes[i];
    }

    for(int i = 0 ; i <= maxSpikes ; i++) 
        cdfSpikes[i] = 0.0;
    for(int i = 0 ; i < trials ; i++)
        for(int j = 0 ; j < spikes[i] ; j++) 
            cdfSpikes[j]++;
    for(int i = 0 ; i <= maxSpikes ; i++) 
        cdfSpikes[i] /= trials;

    // cdfSpikes[i] now contains the proportion of trials at which the maximum spike count
    // is <=i; ie the proportion that would be 'detected' at a threshold of i
    //    ^^^ >=i ???
    // As a check: cdfSpikes[0] should equal 1

    assert(cdfSpikes[0] == 1);

    free(spikes);
    return maxSpikes;
}//generateCdf()

/*
** Generate one point on FOS curve at stim dB, based on responses with no stimulus taken from cdfSpikesNoStim
** AreaUnderROC corresponds with p = P('seen') in ideal 2AFC task
** P('seen') in detection task = 2*(p-0.5)
*/
double generateFOSPoint(double intensity, int maxSpikesNoStim, double * cdfSpikesNoStim) {
    double areaUnderROC = 0;
    int overallMax;

    int maxSpikes = generateCdf(intensity, TRIALS_PER_INTENSITY);
    if(maxSpikesNoStim < maxSpikes) 
        overallMax = maxSpikesNoStim;
    else 
        overallMax = maxSpikes;

    for(int i = 1 ; i <= overallMax ; i++)
        areaUnderROC += 0.5 * (cdfSpikes[i-1] + cdfSpikes[i]) * (cdfSpikesNoStim[i-1] - cdfSpikesNoStim[i]);
    if(!FIT_2AFC) {
        if(areaUnderROC < 0.5) 
            areaUnderROC = 0.5;
        return 2*(areaUnderROC-0.5);      // Note: returning perimetric Prob(detection)
    }
    else 
        return areaUnderROC;           // returning 2AFC Prob(correct)
}//generateFOSPoint

/*
** Calculate the 2-parameter Quick function of a single data point
** i.e. predicted probability of seeing given these parameters
double quickfunction(double contrast, double threshold, double slope) {
    double qf =  1 - pow(2, -pow(contrast / threshold, slope));
    if(FIT_2AFC) qf = 0.5 * (1 + qf);        // Convert into 2AFC probability
    return qf;
}//quickfunction
*/

/*
** Calculate n!/k!(n-k)! for binomial coefficient in pointLikelihood()
double binomialCoef(int n, int k) {
    double Coef = 1.0;
    if(k > (n-k)) k = n-k;    // Make k <= n/2
    for(int i = 1 ; i <= k ; i++) Coef *= (n-k+i)/i;
    return Coef;
}//binomialCoef
*/

/*
** Calculate likelihood of one point when the predicted probability is probPredicted, and the actual is probActual
** assuming binomial sampling on 200 trials per intensity
double pointLikelihood(double probPredicted, double probActual) {
    if(probPredicted > 0.99999) probPredicted = 0.99999;        // Otherwise likelihood function falls apart!
    if(probPredicted < 0.00001) probPredicted = 0.00001;        // Otherwise likelihood function falls apart!

    int actual = (int) round(probActual * TRIALS_PER_INTENSITY); 
    double lHood = binomialCoef(TRIALS_PER_INTENSITY, actual) * pow(probPredicted, actual) * pow(1-probPredicted, TRIALS_PER_INTENSITY-actual);
    return log(lHood);
}
*/

/*
** Calculate loglikelihood of FOS data given the two parameters of the Quick function
double loglikelihood(int noContrasts, double threshold, double slope, int threshEstimate) {
    double sum = 0.0;
    int low = threshEstimate - 8;
    if(low < 0) 
        low = 0;
    int high = threshEstimate + 8;
    if (low > noContrasts) 
        high = noContrasts;
    for(int stim = low ; stim < high ; stim++) 
        if(fosStimuli[stim] != FOS_POINT_NOT_DONE)
            sum += pointLikelihood(quickfunction(CONTRAST(fosStimuli[stim]), CONTRAST(threshold), slope), fosProbs[stim]);
    return sum;
}//loglikelihood
*/

/*
** Search for MLE of threshold and slope of FOS curve
void fitFOS(int noContrasts, int threshEstimate) {
    double slopeTry = 0.0;
    int low = threshEstimate - 4;
    int high = threshEstimate + 4;
    threshold = 0.0;
    slope = 0.0;

    double currentLikelihood = -1000000000;
    double bestLikelihood = -1000000000;
    fprintf(stderr,"threshEstimate=%4i,low=%4i, high=%4i\n",threshEstimate,low,high);
    for(double thresholdTry = low ; thresholdTry <= high ; thresholdTry += 0.1) {
        for(double thresholdTry = low ; thresholdTry <= high ; thresholdTry += 0.1) {
            for(double thresholdTry = low ; thresholdTry <= high ; thresholdTry += 0.1) {
                for(double slopeIndex = -1.0 ; slopeIndex < 2.0 ; slopeIndex += 0.02) {
                    slopeTry = pow(10, slopeIndex);
                    currentLikelihood = loglikelihood(noContrasts, thresholdTry, slopeTry, threshEstimate);
                    if(currentLikelihood > bestLikelihood) {
                        bestLikelihood = currentLikelihood;
                        threshold = thresholdTry;
                        slope = slopeTry;
                    }
                }
            }
        }
    }
}//fitFOS
*/

/*
** Brute force least-squares fit of Cummulative Gaussian
** search is over threshold \in 0:40, slope \in 0:10
**
** SIDEEFECTS: changes global threshold and slope.
*/
void fitFOS(int noContrasts) {
    double min = 1000000000;
    for(int thresholdTry = 0 ; thresholdTry < noContrasts ; thresholdTry++)
        for(int slopeTry = 0 ; slopeTry < 11 ; slopeTry++) {
            for(double fpTry = 0 ; fpTry <= 0.05 ; fpTry += 0.01) {
                for(double fnTry = 0 ; fnTry <= 0.05 ; fnTry += 0.01) {
                    double ls = 0.0;
                    for(int t = 0 ; t < noContrasts ; t++)
                        if (fosStimuli[t] != FOS_POINT_NOT_DONE) {
                            double diff = fosProbs[t] - (fpTry + (1 - fpTry - fnTry)*(1-gsl_cdf_gaussian_P(t - thresholdTry, slopeTry)));
                            ls += diff*diff;
                        }
                    if (ls < min) {
                        threshold = thresholdTry;
                        slope     = slopeTry;
                        fp        = fpTry;
                        fn        = fnTry;
                        min = ls;
                    }
                }
            }
        }
}//fitFOS

/*
** Generate a FOS curve; output mean and standard deviation with Quick function (transformed Weibull)
*/
void
generateFOS() {
    cdfSpikes = (double *)malloc(sizeof(double)*100);
    cdfSpikesNoStim = (double *)malloc(sizeof(double)*100);
    int stimsSoFar = 0;
    double minAccept = 0.1 + 0.45 * FIT_2AFC;
    double maxAccept = 0.9 + 0.05 * FIT_2AFC;
    
    time_t start = time(NULL);
    fprintf(stderr,"Start clock\n");
    fprintf(stderr,"Establishing base firing with no stim...");
    int maxSpikesNoStim = generateCdf(0, TRIALS_NO_STIM);
    fprintf(stderr,"DONE: %ld seconds\n", time(NULL)-start);

    for(int i = 0 ; i <= maxSpikesNoStim ; i++)     // take copy of cdfSpikes
        cdfSpikesNoStim[i] = cdfSpikes[i];

    for(int i = 0 ; i < NO_CONTRASTS ; i++) 
        fosStimuli[i] = FOS_POINT_NOT_DONE;

    #define DO_A_FOS_POINT(_s) do { \
        fosStimuli[_s] = _s; \
        fosProbs[_s] = generateFOSPoint(CONTRAST(_s), maxSpikesNoStim, cdfSpikesNoStim); \
        if ((fosProbs[_s] > minAccept) && (fosProbs[_s] < maxAccept)) \
            stimsSoFar++; \
        fprintf(stderr,"Stimulus %4.1f, Prob of Seeing %5.4f %ld seconds\n", fosStimuli[_s], fosProbs[_s], time(NULL)- start); \
    } while (0);

        // dB values to test regardless of outcome
    #define NUM_INITIAL_CONTRASTS 5
    int initialStim[NUM_INITIAL_CONTRASTS] = { 0, 10, 20, 30, 40};
    for(int i = 0 ; i < NUM_INITIAL_CONTRASTS ; i++)
        DO_A_FOS_POINT(initialStim[i]);

    if (fosProbs[initialStim[NUM_INITIAL_CONTRASTS-1]] > 0.5) {
        fprintf(stdout, "Top dB value not prob <= 0.5\n");
        fprintf(stderr, "Top dB value not prob <= 0.5\n");
        return;
    }
    if (fosProbs[initialStim[0]] < 0.5) {
        fprintf(stdout, "Lowest dB value not prob >= 0.5\n");
        fprintf(stderr, "Lowest dB value not prob >= 0.5\n");
        return;
    }
    int hiStim = 0;
    int loStim = NUM_INITIAL_CONTRASTS-1;
    while(hiStim < NUM_INITIAL_CONTRASTS && fosProbs[initialStim[hiStim]] > 0.5) hiStim++;
    while(loStim >= 0                    && fosProbs[initialStim[loStim]] < 0.5) loStim--;
    if (hiStim == NUM_INITIAL_CONTRASTS) 
        loStim = hiStim = NUM_INITIAL_CONTRASTS - 1;
    if (loStim == -1) 
        loStim = hiStim = 0;
    fprintf(stderr,"hi=%d lo=%d => %d %d\n", hiStim, loStim, initialStim[hiStim],initialStim[loStim]);

    int loDb=-1, hiDb=-1;
    if (hiStim == loStim) {
        if (fosProbs[initialStim[hiStim]] > 0.55)
            loStim++; 
        else if (fosProbs[initialStim[hiStim]] < 0.45) 
            hiStim--;
        else {
            loDb = initialStim[hiStim] - 2;
            hiDb = loDb + 5;
        }
    }
    if (loDb == -1) {  // if haven't set dB values, do interpolation
        if (loStim > NUM_INITIAL_CONTRASTS-1) 
            loStim = NUM_INITIAL_CONTRASTS;
        if (hiStim < 0)
            hiStim = 0;
        loDb = initialStim[hiStim];
        hiDb = initialStim[loStim];
        float mid = (float)loDb + (fosProbs[loDb] - 0.5)/(fosProbs[loDb]-fosProbs[hiDb])*((float)hiDb - (float)loDb);
        loDb = (int)mid - 2;
        hiDb = (int)mid + 2;
    }
/*
        // now search for the 0.5 point, and take +-2 around that
    int loStim = 0;
    int hiStim = NUM_INITIAL_CONTRASTS-1;
    while(fosProbs[initialStim[loStim]] > 0.98 + 0.01*FIT_2AFC) loStim++;
    while(fosProbs[initialStim[hiStim]] < 0.02 + 0.49*FIT_2AFC) hiStim--;
    fprintf(stderr,"lo=%d hi=%d => %d %d\n", loStim, hiStim, initialStim[loStim],initialStim[hiStim]);

    int loDb=-1, hiDb=-1;
    if (loStim == hiStim) {
        if (fosProbs[initialStim[loStim]] > 0.55)
            hiStim++; 
        else if (fosProbs[initialStim[loStim]] < 0.45) 
            loStim--;
        else {
            loDb = initialStim[loStim] - 2;
            hiDb = loDb + 5;
        }
    }
    if (loDb == -1) {  // if haven't set dB values, do interpolation
        if (hiStim > NUM_INITIAL_CONTRASTS-1) 
            hiStim = NUM_INITIAL_CONTRASTS;
        if (loStim < 0)
            loStim = 0;
        loDb = initialStim[loStim];
        hiDb = initialStim[hiStim];
        float mid = (float)loDb + (fosProbs[loDb] - 0.5)/(fosProbs[loDb]-fosProbs[hiDb])*((float)hiDb - (float)loDb);
        loDb = (int)mid - 2;
        hiDb = (int)mid + 2;
    }
*/

    if (loDb < 0)
        loDb = 0;
    if (hiDb > NO_CONTRASTS)
        hiDb = NO_CONTRASTS;

    //fprintf(stderr,"lo=%d hi=%d\n", loDb, hiDb);
    for(int stim = loDb ; stim <= hiDb  ; stim++)
        DO_A_FOS_POINT(stim);
/*
fosStimuli[10] = 10;
fosStimuli[11] = 11;
fosStimuli[12] = 12;
fosProbs[10] = 0.9;
fosProbs[11] = 0.5;
fosProbs[12] = 0.4;
*/

    fitFOS(NO_CONTRASTS);

    printf("Threshold =%4.2f, Spread (sd of CumGauss) =%4.2f, FP=%4.2f FN=%4.2f\n", threshold, slope, fp, fn);

}//generateFOS()

/*
** MAIN!
*/
int
main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr,"Usage: %s filename\n",argv[0]);
        return -1;
    }

    if (!allocateMem()) {
        fprintf(stderr,"Out of memory for something!\n");
        return -1;
    }

    srand(time(NULL));

    poissonTable = poissonInitialise(CORTICAL_NOISE / TIME_UNITS_PER_SECOND);
    createResponseLookup();
    makeUzzellTable();

    FILE *f;
    if (argv[1][0] == '-')
        f = stdin;
    else
        f = fopen(argv[1],"r");
    if (f == NULL) {
        fprintf(stderr,"Cannot open input file %s\n",argv[1]);
        return -1;
    }
    initialiseCells(f);
    fclose(f);

    if (poissonTable == NULL) {
        fprintf(stderr, "Tried to generate a random poisson from empty table\n"); 
        exit(-1);
    }

    for(int fosNumber = 0 ; fosNumber < 1 ; fosNumber++) {
        if (!STOP_AT_RGC) setCorticalTriggerPoints();

    //  printf("# Stimulus %d ",i);
    //  printf(Dysfunction(&d[dys]); 
    //  printf("\n");
    //  double xJitter = (double)rand() / (double)RAND_MAX / 4;
    //  double yJitter = (double)rand() / (double)RAND_MAX / 4;
        //doStimulus(1, xJitter, yJitter);
//    for(int j = 20 ; j < 40 ; j++)
//        for(int i = 0 ; i < NUMBER_OF_STIMULI ; i++)
//            printf("%5i, %5i\n",j,doStimulus((int)CONTRAST(j), 0, 0));

        generateFOS();
    }
}//main
