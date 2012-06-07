#define TRUE 1
#define FALSE 0

    // size of the ganglion cell array (must be odd number)
#define X 19
#define Y 19

#define RECEPTOR_X              101     // number of boundary points (not cells) (must be odd)
#define RECEPTOR_Y              101     // number of boundary points (not cells) (must be odd)
#define RECEPTOR_FIELD_EXTENT   (6.0 * DOG_MCELL.sWidth)    // in degrees
#define GCELL_SPACING           0.16    // MCELLs degrees
#define STIMULUS_RADIUS         0.43    // Size III in degrees

    // be careful not to end up with non integer END_TIME etc... 
#define TIME_UNITS_PER_SECOND  2000.0
#define TRIGGER_IGNORE_TIME (1.0000 * TIME_UNITS_PER_SECOND) // TIME UNITS - ignore for spikes
#define TRIGGER_END_TIME    (11.000 * TIME_UNITS_PER_SECOND) // TIME UNITS - training period
#define END_TIME            (1.0000 * TIME_UNITS_PER_SECOND) // TIME UNITS - must be whole integer
#define STIM_START          (0.4000 * TIME_UNITS_PER_SECOND) // TIME UNITS - must be whole integer

#define STIM_DURATION       (0.2000 * TIME_UNITS_PER_SECOND) // TIME UNITS - must be whole integer
#define ATTENTION_WINDOW    (int) (0.2000 * TIME_UNITS_PER_SECOND) // TIME UNITS - must be whole integer

#define TRIALS_PER_INTENSITY	40
#define TRIALS_NO_STIM		200
#define NO_FOS_POINTS	    	41
#define FIT_2AFC				TRUE

typedef struct corticalCell {
    double timeLastFired;		// in number of seconds from start
    double firingThreshold;		// firing threshold set in code...
    double sumOfEPSP[(int)END_TIME];// accumulates EPSPs for firing.
    double RGCOut;			// sum of instantaneous inputs from RGCs
    double meanEPSP;			// for 1/CORTICAL_BASE_FIRE_RATE seconds
    double maxSumOfEPSP;		// delete me !!!
    int spikesList[(int)END_TIME];	// did cell fire at that timepoint?
    int fireCount;
} CCell;

typedef struct gcell {
    double refractoryPeriod;    // number of seconds before can fire again
    double timeLastFired;       // in number of seconds from start

    double heightScale;         // type of dysfunction, 1.0 for normal
    double spatialScale;        // type of dysfunction, 1.0 for normal
    double fireReduction;       // type of dysfunction, 1.0 for normal, multiper on gcell fire rate after everything 
    double receptorWeights[RECEPTOR_X][RECEPTOR_Y];

    float x, y;                 // centre of cell (hex grid, so not obvious) 

    int numberOfCCells;         // number of CCells that link to this gcell
    CCell  *cCells[X*Y];        // array of pointers to ccells that take input from this gcell
    double  weights[X*Y];       // DoG weight to pass up to CCell[][]
                                // XXX - over allocation of memory for now.
} GCell;

#define RGC_REFRACTORY_PERIOD 0.0015 // seconds
#define RGC_BASE_FIRE_RATE   10.0    // "background noise" in fires per second 
#define RGC_MAX_FIRE_RATE    65.0    // Maximum fires per second
#define RGC_SEMI_SATURATION  100.0   // Bill Swanson's magic number
    // Bill's magic formula for non-linear response from GCell that asymptotes to RGC_MAX_FIRE_RATE
#define RGC_RESPONSE(x)    (RGC_BASE_FIRE_RATE + ((RGC_MAX_FIRE_RATE-RGC_BASE_FIRE_RATE)*(x)/(RGC_SEMI_SATURATION + (x)))) 

#define CORTICAL_BASE_FIRE_RATE 5.0     // "background noise" firing rate
#define CORTICAL_NOISE   (250.0 * CORTICAL_BASE_FIRE_RATE)   // using this as mean of poisson
                                                             // 250 cells at 5 ips feeding into each Ccell
#define CORTICAL_ABS_REF_PER 0.0015         // absolute refractory period
#define CORTICAL_REL_REF_PER 0.0020         // relative refractory period
#define CORTICAL_CORNER_FREQ 55.0		  // "corner frequency" for EPSP response to an incoming spike
#define CORTICAL_EPSP_PERIOD (0.0400 * TIME_UNITS_PER_SECOND) // time period of EPSP to sum over; in time periods

#define THALMIC_CELLS_PER_RGC 4.0  // this is a multiplier on GCell input to CCells
                                   // EPSP of gcell c.f. cortical 

    // assume that uzzellTable only goes up to a point where the value is close enough to 1.0
#define UZZELL(t) ( \
    ((t) <= CORTICAL_ABS_REF_PER) ? \
        0.0 : \
        ((t) > maxUzzellTime ? 1.0 : uzzellTable[(int)((double)(t) * (double)TIME_UNITS_PER_SECOND)])  \
)

#define CONTRAST(dB)    ((100.0 * 10000.0) / (31.5 * pow(10.0, ((dB)/10.0)) ))
#define DB(contrast)    ((10.0 * log10(10000.0 * 100.0 / (contrast) / 31.5 )))

//#define STOP_AT_RGC    TRUE
#define STOP_AT_RGC    FALSE

#define FIXATION_LOSS_X_SD 0.3
#define FIXATION_LOSS_Y_SD 0.2

    // The centre of the RGCs is offset from a regular hex grid 
    // jittered by a gaussian with the following standard dev
#define RGC_SPACING_JITTER_X_SD 0 // 0.02
#define RGC_SPACING_JITTER_Y_SD 0 // 0.02
