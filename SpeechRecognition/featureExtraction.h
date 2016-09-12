#ifndef _FEATUREEXTRACTION_H_
#define _FEATUREEXTRACTION_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include "readwave.h"
#include "fft.h"

using namespace std;


//----------
#define SAMPLE_RATE             (44100)  // The sample rate that mac allows only

#define PREEMPHASIZED_FACTOR    (0.95)   // s[n] = s[n] - alpha * s[n-1]

#define SAMPLE_PER_FRAME        (882)    // in order to make each frame 20ms

#define FRAME_INTERVAL_TIME     (10)      // set the interval time of each frame to be 10ms, so have 10ms overlapping

#define ACTUAL_SAMPLE_PER_FRAME (1024)    // in order to do the fft

#define MEL_POINT               (40)

#define DCT_DIMENSION           (13)

#define fre_MinMel              (133.33)

#define fre_MaxMel              (6855.4976)

// ---------
#define NUM_FILTER              (40)     // set the num of filters to be 40

void featureExtraction();

#endif
