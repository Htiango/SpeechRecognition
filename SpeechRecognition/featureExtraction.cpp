#include "featureExtraction.h"


/**
 *  do the pre-emphasize   s'[n] = s[n] - alpha * s[n-1]
 *
 *  @param waveData    input original wave data
 *  @param factor      the factor alpha in the formula
 *
 *  @return            waveData after pre-emphasized
 */
short* preemphasized(short* waveData, double factor)
{
    long sampleNum = sizeof(waveData);
    
    short* waveDataAfter = new short[sampleNum];
    
    waveDataAfter[0] = waveData[0];
    
    for (long i = 1; i < sampleNum ; i ++) {
        waveDataAfter[i] = waveData[i] - waveData[i-1] * factor;
    }
    
    return waveDataAfter;
}

/**
 *  construct a hamming window
 *
 *  @param hamWin     hamming window function
 *  @param hamWinSize hamming window size
 */
void hamming( double *hamWin, int hamWinSize )
{
    for (int i = 0; i < hamWinSize; i++)
    {
        hamWin[i] = (double)(0.54 - 0.46 * cos(2 * PI * (double)i / ((double)hamWinSize - 1) ));
    }
}

/**
 *  change the sample data into frame data
 *
 *  @param waveData          sample data
 *  @param numPerframe       the sample point of each frame
 *  @param actualNumPerFrame the actual point after zero-padding
 *
 *  @return  2-d frames, each line represents several points in the frame
 */
short** getFrame(short* waveData, int numPerframe, int actualNumPerFrame)
{
    long sampleNum = sizeof(waveData);
    long frameNum;
    double *hamWin = new double[numPerframe];
    hamming(hamWin, numPerframe);
    
    if( sampleNum % numPerframe == 0){
        frameNum = sampleNum / numPerframe;
    }
    else{
        frameNum = sampleNum / numPerframe +1 ;
    }
    
    short** frameData = new short*[frameNum];
    for (int i = 0; i < frameNum; i++) {
        frameData[frameNum] = new short[actualNumPerFrame];
    }
    
    for (int i = 0; i < frameNum; i ++) {
        for (int j = 0; j < actualNumPerFrame; j++) {
            if ((j >= numPerframe) || ((j + i * numPerframe) >= sampleNum)) {
                frameData[i][j] = 0;
            }
            else frameData[i][j] = waveData[j + i * numPerframe] * hamWin[j];
        }
    }
    return frameData;
}



/**
 *  for each frame, do the fft and compute the energy for the first half
 *
 *  @param frameData         sample point in each frame
 *  @param actualNumPerFrame the length of the frame after zero-padding
 *
 *  @return discrete frequency power spectrum
 */
double* getFftEnergy(short* frameData, int actualNumPerFrame ){

    double* frameEnergy = new double[actualNumPerFrame];
    Complex *frameDataComplex;
    
    for (int i = 0; i < actualNumPerFrame; i ++) {
        frameDataComplex[i] = frameData[i] * 1.0;
    }
    
    CArray data(frameDataComplex, actualNumPerFrame);
    
    fft(data);
    
    for (int i = 0; i < (actualNumPerFrame / 2 + 1) ; i ++) {
        frameEnergy[i] = 10 * log10(pow(abs(data[i]),2));
    }
    return frameEnergy;
}


/**
 *  get the mel-fre power spectrum
 *
 *  @param fftSample discrete frequency power spectrum after fft(the first half)
 *
 *  @return mel-fre power spectrum
 */
double* getMelEnergy(double* fftSampleEnergy){
    int filterNum = MEL_POINT;
    double* melEnergy = new double[filterNum];
    int frameSize = sizeof(fftSampleEnergy);
    int sampleRate = SAMPLE_RATE;
    double freMax = sampleRate / 2 ;
    double melFreMax = 1125 * log(1 + freMax / 700);
    
    double melGap = melFreMax / (filterNum + 1); // mel-fre gap
    
    double* melEdge = new double[filterNum + 2];   // each mel-fre edge
    double* freEdge = new double[filterNum + 2];   // each actual-fre edge
    int* freEdgePoint = new int[filterNum + 2];    // each actual-fre edge point
    
    for (int i = 0; i < filterNum + 2; i++) {
        melEdge[i] = melGap * i;
        freEdge[i] = 700 * (exp( melEdge[i] / 1125) - 1);
        freEdgePoint[i] = round((frameSize + 1) * freEdge[i] / sampleRate * 2);
    }
    
    delete [] melEdge;
    delete [] freEdge;
    
    for (int i = 0; i < filterNum; i++) {
        melEnergy[i] = 0;
    }
    
    for (int i = 1; i <= filterNum; i++) {
        double temp = 0;
        
        for (int j = 0; j < frameSize; j++) {
            if (j < freEdgePoint[i -1]) {
                temp = 0;
            }
            else if (j >= freEdgePoint[i - 1] && j <= freEdgePoint[i]){
                temp = (j - freEdgePoint[i - 1]) / (freEdgePoint[i] - freEdgePoint[i - 1]);
            }
            else if (j >= freEdgePoint[i] && j <= freEdgePoint[i+1]){
                temp = (freEdgePoint[i + 1] - j) / (freEdgePoint[i+1] - freEdgePoint[i]);
            }
            else if (j > freEdgePoint[i + 1]){
                temp = 0;
            }
            melEnergy[i - 1] += fftSampleEnergy[j] * temp ;
        }
        
    }
    
    return melEnergy;
}


/**
 *  do the log for the mel power spectrum
 *
 *  @param melEnergy mel-fre power spectrum
 *
 *  @return log mel-fre power spectrum
 */
double* getMelLogEnergy(double* melEnergy){
    int filterNum = sizeof(melEnergy);
    double* melLogEnergy = new double[filterNum];
    for (int i = 0; i < filterNum; i++) {
        melLogEnergy[i] = log(melEnergy[i]);
    }
    
    return melLogEnergy;
}

/**
 *  do the DCT for each frame.
 *
 *  @param melEnergy        the log mel-fre power spectrum of a frame
 *
 *  @return         the 13-D DCT result
 */
double* getDCT(double* melLogEnergy){
    int filterNum = sizeof(melLogEnergy);
    int dctDimension = DCT_DIMENSION;
    
    double* dctResult = new double[dctDimension];
    
    for (int i = 0; i < dctDimension; i++) {
        for (int j = 0; j < filterNum; j++) {
            dctResult[i] += melLogEnergy[j] * cos(PI * i / (2 * filterNum) *  (2 * j + 1));
        }
    }
    return dctResult;
}

void featureExtraction(){
    short* dataWave;    // store the original wave data
    short* waveDataAfter; // store the wave data after pre-emphasize
    short** frameData;   // store the wave data in each frame
    int* numSample;
    int* sampleRate;
    long frameNum;
    

    
    const char *wavFile = "/Users/hty/desktop/record.wav";
    
    // read in the wave data
    dataWave = ReadWavFile(wavFile, numSample, sampleRate);
    // do the preemphasize
    waveDataAfter = preemphasized(waveDataAfter, PREEMPHASIZED_FACTOR);
    // divide the wave data into frame (Here represent as 2-D)
    frameData = getFrame(waveDataAfter, SAMPLE_PER_FRAME, ACTUAL_SAMPLE_PER_FRAME);
    // get the frame number
    frameNum = sizeof(frameData);
    
    // do the fft and get each frame's power spectrum
    // get the mel-fre power spectrum, do the log and DCT
    double** frameEnergy = new double*[frameNum];
    double** melEnergy = new double*[frameNum];
    double** melLogEnergy = new double*[frameNum];
    double** frameDCT = new double*[frameNum];
    for (long i = 0; i < frameNum; i++) {
        frameEnergy[i] = getFftEnergy(frameData[i], ACTUAL_SAMPLE_PER_FRAME);
        melEnergy[i] = getMelEnergy(frameEnergy[i]);
        melLogEnergy[i] = getMelLogEnergy(melEnergy[i]);
        frameDCT[i] = getDCT(melLogEnergy[i]);
    }
    
    
}

