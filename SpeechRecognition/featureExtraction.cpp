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


double* getFftEnergy(short* frameData, int actualNumPerFrame ){

    double* frameEnergy = new double[actualNumPerFrame];
    Complex *frameDataComplex;
    
    for (int i = 0; i < actualNumPerFrame; i ++) {
        frameDataComplex[i] = frameData[i] * 1.0;
    }
    
    CArray data(frameDataComplex, actualNumPerFrame);
    
    fft(data);
    
    for (int i = 0; i < actualNumPerFrame; i ++) {
        frameEnergy[i] = 10 * log10(pow(abs(data[i]),2));
    }
    return frameEnergy;
}





