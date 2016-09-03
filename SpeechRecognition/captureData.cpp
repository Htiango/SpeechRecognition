# include "captureData.h"

/**
 smoothed signal level; initial value typically set to energy of first frame
 **/
double level = 0;

/**
 background signal level; initially set to avg energy of first 10 frames
 **/
double backgroundOf10 = 0;

int frameCount = 0;

int continueSilenceTime = 0;

int totalFrames = 0;

int outTime  = 0;


typedef struct
{
    int          frameIndex;  /* Index into sample array. */
    int          maxFrameIndex;
    SAMPLE      *recordedSamples;
}
paTestData;


/**
 *  compute the energy of each frame, expressed in decibel
 *
 *  @param audioframe One frame.
 *
 *  @return decibel
 */
double EnergyPerSampleInDecibel(SAMPLE *audioframe, long framesToCalc)
{
    double sum = 0;
    double decibel = 0;
    double pointValue = 0; // the value of each point of the frame
    
    for (int i = 0; i < framesToCalc; i++) {
        pointValue = *audioframe;
        sum += pow(pointValue, 2);
        audioframe++;
    }
    
    decibel = 10 * log10(sum);
    
    return decibel;
}

/**
 *  determine whether the frame is a speech
 *
 *  @param audioframe   input voice frame
 *  @param framesToCalc the length of the frame
 *
 *  @return a bool whether it is a frame
 */
bool classifyFrame(SAMPLE *audioframe, long framesToCalc)
{
    bool isSpeech = false;
    double current = EnergyPerSampleInDecibel(audioframe, framesToCalc);
    double background = backgroundOf10 / 10;
    
    level = ((level * FORGET_FACTOR) + current) / (FORGET_FACTOR + 1);
    
    if(current < background){
        background = current;
    }
    else{
        background = (current - background) * ADJUSTMENT;
    }
    
    if (level < background) {
        level = background;
    }
    if (level - background > THRESHOLD) {
        isSpeech = true;
    }
    
    if (!isSpeech) {
        continueSilenceTime += 1;
//        printf("The Continue Silence Frame is %d\n", continueSilenceTime);
    }
    
    if (isSpeech && continueSilenceTime != 0) {
        printf("The Continue Silence Frame is %d\n", continueSilenceTime);
        continueSilenceTime = 0;
    }
    
    return isSpeech;
    
}


/**
 *  use portaudio to record the input voice data, it can automatically define the end of the speech and end it. using the endpointing algorithm
 *  This routine will be called by the PortAudio engine when audio is needed.
 *  It may be called at interrupt level on some machines so don't do anything
 *  that could mess up the system like calling malloc() or free().
 *
 *
 *  @param inputBuffer     inputBuffer description
 *  @param outputBuffer    outputBuffer description
 *  @param framesPerBuffer the number of frames in each buffer
 *  @param timeInfo        Timing information for the buffers passed to the stream callback.
 *  @param statusFlags     Flag bit constants for the statusFlags to PaStreamCallback.
 *  @param userData        the user's input data
 *
 *  @return     Signal that the stream should stop invoking the callback and finish once all output samples have played.
 */
static int recordCallback( const void *inputBuffer, void *outputBuffer,
                          unsigned long framesPerBuffer,
                          const PaStreamCallbackTimeInfo* timeInfo,
                          PaStreamCallbackFlags statusFlags,
                          void *userData )
{
    paTestData *data = (paTestData*)userData;
    const SAMPLE *rptr = (const SAMPLE*)inputBuffer;
    SAMPLE *wptr = &data->recordedSamples[data->frameIndex * NUM_CHANNELS];
    long framesToCalc;
    long i;
    int finished;
    unsigned long framesLeft = data->maxFrameIndex - data->frameIndex;
    
    (void) outputBuffer; /* Prevent unused variable warnings. */
    (void) timeInfo;
    (void) statusFlags;
    (void) userData;
    
    if( framesLeft < framesPerBuffer )
    {
        framesToCalc = framesLeft;
        finished = paComplete;
    }
    else
    {
        framesToCalc = framesPerBuffer;
        finished = paContinue;
    }
    
    if( inputBuffer == NULL )
    {
        for( i=0; i<framesToCalc; i++ )
        {
            *wptr++ = SAMPLE_SILENCE;  /* left */
            if( NUM_CHANNELS == 2 ) *wptr++ = SAMPLE_SILENCE;  /* right */
        }
    }
    else
    {
        for( i=0; i<framesToCalc; i++ )
        {
            *wptr++ = *rptr++;  /* left */
            if( NUM_CHANNELS == /* DISABLES CODE */ (2) ) *wptr++ = *rptr++;  /* right */
        }
    }
    data->frameIndex += framesToCalc;
    
    if (frameCount < 10) {
        frameCount += 1;
        if (frameCount == 1) {
            level = EnergyPerSampleInDecibel(wptr - framesToCalc, framesToCalc);
            printf("The first energy is %f \n", level);
        }
        else{
            double currentTemp = EnergyPerSampleInDecibel(wptr - framesToCalc, framesToCalc);
            level = ((level * FORGET_FACTOR) + currentTemp) / (FORGET_FACTOR + 1);
        }
        backgroundOf10 += EnergyPerSampleInDecibel(wptr - framesToCalc, framesToCalc);
    }
    else{
        classifyFrame(wptr - framesToCalc, framesToCalc);
    }
    
    if(continueSilenceTime > SILENCETHRESHOLD){
        finished = paComplete;
        cout << "\n Find Endpoint and top! \n";
        totalFrames = data->frameIndex;
    }
    
    
    
    
    return finished;
}


static int playCallback( const void *inputBuffer, void *outputBuffer,
                        unsigned long framesPerBuffer,
                        const PaStreamCallbackTimeInfo* timeInfo,
                        PaStreamCallbackFlags statusFlags,
                        void *userData )
{
    paTestData *data = (paTestData*)userData;
    SAMPLE *rptr = &data->recordedSamples[data->frameIndex * NUM_CHANNELS];
    SAMPLE *wptr = (SAMPLE*)outputBuffer;
    unsigned int i;
    int finished;
    unsigned int framesLeft = totalFrames - data->frameIndex;
    
    (void) inputBuffer; /* Prevent unused variable warnings. */
    (void) timeInfo;
    (void) statusFlags;
    (void) userData;
    
    if( framesLeft < framesPerBuffer )
    {
        /* final buffer... */
        for( i=0; i<framesLeft; i++ )
        {
            *wptr++ = *rptr++;  /* left */
            if( NUM_CHANNELS == /* DISABLES CODE */ (2) ) *wptr++ = *rptr++;  /* right */
        }
        for( ; i<framesPerBuffer; i++ )
        {
            *wptr++ = 0;  /* left */
            if( NUM_CHANNELS == /* DISABLES CODE */ (2) ) *wptr++ = 0;  /* right */
        }
        data->frameIndex += framesLeft;
        finished = paComplete;
    }
    else
    {
        for( i=0; i<framesPerBuffer; i++ )
        {
            *wptr++ = *rptr++;  /* left */
            if( NUM_CHANNELS == /* DISABLES CODE */ (2) ) *wptr++ = *rptr++;  /* right */
        }
        data->frameIndex += framesPerBuffer;
        finished = paContinue;
    }
    return finished;
}

void capture(){
    PaStreamParameters  inputParameters,
    outputParameters;
    PaStream*           stream;
    PaError             err = paNoError;
    paTestData          data;
    int                 i;

    int                 numSamples;
    int                 numBytes;
    SAMPLE              max, val;
    double              average;
    
    printf("patest_record.c\n"); fflush(stdout);
    
    data.maxFrameIndex = totalFrames = NUM_SECONDS * SAMPLE_RATE; /* Record for a few seconds. */
    data.frameIndex = 0;
    numSamples = totalFrames * NUM_CHANNELS;
    numBytes = numSamples * sizeof(SAMPLE);
    data.recordedSamples = (SAMPLE *) malloc( numBytes ); /* From now on, recordedSamples is initialised. */
    if( data.recordedSamples == NULL )
    {
        printf("Could not allocate record array.\n");
        goto done;
    }
    for( i=0; i<numSamples; i++ ) data.recordedSamples[i] = 0;
    
    err = Pa_Initialize();
    if( err != paNoError ) goto done;
    
    inputParameters.device = Pa_GetDefaultInputDevice(); /* default input device */
    if (inputParameters.device == paNoDevice) {
        fprintf(stderr,"Error: No default input device.\n");
        goto done;
    }
    inputParameters.channelCount = 1;                    /* stereo input */
    inputParameters.sampleFormat = PA_SAMPLE_TYPE;
    inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultLowInputLatency;
    inputParameters.hostApiSpecificStreamInfo = NULL;
    
    
    cout << "Press any key to continue." << endl;
    getchar();
    
    /* Record some audio. -------------------------------------------- */
    err = Pa_OpenStream(
                        &stream,
                        &inputParameters,
                        NULL,                  /* &outputParameters, */
                        SAMPLE_RATE,
                        FRAMES_PER_BUFFER,
                        paClipOff,      /* we won't output out of range samples so don't bother clipping them */
                        recordCallback,
                        &data );
    if( err != paNoError ) goto done;
    
    err = Pa_StartStream( stream );
    if( err != paNoError ) goto done;
    printf("\n=== Now recording!! Please speak into the microphone. ===\n"); fflush(stdout);
    
    while( ( err = Pa_IsStreamActive( stream ) ) == 1 )
    {
        Pa_Sleep(1000);
        outTime += 1;
        printf("index = %d, time = %d\n", data.frameIndex, outTime ); fflush(stdout);
    }
    if( err < 0 ) goto done;
    
    err = Pa_CloseStream( stream );
    if( err != paNoError ) goto done;
    
    if(totalFrames != data.maxFrameIndex){
        data.maxFrameIndex = totalFrames;
        numSamples = totalFrames * NUM_CHANNELS;
        numBytes = numSamples * sizeof(SAMPLE);
    }
    
    /* Measure maximum peak amplitude. */
    max = 0;
    average = 0.0;
    for( i=0; i<numSamples; i++ )
    {
        val = data.recordedSamples[i];
        if( val < 0 ) val = -val; /* ABS */
        if( val > max )
        {
            max = val;
        }
        average += val;
    }
    
    average = average / (double)numSamples;
    
    printf("sample max amplitude = " PRINTF_S_FORMAT"\n", max );
    printf("sample average = %lf\n", average );
    
    /* Write recorded data to a file. */
#if WRITE_TO_FILE
    {
        char wavFile[] = "/Users/hty/desktop/record.wav";
        WriteWave(wavFile, data.recordedSamples , numSamples, SAMPLE_RATE);
    }
#endif

    
    /* Playback recorded data.  -------------------------------------------- */
    data.frameIndex = 0;
    
    outputParameters.device = Pa_GetDefaultOutputDevice(); /* default output device */
    if (outputParameters.device == paNoDevice) {
        fprintf(stderr,"Error: No default output device.\n");
        goto done;
    }
    outputParameters.channelCount = 1;                     /* stereo output */
    outputParameters.sampleFormat =  PA_SAMPLE_TYPE;
    outputParameters.suggestedLatency = Pa_GetDeviceInfo( outputParameters.device )->defaultLowOutputLatency;
    outputParameters.hostApiSpecificStreamInfo = NULL;
    
    printf("\n=== Now playing back. ===\n"); fflush(stdout);
    err = Pa_OpenStream(
                        &stream,
                        NULL, /* no input */
                        &outputParameters,
                        SAMPLE_RATE,
                        FRAMES_PER_BUFFER,
                        paClipOff,      /* we won't output out of range samples so don't bother clipping them */
                        playCallback,
                        &data );
    if( err != paNoError ) goto done;
    
    if( stream )
    {
        err = Pa_StartStream( stream );
        if( err != paNoError ) goto done;
        
        printf("Waiting for playback to finish.\n"); fflush(stdout);
        
        while( ( err = Pa_IsStreamActive( stream ) ) == 1 ) Pa_Sleep(100);
        if( err < 0 ) goto done;
        
        err = Pa_CloseStream( stream );
        if( err != paNoError ) goto done;
        
        printf("Done.\n"); fflush(stdout);
    }
    
done:
    Pa_Terminate();
    if( data.recordedSamples )       /* Sure it is NULL or valid. */
        free( data.recordedSamples );
    if( err != paNoError )
    {
        fprintf( stderr, "An error occured while using the portaudio stream\n" );
        fprintf( stderr, "Error number: %d\n", err );
        fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
        err = 1;          /* Always return 0 or 1, but no other return codes. */
    }
    
}

