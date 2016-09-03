#ifndef _READWAV_H_
#include "stdio.h"

struct WavFileHead
{
	//Resource Interchange File Flag (0-3) "RIFF"
	char RIFF[4];
	//File Length ( not include 8 bytes from the beginning ) (4-7)
	long FileLength;
	//WAVE File Flag (8-15) "WAVEfmt "
	char WAVEfmt_[8];
	//Transitory Byte ( normally it is 10H 00H 00H 00H ) (16-19)
	unsigned long noUse;
	//Format Category ( normally it is 1 means PCM-u Law ) (20-21)
	short FormatCategory;
	//NChannels (22-23)
	short NChannels;
	//Sample Rate (24-27)
	long SampleRate;
	//l=NChannels*SampleRate*NBitsPersample/8 (28-31)
	long SampleBytes;
	//i=NChannels*NBitsPersample/8 (32-33)
	short BytesPerSample;
	//NBitsPersample (34-35)
	short NBitsPersample;
	//Data Flag (36-39) "data"
	char data[4];
	//Raw Data File Length (40-43)
	long RawDataFileLength;
};

// original functions
bool	WaveRewind(FILE *wav_file, WavFileHead *wavFileHead);
short	*ReadWave(char *wavFile, int *numSamples, int *sampleRate);
void	WriteWave(char *wavFile, short *waveData, int numSamples, int sampleRate);
void	FillWaveHeader(void *buffer, int raw_wave_len, int sampleRate);

// additive functions
void    GetWavHeader(char *wavFile,short *Bits,int *Rate,short *Format,int *Length,short *Channels); 
short   *ReadWavFile(char *wavFile, int *numSamples, int *sampleRate); 
void    readwav_t(char *wavFile, short *waveData, long times, int *numSamples, int *sampleRate);
void    GetWavTime(char *wavFile, double *duration);
void    ReadWav(char *wavFile, short *waveData, int *numSamples, int *sampleRate);

#endif //_READWAV_H_
