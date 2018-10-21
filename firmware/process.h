/** \file process.h
 *
 * Header file of process module of the ALC-D project.
 *
 * <pre>
 *  ------------------- Function Description ------------------
 * The process.h file contains all defines, structures and
 * prototypes needed for the process.c file.
 *  ------------------------- Updates -------------------------
 *
*  2015-10-01 	/ 	Enzo Mastinu 		/ Creation
 * </pre>
 */

#ifndef PROCESS_H_
#define PROCESS_H_


/* libraries */
#include <stdint.h>
#include <stdbool.h>
#include "arm_math.h"
#include "driverlib/uart.h"
#include "driverlib/pwm.h"
#include "driverlib/adc.h"
#include "inc/hw_ints.h"
#include "driverlib/interrupt.h"
#include "inc/hw_timer.h"
#include "driverlib/timer.h"
#include "inc/hw_memmap.h"

#include "ADS1299.h"


/* defines */
// Commands' opcode definition for communication protocol
#define TEST_CONNECTION			 'A'
#define DEVNAME_FIRMWAREV_READ	 0xA0
#define ADS1299_REGS_READ		 'S'
#define FILTERS_ENABLE_SET		 'H'
#define WAVELET_ENABLE_SET 		 'E'
#define COMPRESSION_ENABLE_SET   'Q'
#define ADS1299_GAIN_SET		 'K'
#define ADS1299_DATARATE_SET	 'F'
#define ADS1299_LEADOFF_SET      'I'
#define ADS1299_LEADOFF_GET      'g'
#define SAMPLING_FREQ_SET		 'r'
#define TEST_SIGNAL_ENABLE_SET   'P'
#define FEATURES_PARAMETERS_SET	 'f'
#define CONTROL_ALGORITHM_SET	 'a'
#define TEST_FEATURES_EXTRACTION 'x'
#define TEST_CLASSIFIER			 'c'
#define TEST_SIG_CLASSIFIER	     'C'
#define UPDATE_LDA_COEFF		 'U'
#define UPDATE_SVM_COEFF		 'u'
#define BATTERY_CHECK            'B'
#define START_ACQ			 	 'G'
#define STOP_ACQ				 'T'

#define COMM_UARTPORT			  UART0_BASE
#define COMM_UARTINT			  INT_UART0

#define TIME_WINDOW_SAMPLES      256	// 256 ms of time window
#define WINDOW_INCREMENT_SAMPLES 64		// 64 ms of increment
#define OVERLAPPING_SAMPLES      192	// 192 ms of overlapping
#define FEATURES_SIZE 	  		 40		// 8nCh * 5nFeat
#define NO_MOVEMENTS			 7		// 7movs
#define NO_FEATURES			     5		//

// control modes
#define LDA_CONTROL				1
#define SVM_CONTROL				2

// SWT parameters
#define FILTER_LEN 				4
#define WAVE_LEVEL 				4
#define WAVE_OUT_LEN            287 // TIME_WINDOW_SAMPLES + 32 - 1

#define WAVE_TYPE_NONE          0
#define WAVE_TYPE_MOVE          1
#define WAVE_TYPE_SOFT          2
#define WAVE_TYPE_HARD          3
#define WAVE_TYPE_HYPER         4
#define WAVE_TYPE_ADAPT         5
#define WAVE_TYPE_NNEG          6

// BT module
#define WT12_NOTPRESENT         0x00
#define WT12_PRESENT            0x01
#define WT12_CONNECTED          0x02

/* structures */
typedef struct
{
	float32_t samples[8][TIME_WINDOW_SAMPLES];
} timeWindow;
typedef struct
{
	float32_t samples[8][WINDOW_INCREMENT_SAMPLES];
} incrementalWindow;
typedef struct
{
	float32_t tmabs[8];
	float32_t twl[8];
	float32_t tzc[8];
	float32_t tslpch[8];
	float32_t std[8];
	unsigned char featureEnables[5];
} features;
typedef struct
{
	float32_t w[NO_MOVEMENTS][FEATURES_SIZE];
	float32_t bias[NO_MOVEMENTS];
} ldaCoeff;
typedef struct
{
	float32_t w[NO_MOVEMENTS][FEATURES_SIZE];
	float32_t bias[NO_MOVEMENTS];
	float32_t scale[FEATURES_SIZE];
	float32_t shift[FEATURES_SIZE];
} svmCoeff;
typedef struct
{
	float32_t range[FEATURES_SIZE];
	float32_t midrange[FEATURES_SIZE];
} normalizationParameters;
typedef struct
{
    // Approximation and Detail coefficients for SWT
	float32_t 	cA[TIME_WINDOW_SAMPLES];
	float32_t 	cD[WAVE_LEVEL][TIME_WINDOW_SAMPLES];

	uint8_t 	noiseWlevel;
	uint8_t 	dominantWlevel;
	float32_t	minimaxiFactor;

	// Buffers used for SWT/ISWT processing
	float32_t   s_tmp[WAVE_OUT_LEN];
	float32_t 	o_tmp[TIME_WINDOW_SAMPLES];
	float32_t 	c_tmp[TIME_WINDOW_SAMPLES];
	float32_t	x_tmp[TIME_WINDOW_SAMPLES];
} waveletCoefficients;


/* prototypes */
void ResetFilterArrays(unsigned char nCh);
void UpdateWaveletCoefficients(void);
void ExtractFeatures(volatile timeWindow *ptrTWs, features *ptrFeat, float32_t *ptrFeatV, unsigned char nCh);
void UARTSendFeatures(float32_t *ptrFeatV, unsigned char nCh);
void ProcessSaveNewSampleIncrWs(volatile incrementalWindow *ptrIncrWs, unsigned char nCh, unsigned int sampleIndex);
void SlideDownTimeWindow(volatile timeWindow *ptrTWs, unsigned char nCh);
float32_t FilterSample(float32_t inData, unsigned char filterIndex);
uint16_t CompressSample(float32_t inData);
void UARTReceive4Bytes(uint32_t *ptr);
void UARTSend4Bytes(unsigned char *ptr);
void UARTSendByte(unsigned char byte);

unsigned char SVMclassify(unsigned char nCh, unsigned char nMov);
unsigned char DiscriminantTestLDA(unsigned char nCh, unsigned char nMov);
void NormalizeFeatures(unsigned char nFeat, unsigned char nCh);
void IdentifyMovement(unsigned char mode, unsigned char *outIdx, unsigned char *strength, features *ptrFeat);

/* global variables */


#endif /* PROCESS_H_ */

