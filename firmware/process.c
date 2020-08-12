/** \file process.c
 *
 * Source file of process module of the NEUROMOTUS project.
 *
 * <pre>
 *  ------------------- Function Description ------------------
 * The process.c file contains all routines needed to properly
 * process the EMG signals acquired from the ADS1299. It includes
 * also helpful routines for communication.
 *
 *  ------------------------- Updates -------------------------
 *
*  2015-10-01 	/ 	Enzo Mastinu 		/ Creation
 * </pre>
 */


/* libraries */
#include "process.h"


/* defines */
// filters parameter
#define NZEROS 			  		2
#define NPOLES 			  		2
#define GAIN_HP_2000	  		1.092933031f
#define GAIN_N_2000            	1.002243999f
#define GAIN_HP_1000			1.045431062f
#define GAIN_N_1000             1.004488020f
#define GAIN_HP_500			  	1.194615832f
#define GAIN_N_500            	1.008976220f


/* structures */


/* prototypes */
float32_t FilterHP(float32_t inData, unsigned char ch, uint32_t frequency);
float32_t FilterNOTCH(float32_t inData, unsigned char ch, uint32_t frequency);
void SWTTimeWindow(float32_t *ptrTWs);
void ISWTTimeWindow(float32_t *ptrTWs);
void DenoiseCoefficients(void);
void DeartifactCoefficients(void);
void circ_conv_f32(float32_t *pSrcA, uint32_t srcALen,
        float32_t *pSrcB, uint32_t srcBLen, float32_t *pDst);


/* global variables */
const float32_t lsbVolt = 0.000000536442f;
float32_t xv[NZEROS+1][8];
float32_t yv[NZEROS+1][8];
float32_t xvN[NPOLES+1][8];
float32_t yvN[NPOLES+1][8];
volatile unsigned char nChannels = 8;
volatile unsigned char filterEnable = 1;
volatile unsigned char compressionEnable = 0;
volatile unsigned char waveletMode;
extern volatile uint32_t sF;
volatile unsigned char batteryPresent = 1;
waveletCoefficients coeffs;


/* functions */

void FilterTimeWindow(float32_t **samples) {
    uint8_t channel;

    for (channel=0; channel<nChannels; channel++) {
        SWTTimeWindow(samples[channel]);
        if(waveletMode) {
            DenoiseCoefficients();
            DeartifactCoefficients();
        }
        ISWTTimeWindow(samples[channel]);
    }
}

/**
 * \brief Needed to change the wavelet parameters when changing sampling
 * frequency or window size.
 * @param[out]  *coeffs structure to store wavelet coefficients in
 * @return none.
 */
void UpdateWaveletCoefficients(void) {
	switch(sF) {
	case 500:
		coeffs.noiseWlevel 	    = 1;
		coeffs.dominantWlevel 	= 0;
	break;
	case 2000:
		coeffs.noiseWlevel		= 0;
		coeffs.dominantWlevel	= 2;
	break;
	default:    // 1000
		coeffs.noiseWlevel 	    = 0;
		coeffs.dominantWlevel 	= 1;
	}
	coeffs.minimaxiFactor = 0.3936 + 0.1829*(log(timeWindowSamples)/log(2));
}

/**
 * \brief Perform denoising on wavelet coefficients
 * @return none.
 *
 * This function should not be called directly. It is a helper function
 * for FilterTimeWindow. It performs signal denoising and Wiener wavelet 
 * filtering on a set of wavelet coefficients.
 */
void DenoiseCoefficients(void) {
	float32_t   sigma, thr, thr2, wcf, var;
	uint32_t    lvl, sampleIdx;
	uint32_t 	windowLen = timeWindowSamples;

	// Use the highest non-dominant frequency band to find STD_DEV and variance
	arm_std_f32(coeffs.cD[coeffs.noiseWlevel], windowLen, &sigma);
	var = sigma*sigma;

	// Calculate noise threshold
	thr = sigma*coeffs.minimaxiFactor;
	thr2 = thr*thr;

	switch(waveletMode) {
	case WAVE_TYPE_SOFT:	// Soft Thresholding
		for(lvl=0;lvl<WAVE_LEVEL;lvl++) {
			for(sampleIdx=0;sampleIdx<windowLen;sampleIdx++) {
				wcf = fabsf(coeffs.cD[lvl][sampleIdx]);
				wcf = (wcf>thr)*(wcf-thr);
				wcf = wcf*wcf;
				coeffs.cD[lvl][sampleIdx] *= wcf/(wcf+var);
			}
		}

		for(sampleIdx=0;sampleIdx<windowLen;sampleIdx++) {
			wcf = fabsf(coeffs.cA[sampleIdx]);
			wcf = (wcf>thr)*(wcf-thr);
			wcf = wcf*wcf;
			coeffs.cA[sampleIdx] *= wcf/(wcf+var);
		}
	break;
	case WAVE_TYPE_HARD: // Hard Thresholding
		for(lvl=0;lvl<WAVE_LEVEL;lvl++) {
			for(sampleIdx=0;sampleIdx<windowLen;sampleIdx++) {
				wcf = fabsf(coeffs.cD[lvl][sampleIdx]);
				wcf = (wcf>thr)*wcf*wcf;
				coeffs.cD[lvl][sampleIdx] *= wcf/(wcf+var);
			}
		}

		for(sampleIdx=0;sampleIdx<windowLen;sampleIdx++) {
			wcf = fabsf(coeffs.cA[sampleIdx]);
			wcf = (wcf>thr)*wcf*wcf;
			coeffs.cA[sampleIdx] *= wcf/(wcf+var);
		}
	break;
	case WAVE_TYPE_HYPER: // Semi-Hyperbolic Shrinkage
		for(lvl=0;lvl<WAVE_LEVEL;lvl++) {
			for(sampleIdx=0;sampleIdx<windowLen;sampleIdx++) {
				wcf = fabsf(coeffs.cD[lvl][sampleIdx]);
				wcf = (wcf>thr)*(wcf*wcf-thr2);
				coeffs.cD[lvl][sampleIdx] *= wcf/(wcf+var);
			}
		}

		for(sampleIdx=0;sampleIdx<windowLen;sampleIdx++) {
			wcf = fabsf(coeffs.cA[sampleIdx]);
			wcf = (wcf>thr)*(wcf*wcf-thr2);
			coeffs.cA[sampleIdx] *= wcf/(wcf+var);
		}
	break;
	case WAVE_TYPE_ADAPT: // Adaptive Denoising Shrinkage
		for(lvl=0;lvl<WAVE_LEVEL;lvl++) {
			for(sampleIdx=0;sampleIdx<windowLen;sampleIdx++) {
				wcf = coeffs.cD[lvl][sampleIdx];
				wcf = wcf - thr + (2*thr)/(1+exp(2.1*wcf/thr));
				wcf = wcf*wcf;
				coeffs.cD[lvl][sampleIdx] *= wcf/(wcf+var);
			}
		}

		for(sampleIdx=0;sampleIdx<windowLen;sampleIdx++) {
			wcf = (coeffs.cA[sampleIdx]);
			wcf = wcf - thr + (2*thr)/(1+exp(2.1*wcf/thr));
			wcf = wcf*wcf;
			coeffs.cA[sampleIdx] *= wcf/(wcf+var);
		}
	break;
	case WAVE_TYPE_NNEG: // Non-Negative
		for(lvl=0;lvl<WAVE_LEVEL;lvl++) {
			for(sampleIdx=0;sampleIdx<windowLen;sampleIdx++) {
				wcf = coeffs.cD[lvl][sampleIdx];
				wcf = (fabsf(wcf)>thr) ? (wcf-(thr2/wcf)) : (0);
				wcf = wcf*wcf;
				coeffs.cD[lvl][sampleIdx] *= wcf/(wcf+var);
			}
		}

		for(sampleIdx=0;sampleIdx<windowLen;sampleIdx++) {
			wcf = (coeffs.cA[sampleIdx]);
			wcf = (fabsf(wcf)>thr) ? (wcf-(thr2/wcf)) : (0);
			wcf = wcf*wcf;
			coeffs.cA[sampleIdx] *= wcf/(wcf+var);;
		}
	break;
	}
}

/**
 * \brief Perform artifact reduction on wavelet coefficients
 * @return none.
 *
 * This function should not be called directly. It is a helper function
 * for FilterTimeWindow. It performs motion artifact reduction on a set
 * of wavelet coefficients.
 */
void DeartifactCoefficients(void) {
	float32_t   sigma, thr_a, mabs;
	uint32_t    lvl, sampleIdx;
	uint32_t 	windowLen = timeWindowSamples;
	uint8_t 	domWlevel = coeffs.dominantWlevel;

	// Use the dominant frequency band to find MABS and STD_DEV
	arm_std_f32(coeffs.cD[domWlevel], windowLen, &sigma);

	mabs = 0;
	for(sampleIdx=0;sampleIdx<windowLen;sampleIdx++) {
		mabs += fabs(coeffs.cD[domWlevel][sampleIdx]);
	}

	mabs /= windowLen;
	thr_a = mabs + sigma;

	for(sampleIdx=0;sampleIdx<windowLen;sampleIdx++) {
		if(fabs(coeffs.cA[sampleIdx]) > thr_a) {
            coeffs.cA[sampleIdx] = 0;

			for(lvl=0;lvl<WAVE_LEVEL;lvl++) {
				if(fabs(coeffs.cD[lvl][sampleIdx]) > mabs) {
					coeffs.cD[lvl][sampleIdx] = 0;
				}
			}
		}
	}
}

/**
 * \brief Performs circular convolution.
 * @param[in] *pSrcA Pointer to signal to convolve
 * @param[in] srcALen Length of signal
 * @param[in] *pSrcB Pointer to filter
 * @param[in] srcBLen Length of filter
 * @param[in] *pDst Pointer to convolution destination. Must be at least srcALen long
 * @return none.
 *
 * This function performs convolution with periodic padding to
 * evaluate the SWT/ISWT correctly.
 */
void circ_conv_f32(float32_t *pSrcA, uint32_t srcALen,
        float32_t *pSrcB, uint32_t srcBLen, float32_t *pDst)
{
    uint32_t    offset = (srcBLen>>1u)-1u;
    uint32_t    scratchLen = srcALen+srcBLen-1u;

    arm_copy_f32(&pSrcA[srcALen-offset], coeffs.s_tmp, offset);
    arm_copy_f32(pSrcA, &coeffs.s_tmp[offset], srcALen);
    arm_copy_f32(pSrcA, &coeffs.s_tmp[srcALen+offset], srcBLen-offset-1);
    // In this function, we lie about the pointer to pDst.
    // Elements [pDst-srcBLen+1 : pDst-1] are unused, so it is safe to do this.
    arm_conv_partial_f32(coeffs.s_tmp, scratchLen,
                         pSrcB, srcBLen, pDst-(srcBLen-1),
                         srcBLen-1, srcALen);
}

/**
 * \brief Performs the SWT on a set of data.
 * @param[in] *ptrSamples data window to transform
 * @param[in] *ptrCoeffs structure to store wavelet coefficients in
 * @return none.
 *
 * This function should not be called directly. It is a helper function
 * for FilterTimeWindow. It is necessary to perform denoising and
 * artifact reduction.
 */
void SWTTimeWindow(float32_t *ptrSamples) {
	float32_t 	hi_d[32] = {-0.4830, 0.8365, -0.2241, -0.1294};
	float32_t 	lo_d[32] = {-0.1294, 0.2241,  0.8365,  0.4830};
	uint32_t    lvl, fIdx;
	uint32_t 	windowLen = timeWindowSamples;
	uint32_t    fSize = FILTER_LEN;

	circ_conv_f32(ptrSamples,windowLen,hi_d,FILTER_LEN,coeffs.cD[0]);
	circ_conv_f32(ptrSamples,windowLen,lo_d,FILTER_LEN,coeffs.cA);

	for(lvl=1;lvl<WAVE_LEVEL;lvl++) {
		// Upsample filters
		for(fIdx=FILTER_LEN-1;fIdx>0;fIdx--) {
			hi_d[fIdx<<lvl] = hi_d[fIdx<<(lvl-1)];
			hi_d[fIdx<<(lvl-1)] = 0;

			lo_d[fIdx<<lvl] = lo_d[fIdx<<(lvl-1)];
			lo_d[fIdx<<(lvl-1)] = 0;
		}

		fSize   =  FILTER_LEN*(1 << lvl);
		circ_conv_f32(coeffs.cA,windowLen,hi_d,fSize,coeffs.cD[lvl]);
		circ_conv_f32(coeffs.cA,windowLen,lo_d,fSize,coeffs.o_tmp);
		arm_copy_f32(coeffs.o_tmp, coeffs.cA, windowLen);
	}
}

/**
 * \brief Performs the ISWT on a set of coefficients
 * @param[in] pDestination pointer to store time-domain data in
 * @param[in] ptrCoeffs structure containing wavelet coefficients
 * @return none.
 *
 * This function should not be called directly. It is a helper function
 * for FilterTimeWindow. It reconstructs a time window from wavelet
 * coefficients.
 */
void ISWTTimeWindow(float32_t *pDestination) {
	float32_t   hi_r[FILTER_LEN] = {-0.1294, -0.2241,  0.8365, -0.4830};
	float32_t   lo_r[FILTER_LEN] = { 0.4830,  0.8365,  0.2241, -0.1294};
	int32_t     lvl; // This has to be signed, or we'll hit an infinite loop!
	uint32_t    dstIdx, srcIdx, first, step, sub;
	uint32_t 	windowLen = timeWindowSamples;

	for(lvl=WAVE_LEVEL-1;lvl>=0;lvl--) {
		step = (1 << lvl);
		// Loop over all possible offsets for given step size
		for(first=0;first<step;first++) {
			arm_fill_f32(0,coeffs.c_tmp,windowLen);
			sub = (step << 1);

			// Even indexes
			for(srcIdx=first,dstIdx=0;srcIdx<windowLen;srcIdx+=sub,dstIdx+=2) {
				coeffs.c_tmp[dstIdx] = coeffs.cA[srcIdx];
			}
			circ_conv_f32(coeffs.c_tmp,dstIdx,lo_r,FILTER_LEN,coeffs.o_tmp);
			arm_copy_f32(coeffs.o_tmp, &coeffs.x_tmp[1], dstIdx-1);
			coeffs.x_tmp[0] = coeffs.o_tmp[dstIdx-1];

			for(srcIdx=first,dstIdx=0;srcIdx<windowLen;srcIdx+=sub,dstIdx+=2) {
				coeffs.c_tmp[dstIdx] = coeffs.cD[lvl][srcIdx];
			}
			circ_conv_f32(coeffs.c_tmp,dstIdx,hi_r,FILTER_LEN,coeffs.o_tmp);
			arm_add_f32(coeffs.o_tmp, &coeffs.x_tmp[1], &coeffs.x_tmp[1], dstIdx-1);
			coeffs.x_tmp[0] += coeffs.o_tmp[dstIdx-1];

			// Odd indexes
			for(srcIdx=first+step,dstIdx=0;srcIdx<windowLen;srcIdx+=sub,dstIdx+=2) {
				coeffs.c_tmp[dstIdx] = coeffs.cA[srcIdx];
			}
			circ_conv_f32(coeffs.c_tmp,dstIdx,lo_r,FILTER_LEN,coeffs.o_tmp);
			arm_add_f32(coeffs.o_tmp, &coeffs.x_tmp[2], &coeffs.x_tmp[2], dstIdx-2);
			coeffs.x_tmp[0] += coeffs.o_tmp[dstIdx-2];
			coeffs.x_tmp[1] += coeffs.o_tmp[dstIdx-1];

			for(srcIdx=first+step,dstIdx=0;srcIdx<windowLen;srcIdx+=sub,dstIdx+=2) {
				coeffs.c_tmp[dstIdx] = coeffs.cD[lvl][srcIdx];
			}
			circ_conv_f32(coeffs.c_tmp,dstIdx,hi_r,FILTER_LEN,coeffs.o_tmp);
			arm_add_f32(coeffs.o_tmp, &coeffs.x_tmp[2], &coeffs.x_tmp[2], dstIdx-2);
			coeffs.x_tmp[0] += coeffs.o_tmp[dstIdx-2];
			coeffs.x_tmp[1] += coeffs.o_tmp[dstIdx-1];

			// Average entries
			for(dstIdx=first,srcIdx=0;dstIdx<windowLen;dstIdx+=step,srcIdx++) {
				coeffs.cA[dstIdx] = 0.5*coeffs.x_tmp[srcIdx];
			}
		}
	}

	// Results stored in coeffs.cA
}

/**
 * \brief It filters a new sample applying a IIR high-pass filter in real time.
 * @param[in] inData float input value
 * @param[in] ch index of the channel
 * @param[in] frequency is the sampling frequency
 * @return float filtered value.
 *
 * Butterworth HP 2nd order @ 20Hz IIR Filter
 * this function takes a float input and applies a real time filter
 * process. The filter needs information about the previous states
 * of input and output values, they are stored into
 * xv[order+1][nChannels] and yv[order+1][nChannels] global arrays.
 * The filters used are IIR:
 * - from current input, past inputs and past outputs
 *   we generate the current value. In FIR you use only the current
 *   and past inputs.
 * IIR parameters are calculated via Matlab:
 * Fc = x;
 * Wn = 20 / (Fc/2);
 * Order = 2;
 * [num den] = butter(Order, Wn, 'high');
 * we normalize the coefficients
 * num = num/num(1)
 * den = den/den(1)
 * then, the coefficients the multiply the inputs are
 * x0 * num(1), x1 * num(2), x2 * num(3)
 * while for the outputs is
 * y2 * den(1), y1 * den(2), y0 * den(3)
 *
 *  Measured Execution Time @80MHz:
 *  13/07/2016 -> 1.31 us
 */
float32_t FilterHP(float32_t inData, unsigned char ch, uint32_t frequency) {

	switch(frequency)
	{
		case 500:
			xv[0][ch] = xv[1][ch]; xv[1][ch] = xv[2][ch];
			xv[2][ch] = inData / GAIN_HP_500;
			yv[0][ch] = yv[1][ch]; yv[1][ch] = yv[2][ch];
			yv[2][ch] =   (xv[0][ch] + xv[2][ch]) - 2 * xv[1][ch] + ( -0.7008967812f * yv[0][ch]) + (  1.6474599811f * yv[1][ch]);
		break;
		case 1000:
			xv[0][ch] = xv[1][ch]; xv[1][ch] = xv[2][ch];
			xv[2][ch] = inData / GAIN_HP_1000;
			yv[0][ch] = yv[1][ch]; yv[1][ch] = yv[2][ch];
			yv[2][ch] =   (xv[0][ch] + xv[2][ch]) - 2 * xv[1][ch] + ( -0.8371816513f * yv[0][ch]) + (  1.8226949252f * yv[1][ch]);
		break;
		case 2000:
			xv[0][ch] = xv[1][ch]; xv[1][ch] = xv[2][ch];
			xv[2][ch] = inData / GAIN_HP_2000;
			yv[0][ch] = yv[1][ch]; yv[1][ch] = yv[2][ch];
			yv[2][ch] =   (xv[0][ch] + xv[2][ch]) - 2 * xv[1][ch] + ( -0.9149758348f * yv[0][ch]) + (  1.9111970674f * yv[1][ch]);
		break;
		default:
			return inData;
	}
	return yv[2][ch];

}

/**
 * \brief It filters a new sample applying a IIR notch filter in real time.
 * @param[in] inData float input value
 * @param[in] ch index of the channel
 * @param[in] frequency is the sampling frequency
 * @return float filtered value.
 *
 * Butterworth NOTCH 1st order @ 45-55 Hz IIR Filter
 * this function takes a float input and applies a real time filter
 * process. The filter needs information about the previous states
 * of input and output values, they are stored into
 * xvN[order+1][nChannels] and yvN[order+1][nChannels] global arrays.
 * IIR parameters are calculated via Matlab:
 * Fc = x;
 * W0 = 50 / (Fc/2);
 * Q = 35;
 * Bw = W0/Q;
 * [num den] = iirnotch(W0, Bw);
 * we normalize the coefficients
 * num = num/num(1)
 * den = den/den(1)
 * then, the coefficients the multiply the inputs are
 * x0 * num(1), x1 * num(2), x2 * num(3)
 * while for the outputs is
 * y2 * den(1), y1 * den(2), y0 * den(3)
 *
 *  Measured Execution Time @80MHz:
 *  13/07/2016 -> 1.38 us
 */
float32_t FilterNOTCH(float32_t inData, unsigned char ch, uint32_t frequency) {

	switch(frequency)
	{
		case 500:
			xvN[0][ch] = xvN[1][ch]; xvN[1][ch] = xvN[2][ch];
			xvN[2][ch] = inData / GAIN_N_500;
			yvN[0][ch] = yvN[1][ch]; yvN[1][ch] = yvN[2][ch];
			yvN[2][ch] = (xvN[0][ch] + xvN[2][ch]) - 1.6180339887f * xvN[1][ch] + ( -0.9822072713f * yvN[0][ch]) + (  1.6036393689f * yvN[1][ch]);
		break;
		case 1000:
			xvN[0][ch] = xvN[1][ch]; xvN[1][ch] = xvN[2][ch];
			xvN[2][ch] = inData / GAIN_N_1000;
			yvN[0][ch] = yvN[1][ch]; yvN[1][ch] = yvN[2][ch];
			yvN[2][ch] = (xvN[0][ch] + xvN[2][ch]) - 1.9021130326f * xvN[1][ch] + ( -0.9910640654f * yvN[0][ch]) + (  1.8936144537f * yvN[1][ch]);
		break;
		case 2000:
			xvN[0][ch] = xvN[1][ch]; xvN[1][ch] = xvN[2][ch];
			xvN[2][ch] = inData / GAIN_N_2000;
			yvN[0][ch] = yvN[1][ch]; yvN[1][ch] = yvN[2][ch];
			yvN[2][ch] = (xvN[0][ch] + xvN[2][ch]) - 1.9753766812f * xvN[1][ch] + ( -0.9955220515f * yvN[0][ch]) + (  1.9709538636f * yvN[1][ch]);
		break;
		default:
			return inData;
	}
	return yvN[2][ch];

}

/**
 * \brief Needed to reset the history of the filters and delete the
 * previous states.
 * @param[in] nCh number of channels
 * @return none.
 *
 *  Measured Execution Time @80MHz:
 *  13/07/2016 -> 16 us
 */
void ResetFilterArrays(unsigned char nCh) {

	unsigned char iCh, k;

	for(iCh=0; iCh<nCh; iCh++) {
		for(k=0; k<3; k++) {
			xv[k][iCh] = 0;
			yv[k][iCh] = 0;
			xvN[k][iCh] = 0;
			yvN[k][iCh] = 0;
		}
	}

}

/**
 * \brief It calls the routines for filtering a new sample.
 * @param[in] inData float input value
 * @param[in] filterIndex index of the channel which the sample belongs
 * @return float output value.
 *
 * This function contains the routines needed to apply digital
 * filters to the data in a on-fly technique. This means that the
 * latest acquired sample is passed through the filters and the
 * filtered version is returned by this function.
 *
 *  Measured Execution Time @80MHz:
 *  13/07/2016 -> 3 us
 */
float32_t FilterSample(float32_t inData, unsigned char filterIndex) {

	float32_t tempFloatValue;
	float32_t tempFloatValue2 = inData;

	if(filterEnable) {
		// filter data
		tempFloatValue = FilterNOTCH(inData, filterIndex, sF);
		tempFloatValue2 = FilterHP(tempFloatValue, filterIndex, sF);
	}
	return tempFloatValue2;

}

/**
 * \brief This function can be used to easily handle the reception on the
 * UART of 4bytes sized variables, e. g. float.
 * @param[in] *ptr points the variable to receive
 * @return none.
 */
bool UARTReceive4Bytes(uint32_t *ptr) {
    uint32_t    bytes[4];
    bool        error = true;

    if (!UARTReceiveByte((uint8_t *) &bytes[0]) && !UARTReceiveByte((uint8_t *) &bytes[1]) &&
            !UARTReceiveByte((uint8_t *) &bytes[2]) && !UARTReceiveByte((uint8_t *) &bytes[3])) {
        *ptr = ((0xFF & bytes[0]) << 24) | ((0xFF & bytes[1]) << 16) |
                ((0xFF & bytes[2]) << 8) | (0xFF & bytes[3]);
        error = false;
    }

    return error;
}

/**
 * \brief This function can be used to easily handle the transmission on the
 * UART of one byte. It is meant to prevent miscommunications from freezing
 * the device.
 * @param[in] *ptr points the variable to receive
 * @return none.
 */
bool UARTReceiveByte(uint8_t *ptr) {
    int32_t     byte;

    bool        error = true;
    uint32_t    lastTick, currentTick, elapsedTicks;
    uint32_t    timeout = COMM_TIMEOUT;

    currentTick = SysTickValueGet();
    lastTick = currentTick;

    while (timeout > 0) {
        byte = UARTCharGetNonBlocking(COMM_UARTPORT);
        if (byte != -1) {
            // Communication successful
            *ptr = (int8_t) byte;
            error = false;
            break;
        }

        elapsedTicks = lastTick-currentTick;
        if(timeout <= elapsedTicks) {
            // Communication timeout
            break;
        } else {
            // Still waiting...
            timeout = timeout - elapsedTicks;
        }

        // Handle update and SysTick looping
        lastTick = currentTick;
        currentTick = SysTickValueGet();
        if (lastTick < currentTick) {
            lastTick = lastTick + TICK_PERIOD;
        }
    }

    return error;
}

/**
 * \brief This function can be used to easily handle the transmit on the
 * UART of 4bytes sized variables, e. g. float.
 * @param[in] *ptr points the variable to send
 * @return none.
 */
void UARTSend4Bytes(unsigned char *ptr) {

	unsigned char k, byte;

	ptr += 3;
	for(k=0; k<sizeof(float32_t); k++) {
		byte = *ptr;
		UARTSendByte(byte);
		ptr--;
	}

}

/**
 * \brief This function can be used to easily handle the transmit on the
 * UART of 4bytes sized variables, e. g. int16.
 * @param[in] *ptr points the variable to send
 * @return none.
 */
void UARTSend2Bytes(unsigned char *ptr) {

	unsigned char k, byte;

	ptr += 1;
	for(k=0; k<sizeof(int16_t); k++) {
		byte = *ptr;
		UARTSendByte(byte);
		ptr--;
	}

}

/**
 * \brief This function can be used to easily handle the transmission on the
 * UART of one byte. It is meant for the application of flow control.
 * @param[in] byte to send
 * @return none.
 */
void UARTSendByte(unsigned char byte) {

	UARTCharPut(COMM_UARTPORT,byte);

}

/**
 * \brief This function is used to handle I2C timeouts if the BATPACKMKII board
 * is not responding or not present.
 */
bool I2CMasterTimeout(uint32_t ui32Base) {
    if(HWREG(ui32Base + I2C_O_MCS) & I2C_MCS_CLKTO) {
        return true;
    } else {
        return false;
    }
}

/**
 * \brief This function is used to read a 16 bit value from a BATPACKMKII
 * battery board via the I2C bus.
 * @param[in] command to send
 * @param[out] *dst points to the variable to recieve
 * @return none.
 */
bool I2CReadBattery(uint8_t cmd, uint8_t *dst) {
    uint8_t     data;

    I2CMasterSlaveAddrSet(I2C1_BASE, BATT_I2C_ADDR, false);
    I2CMasterDataPut(I2C1_BASE, cmd);
    I2CMasterControl(I2C1_BASE, I2C_MASTER_CMD_SINGLE_SEND);
    while(I2CMasterBusy(I2C1_BASE) && !I2CMasterTimeout(I2C1_BASE));
    if (I2CMasterTimeout(I2C1_BASE)) {
        batteryPresent = false;
        return false;
    }
    I2CMasterSlaveAddrSet(I2C1_BASE, BATT_I2C_ADDR, true);
    I2CMasterControl(I2C1_BASE, I2C_MASTER_CMD_SINGLE_RECEIVE);
    while(I2CMasterBusy(I2C1_BASE) && !I2CMasterTimeout(I2C1_BASE));
    if (I2CMasterTimeout(I2C1_BASE)) {
        batteryPresent = false;
        return false;
    }
    data = I2CMasterDataGet(I2C1_BASE);
    *dst = data;

    return true;
}

/*
 * \brief This function is used to read the battery voltage from a BATPACKMKII
 * battery board.
 * @return battery volgate in millivolts, or 0xFFFF if read unsuccessfull
 */
uint16_t BattReadCharge(void) {
    uint8_t     volt_8[2];
    uint16_t    volt_16;

    if (!batteryPresent) {
        return 0xFFFF;
    }

    I2CReadBattery(BATT_VOLT_READ_LOW, &volt_8[0]);
    I2CReadBattery(BATT_VOLT_READ_HIGH, &volt_8[1]);

    volt_16 = volt_8[0] | (volt_8[1] << 8);

    return volt_16;
}
