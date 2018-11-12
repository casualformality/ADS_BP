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
#include "inc/hw_i2c.h"
#include "driverlib/i2c.h"
#include "inc/hw_timer.h"
#include "driverlib/timer.h"
#include "inc/hw_memmap.h"
#include "inc/hw_types.h"

#include "ADS1299.h"


/* defines */
// Commands' opcode definition for communication protocol
#define TEST_CONNECTION			 'A'
#define DEVNAME_FIRMWAREV_READ	 0xA0
#define ADS1299_REGS_READ		 'S'
#define FILTERS_ENABLE_SET		 'H'
#define ADS1299_GAIN_SET		 'K'
#define ADS1299_DATARATE_SET	 'F'
#define ADS1299_LEADOFF_SET      'I'
#define ADS1299_LEADOFF_GET      'g'
#define SAMPLING_FREQ_SET		 'r'
#define BATTERY_CHECK            'B'
#define START_ACQ_COMP           'Q'
#define START_ACQ			 	 'G'
#define STOP_ACQ				 'T'

#define ERROR_BATTERY_LOW        'B'

#define COMM_UARTPORT			  UART0_BASE
#define COMM_UARTINT			  INT_UART0

// BT module
#define WT12_NOTPRESENT         0x00
#define WT12_PRESENT            0x01
#define WT12_CONNECTED          0x02

// BATTPACK module
#define BATT_I2C_ADDR           0x55
#define BATT_SOC_READ_LOW       0x1C
#define BATT_SOC_READ_HIGH      0x1D

/* structures */


/* prototypes */
void ResetFilterArrays(unsigned char nCh);
float32_t FilterSample(float32_t inData, unsigned char filterIndex);
void CompressSamples(float32_t *src, uint8_t len, int16_t *dst);
void UARTReceive4Bytes(uint32_t *ptr);
void UARTSend4Bytes(unsigned char *ptr);
void UARTSend2Bytes(unsigned char *ptr);
void UARTSendByte(unsigned char byte);
bool I2CMasterTimeout(uint32_t ui32Base);
uint16_t BattReadCharge(void);

/* global variables */


#endif /* PROCESS_H_ */

