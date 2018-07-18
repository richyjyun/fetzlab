/* nc323save.c	
 * Copyright (C) 2016 by the University of Washington Regional Primate Center.
 * All rights reserved.
 * Writty by Larry Shupe
 * October 26, 2016
 *
 * Reads NC3B data stored on an SDCARD using nc320.m interface.
 *
 * The calling syntax is:
 *
 *	errcode = nc323save(Drive, Sector, Duration, Pathname, Samplerates, button_handle, 'String', '00:00:00');
 *
 *      Drive -- SDCard file system drive letter (e.g. 'E').
 *      Sector -- Starting sector (512 byte sectors).
 *      Duration -- maximum time in Seconds. (0 for upto end of file.)
 *      Pathname -- String containing file system path name for output files.
 *      Samplerates -- Sample rate of each channel  (Samplerates[0] => Ch01)
 *      button_handle -- Matlab UI object handle to an GUI object to receive time updates.
 *
 *  Creates output files for events, channels, .
 *      <Pathname>_Chan01.i16   to   <Pathname>_Chan32.i16 -- Analog channels
 *      <Pathname>_Chan33.u16   to   <Pathname>_Chan35.u16 -- Aux channels
 *      <Pathname>_AccelX.i16   to   <Pathname>_AccelZ.i16 -- Accelerometer XYZ axis channels
 *      <Pathname>_AccelM.i16   and  <Pathname>_AccelT.i16 -- Accelerometer Magnitude and Temperature channels
 *      <Pathname>_Digi00.u16  -- Contains header word for each sample packet
 *      <Pathname>_Events.u32  -- Event files contain (ID, Value, Timestamp) triplets.
 *              ID = uint32, Value = int32, Timestamp = uint32 at 5000 sample rate)
 *
 * v308 -- support for the 32 channel Intan chip.  Aux channels moved to Channes 33rd,34th,35th
 * v320 -- support for new stim_type/condition_id method of storing stimulus events
 *          and condition switches.  Accelerometer support.
 * V323 -- Matches GUI version nc323.m.  End Condition is condition id 9.
 ***************************************************************************/

#include <windows.h>
#include "mex.h"
#include "matrix.h"
#include "math.h"
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <io.h>
#include <stdio.h>
#include <string.h>

/* Input Arguments */

#define	DRIVE_ARG	prhs[0]
#define	SECTOR_ARG	prhs[1]
#define	TIME_ARG	prhs[2]
#define	PATH_ARG	prhs[3]
#define	RATES_ARG	prhs[4]
#define	SETP1_ARG prhs[5]
#define	SETP2_ARG prhs[6]
#define	SETP3_ARG prhs[7]

/* Output Arguments */

#define	ERRCODE	plhs[0]

/* Constants */

// Size of our read buffer in sectors.
#define INBUF_SECTORS (500)

// Size of file output buffers.
// This needs to be large enough to cover a practical limit on the number
// of events and data values per second.  20004 covers 1 second of
// 20kHz data and is the minimum value.
#define OUTBUF_SIZE (21000)
#define OUTBUF_ACCEL_SIZE (110)

// Maximum number of channels
#define MAX_CHANNELS (35)

// Number of Accelerometer channels (Magnitude,X,Y,Z,Temperature)
#define MAX_ACCEL_CHANNELS (5)

// Number of low rate data channels
#define MAX_LOWRATE_CHANNELS (48)

/* Prototypes */

int makeFile(char *pathName, char *suffix)
{
    int fdout;
    char filename[2048];
    char *extension;
    
    strncpy(filename, pathName, 2048);
    extension = strrchr(filename, '.'); 
    if (extension) // Remove extension
        *extension = '\0';
    strcat(filename, suffix);
	fdout = _open(filename, _O_CREAT | _O_TRUNC | _O_BINARY | _O_RDWR, _S_IREAD | _S_IWRITE);
    if (fdout == -1) 
		mexPrintf("Error: nc3save.c could not create file: %s\n", filename);
    return fdout;
}

void errlog(char *text)
{
	FILE *fd = fopen("c:\\\\nc3\\nc3save_errorlog.txt", "a+");
	fprintf(fd, "%s\n", text);
	fclose(fd);
}

#define POLY 0x8408
unsigned short crc16i(char *data_p, unsigned short length)
{
    unsigned char i;
    unsigned int data;
    unsigned int crc;
    unsigned int bit;
    
    crc = 0xffff;
    if (length == 0)
        return (~crc);
    
    do {
        data = (unsigned int)0xff & *data_p++;
        bit = 128;
        for (i = 0; i < 8; i++, bit >>= 1)
        {
            if ((crc & 0x0001) ^ ((data & bit) / bit))
                crc = (crc >> 1) ^ POLY;
            else
                crc >>= 1;
        }
    } while (--length);
    
    // Reverse bits for compatibility
    
    data = 0;
    for (i=0; i<16; i++)
        if (crc & (1 << i))
            data |= (1 << (15 - i));
    
    return(data);
}

unsigned short crc16(char *data_p, unsigned short length)
{
    unsigned char i;
    unsigned int data;
    unsigned int crc;
    
    crc = 0xffff;
    if (length == 0)
        return (~crc);
    
    do {
        data = (unsigned int)0xff & *data_p++;
        for (i = 0; i < 8; i++, data >>= 1)
        {
            if ((crc & 0x0001) ^ (data & 0x0001))
                crc = (crc >> 1) ^ POLY;
            else
                crc >>= 1;
        }
    } while (--length);
    
    crc = ~crc;
    data = crc;
    crc = (crc << 8) | (data >> 8 & 0xFF);
    return (crc);
}

static const unsigned short crc16tab[256]= {
	0x0000,0x1021,0x2042,0x3063,0x4084,0x50a5,0x60c6,0x70e7,
	0x8108,0x9129,0xa14a,0xb16b,0xc18c,0xd1ad,0xe1ce,0xf1ef,
	0x1231,0x0210,0x3273,0x2252,0x52b5,0x4294,0x72f7,0x62d6,
	0x9339,0x8318,0xb37b,0xa35a,0xd3bd,0xc39c,0xf3ff,0xe3de,
	0x2462,0x3443,0x0420,0x1401,0x64e6,0x74c7,0x44a4,0x5485,
	0xa56a,0xb54b,0x8528,0x9509,0xe5ee,0xf5cf,0xc5ac,0xd58d,
	0x3653,0x2672,0x1611,0x0630,0x76d7,0x66f6,0x5695,0x46b4,
	0xb75b,0xa77a,0x9719,0x8738,0xf7df,0xe7fe,0xd79d,0xc7bc,
	0x48c4,0x58e5,0x6886,0x78a7,0x0840,0x1861,0x2802,0x3823,
	0xc9cc,0xd9ed,0xe98e,0xf9af,0x8948,0x9969,0xa90a,0xb92b,
	0x5af5,0x4ad4,0x7ab7,0x6a96,0x1a71,0x0a50,0x3a33,0x2a12,
	0xdbfd,0xcbdc,0xfbbf,0xeb9e,0x9b79,0x8b58,0xbb3b,0xab1a,
	0x6ca6,0x7c87,0x4ce4,0x5cc5,0x2c22,0x3c03,0x0c60,0x1c41,
	0xedae,0xfd8f,0xcdec,0xddcd,0xad2a,0xbd0b,0x8d68,0x9d49,
	0x7e97,0x6eb6,0x5ed5,0x4ef4,0x3e13,0x2e32,0x1e51,0x0e70,
	0xff9f,0xefbe,0xdfdd,0xcffc,0xbf1b,0xaf3a,0x9f59,0x8f78,
	0x9188,0x81a9,0xb1ca,0xa1eb,0xd10c,0xc12d,0xf14e,0xe16f,
	0x1080,0x00a1,0x30c2,0x20e3,0x5004,0x4025,0x7046,0x6067,
	0x83b9,0x9398,0xa3fb,0xb3da,0xc33d,0xd31c,0xe37f,0xf35e,
	0x02b1,0x1290,0x22f3,0x32d2,0x4235,0x5214,0x6277,0x7256,
	0xb5ea,0xa5cb,0x95a8,0x8589,0xf56e,0xe54f,0xd52c,0xc50d,
	0x34e2,0x24c3,0x14a0,0x0481,0x7466,0x6447,0x5424,0x4405,
	0xa7db,0xb7fa,0x8799,0x97b8,0xe75f,0xf77e,0xc71d,0xd73c,
	0x26d3,0x36f2,0x0691,0x16b0,0x6657,0x7676,0x4615,0x5634,
	0xd94c,0xc96d,0xf90e,0xe92f,0x99c8,0x89e9,0xb98a,0xa9ab,
	0x5844,0x4865,0x7806,0x6827,0x18c0,0x08e1,0x3882,0x28a3,
	0xcb7d,0xdb5c,0xeb3f,0xfb1e,0x8bf9,0x9bd8,0xabbb,0xbb9a,
	0x4a75,0x5a54,0x6a37,0x7a16,0x0af1,0x1ad0,0x2ab3,0x3a92,
	0xfd2e,0xed0f,0xdd6c,0xcd4d,0xbdaa,0xad8b,0x9de8,0x8dc9,
	0x7c26,0x6c07,0x5c64,0x4c45,0x3ca2,0x2c83,0x1ce0,0x0cc1,
	0xef1f,0xff3e,0xcf5d,0xdf7c,0xaf9b,0xbfba,0x8fd9,0x9ff8,
	0x6e17,0x7e36,0x4e55,0x5e74,0x2e93,0x3eb2,0x0ed1,0x1ef0
};
  
unsigned short crc16_ccitt(unsigned char *buf, int len)
{
	register int counter;
	register unsigned short crc = 0xFFFF;
	for( counter = 0; counter < len; counter++)
		crc = (crc<<8) ^ crc16tab[(crc>>8) ^ *buf++];
//		crc = (crc<<8) ^ crc16tab[((crc>>8) ^ *buf++)&0x00FF];
	return crc;
}
/* Main function */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])    
{ 
    int endTime;
    unsigned int clock5k;  // Keep track of 5kHz sample packets
    unsigned int endClock; // Number of sample packets to read.
	double *valuep;   // temporary variable for returning a matlab value.
    double *rates;    // Sample rates sent in RATES_ARG
    double *retvals;  // Return values array.
    int nrates;       // Number rates sent in RATES_ARG
    char *strp;
	char path[1024];
    char suffix[16];
	char sdcard[8];
    char sddrive[8];
	mxArray *intArray;
    mwSize dims[2];
    HANDLE hRawDisk;
    DWORD bytesRead;
    DWORD sector;
    DWORD r;
    DWORD sectorsPerCluster;
    DWORD bytesPerSector;
    DWORD numberOfFreeClusters;
    DWORD totalNumberOfClusters;
    LPVOID lpMsgBuf;
    unsigned char byte; // Current packet type byte from sector buffer.
    int chan;
    int iBuf;
    int crcstart;       // Index into buf[] where the next crc calculation should start.
    int nchanOrder;
    int errorCode;
    int sessionID;
    int hardwareID;
    unsigned short last_diginfo;  // Last diginfo header for 5k sample packet.
    unsigned int last_packet;     // last low bits of packet counter.
    unsigned int high_packet;     // high bits of packet counter.
    unsigned int packet_count;    // Number of packets since last 1 second marker.
    int fidEvents;          // File ID for event data
    int fidDigitals;        // File ID for digital data
    int fid[MAX_CHANNELS];  // File IDs for channel data.
    int fidAccel[MAX_ACCEL_CHANNELS]; // File IDs for accelerometer data.
    int chanOrder[4 * MAX_CHANNELS];
    int low_rate_value[MAX_LOWRATE_CHANNELS];  // Last value for each low rate channel.  8-bit and 16-bit channels will be unsigned values.
    char timestr[16] = "hh:mm:00"; // Access to the HH:MM string to put on the open button
    long hours;
    long minutes;
    long secs;
    unsigned long updatetime = 0;
    int rate_table[4] = {0, 5000, 10000, 20000};  // Conversion between rate code (0..3) and sampling rate.
    int niter=0;
    
    // Buffers for output files.
    unsigned int nevents;                 // Number of 32-bit words stored in event output buffers.
    unsigned int ndigitals;               // Number of elements stored in digital output buffers.
    unsigned int nsamples[MAX_CHANNELS];  // Number of samples stored in each sample output buffer.
    unsigned int naccel;                  // Number of accelerometers samples stored in the accel output buffer.
    
    short samples[MAX_CHANNELS][OUTBUF_SIZE];  // Analog data max rate is 20kHz
    unsigned short digitals[OUTBUF_SIZE / 4];  // Digitals are at 5kHz
    unsigned int events[OUTBUF_SIZE * 4];      // Event data rates can vary widely.  Make sure this is big enough.
    short accel[MAX_ACCEL_CHANNELS][OUTBUF_ACCEL_SIZE];  // Accelerometer channels are at 100Hz
    
    // Buffer for input sectors read from SDCard
    unsigned char buf[INBUF_SECTORS * 512];

    // Check for correct number of input and output arguments.

    if ((nlhs != 1) || (nrhs < 5))
    {
		mexErrMsgTxt("Usage: errcode = nc3save(Drive, Sector, Duration, Pathname, Samplerates)"); 
        return;
    }
        
    // Get input arguments
      
    sdcard[0] = 0;
    sddrive[0] = 0;
    strp = mxArrayToString(DRIVE_ARG);
    if (strp)
    {
        char ch = strp[0];
        if (((ch > 'C') && (ch < 'Z')) || ((ch > 'c') && (ch < 'z')))
        { // Good drive letter
            sprintf(sdcard, "\\\\.\\%c:", ch); // Raw disk name
            sprintf(sddrive, "%c:\\\\", ch);   // Drive prefix
        }
        mxFree(strp);
    }
    
	valuep = mxGetPr(SECTOR_ARG); 
	sector = (int)valuep[0];  
    
	valuep = mxGetPr(TIME_ARG); 
	endTime = valuep[0];
    if (endTime == 0)
        endTime = 345600;
    
    path[0] = 0;
	strp = mxArrayToString(PATH_ARG);
    if (strp)
    {
        strncpy(path, strp, 1024);
        mxFree(strp);
    }
 
    rates = mxGetPr(RATES_ARG);
    dims[0] = mxGetM(RATES_ARG);
    dims[1] = mxGetN(RATES_ARG);
    nrates = dims[0] * dims[1];
    
    // Setup return value
    
    errorCode = -1;
    dims[0] = 1;
    dims[1] = 1 + nrates;
    ERRCODE = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, 0);
    retvals = mxGetPr(ERRCODE);
    retvals[0] = errorCode; // negative return indicates parameter error. 

    // Check for bad arguments
	
    if (sdcard[0] == 0)
        mexErrMsgTxt("invalid drive letter");
    if (endTime < 0)
		mexErrMsgTxt("duration < 0"); 
    if (endTime > 345600) // 4 days.  extend this if needed.
		mexErrMsgTxt("duration > 345600 (4 days)"); 
    if (path[0] == 0)
		mexErrMsgTxt("invalid file path"); 
    if (nrates > MAX_CHANNELS)
		mexErrMsgTxt("Maximum channels exceeded"); 
    for (chan=0; chan<nrates; chan++)
    {  
        if ((rates[chan] != 0) && (rates[chan] != 5000) && (rates[chan] != 10000) && (rates[chan] != 20000))
        {
            mexErrMsgTxt("Unsported sample rate");
            return;
        }
    }
    
    if ((sdcard[0] == 0) || (endTime < 0) || (endTime > 345600) || (path[0] == 0) || (nrates > MAX_CHANNELS))
    {
        mexErrMsgTxt("Parameter Error: errcode = nc3save(Drive, Sector, Duration, Pathname, Samplerates)"); 
        return;
    }
    
    // Open the SDCard in raw disk mode
    
    hRawDisk = CreateFile(sdcard, // Must be exactly in this format
        GENERIC_READ,
        FILE_SHARE_READ | FILE_SHARE_WRITE, // Must have FILE_SHARE_WRITE
        NULL,
        OPEN_EXISTING, // Must be OPEN_EXISTING
        FILE_FLAG_NO_BUFFERING | FILE_FLAG_RANDOM_ACCESS, // Access will not be buffered even if this flag
        // is set so it is less confusing to specify it explicitly
        NULL);

    if (hRawDisk == INVALID_HANDLE_VALUE)
 		mexErrMsgTxt("Failed to open drive in raw access mode"); 

// This doesn't work if the SDCard is not initialized with a file system
//    if (!GetDiskFreeSpace(sddrive, &sectorsPerCluster, &bytesPerSector, &numberOfFreeClusters, &totalNumberOfClusters))
//    {
//       CloseHandle(hRawDisk);
// 	  mexErrMsgTxt("Failed to get drive info"); 
//    }
//     
//    mexPrintf("bytesPerSector %d", bytesPerSector);
//    if (bytesPerSector != 512)
//    {
//        CloseHandle(hRawDisk);
//        mexErrMsgTxt("Bytes per sector is not 512");
//    }

    //mexPrintf("Start sector = %d\n", sector);
    r = SetFilePointer(hRawDisk, sector * 512, NULL, FILE_BEGIN);
    if (r != sector * 512)
    {
        CloseHandle(hRawDisk);
        mexErrMsgTxt("Failed to set file pointer");
    }
    
    // Read first buffer of sectors from the SDCard
    
    if (!ReadFile(hRawDisk, (LPVOID)buf, ((DWORD)512) * ((DWORD)INBUF_SECTORS), &bytesRead, NULL))
    {
        CloseHandle(hRawDisk);
        if (FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM |
                FORMAT_MESSAGE_IGNORE_INSERTS, NULL, GetLastError(),
                MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPTSTR)&lpMsgBuf, 0, NULL))
        {
            mexPrintf("%s\n", (char*)lpMsgBuf);
            LocalFree(lpMsgBuf);
        }
        mexErrMsgTxt("Could not read from SDCard");
    }

    // Clear nsamples for first block.
    
    ndigitals = 0;
    nevents = 0;
    nchanOrder = 0;  // Number of elements in the chanOrder list.
    for (chan=0; chan<MAX_CHANNELS; chan++)
    {
        int i;
        nsamples[chan] = 0;
        fid[chan] = -1; // Clear output file list for analog channels
        for (i=0; i<OUTBUF_SIZE; i++) {
        	samples[chan][i] = 0; 
        }
    }
                
    for (chan=0; chan<MAX_LOWRATE_CHANNELS; chan++)
        low_rate_value[chan] == 0;
    
    naccel = 0;
    for (chan=0; chan<MAX_ACCEL_CHANNELS; chan++)
    {   int i;
        for (i=0; i<OUTBUF_ACCEL_SIZE; i++)
            accel[chan][i] == 0;
    }
    
    fidDigitals = makeFile(path, "_Digi00.u16");   // File ID for event data
    fidEvents = makeFile(path, "_Events.u32");  // File ID for digital data         
    fidAccel[0] = makeFile(path, "_AccelM.i16"); // Accelerometer Movement (squared magnitude)
    fidAccel[1] = makeFile(path, "_AccelX.i16"); // Accelerometer X Axis
    fidAccel[2] = makeFile(path, "_AccelY.i16"); // Accelerometer Y Axis
    fidAccel[3] = makeFile(path, "_AccelZ.i16"); // Accelerometer Z Axis
    fidAccel[4] = makeFile(path, "_AccelT.i16"); // Accelerometer Temperature

    // Process sector data until we reach endTime
    
    iBuf = 0;
    clock5k = 0;  // Number of 5kHz sample packets read so far.
    endClock = endTime;  // Convert to unsigned int.
    endClock *= 5000;    // Convert to total number of sample packets to read.
    sessionID = -1;      // Illegal session ID.
    hardwareID = -1;     // Illegal hardware ID.
    errorCode = 0;
    last_diginfo = 0x3FFF; // Set to catch starting Condition (0) and first data start event, while ignoring on-going stims and discrims.
    last_packet = 0;
    high_packet = 1;     // Illegal high_packet value.
    packet_count = 0;   // Clear packet counter.

// For testing, dump the first 2 SDCard sectors
//    for (iBuf=0; iBuf< 1023; iBuf++) {
//        mexPrintf("%d) 0x%02x\n", iBuf, buf[iBuf]);
//    }
//    iBuf = 0;
    
//Test code.  use rates from parameter while trouble shooting
//                 for (chan=0; chan<nrates; chan++)
//                 { // Open a file for each active channel.
//                     if (rates[chan] != 0)
//                     {
//                         chanOrder[nchanOrder++] = chan;
//                         sprintf(suffix, "_Chan%02d.i16", chan+1);
//                         fid[chan] = makeFile(path, suffix);
//                     }
//                   // Return the sample rate for each channel
//                     retvals[chan + 1] = rates[chan];
//                 }
// 
//                 for (chan=0; chan<nrates; chan++) 
//                     if (rates[chan] == 20000)
//                         chanOrder[nchanOrder++] = chan;
// 
//                 for (chan=0; chan<nrates; chan++) 
//                     if ((rates[chan] == 10000) || (rates[chan] == 20000))
//                         chanOrder[nchanOrder++] = chan;
// 
//                 for (chan=0; chan<nrates; chan++) 
//                     if (rates[chan] == 20000)
//                         chanOrder[nchanOrder++] = chan;
//                 
//                 // Make sure proper number of zero samples are in each output buffer.
//                 for (chan=0; chan<nrates; chan++)                    
//                     nsamples[chan] = clock5k * (rates[chan] / 5000);
    
    while (errorCode == 0)
    {
        // Check if we need to read the next set of sectors for the SDCard.
        // We do this when iBuf crosses into the last sector read so we
        // don't need to check every time iBuf increments.
        
        if (iBuf >= ((INBUF_SECTORS - 1) * 512))
        {
            // Copy last sector in buffer to start of buffer and
            // update indexes.
            int i;
            for (i=0; i<512; i++)
            {
                buf[i] = buf[i + ((INBUF_SECTORS - 1) * 512)];
            }
            sector += (INBUF_SECTORS - 1);
            iBuf -= ((INBUF_SECTORS - 1) * 512);
            //mexPrintf("Current sector = %d\n", sector);
            
            // Read next set of sectors
            
            if (!ReadFile(hRawDisk, (LPVOID)(buf + 512), ((INBUF_SECTORS - 1) * 512), &bytesRead, NULL))
            {
                CloseHandle(hRawDisk);
                if (FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM |
                FORMAT_MESSAGE_IGNORE_INSERTS, NULL, GetLastError(),
                MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPTSTR)&lpMsgBuf, 0, NULL))
                {
                    mexPrintf("%s\n", (char*)lpMsgBuf);
                    LocalFree(lpMsgBuf);
                }
                mexPrintf("End of SDCard reached.");
                errorCode = 1;
            }
        }
        
        // Read next byte from buffer.
        
//     mexPrintf("Bytes @%d: %d %d %d %d    %d %d %d %d\n", iBuf,
//             buf[iBuf], buf[iBuf+1], buf[iBuf+2], buf[iBuf+3],
//             buf[iBuf+4], buf[iBuf+5], buf[iBuf+6], buf[iBuf+7]
//     );
        crcstart = iBuf;
        byte = buf[iBuf++];
        //mexPrintf("iBuf %d = %d\n", iBuf, byte);
        if (byte == 1)
        {  // This requires byte alignment matching Intel pentium.
            unsigned int pvalue;
            unsigned int pid = *(unsigned short *)(buf + iBuf);
            iBuf += 2;
            pvalue = *(unsigned int *)(buf + iBuf);
            iBuf += 4;
            //mexPrintf("Param %d %d\n", pid, pvalue);

            if (pid == 0)
            { // Hardware version ID
                hardwareID = pvalue;
            }
            else if (pid == 1)
            { // Encoded first 16 biophys channel rates.
                for (chan=0; chan<16; chan++)
                    rates[chan] = rate_table[0x3 & (pvalue >> (chan << 1))];
            }
            else if (pid == 2)
            { // Encoded last 16 biophys channel rates.
                for (chan=0; chan<16; chan++)
                    rates[chan+16] = rate_table[0x3 & (pvalue >> (chan << 1))];
            }
            else if (pid == 3)
            { // Encoded AUX channel rates.
                for (chan=0; chan<3; chan++)
                    rates[chan+32] = rate_table[0x3 & (pvalue >> (chan << 1))];
            }
            else if (pid == 10000)
            { // Session ID Param is always the last param.
              // Look through channel rates list to match channel order produced
              // by the Neurochip3 CPU code.  Open data files for output.
              // Fix up number of zero samples that should currently be in each of the output buffers.
                sessionID = pvalue;
                mexPrintf("Param Session ID = %d\n", sessionID);
                nchanOrder = 0;  // Number of elements in the chanOrder list.
                
                for (chan=0; chan<nrates; chan++)
                { // Open a file for each active channel.
                    if (rates[chan] != 0)
                    {
                        chanOrder[nchanOrder++] = chan;
                        
                        if (chan < 32)
                            sprintf(suffix, "_Chan%02d.i16", chan+1); // Normal Intan channels are signed
                        else
                            sprintf(suffix, "_Chan%02d.u16", chan+1); // AUX channels are unsigned
                            
                        fid[chan] = makeFile(path, suffix);
                    }
                  // Return the sample rate for each channel
                    retvals[chan + 1] = rates[chan];
                }

                for (chan=0; chan<nrates; chan++) 
                    if (rates[chan] == 20000)
                        chanOrder[nchanOrder++] = chan;

                for (chan=0; chan<nrates; chan++) 
                    if ((rates[chan] == 10000) || (rates[chan] == 20000))
                        chanOrder[nchanOrder++] = chan;

                for (chan=0; chan<nrates; chan++) 
                    if (rates[chan] == 20000)
                        chanOrder[nchanOrder++] = chan;
                
                // Make sure proper number of zero samples are in each output buffer.
                for (chan=0; chan<nrates; chan++)                    
                    nsamples[chan] = clock5k * (rates[chan] / 5000);
            }
        }
        else if (byte == 0x02)
        {
            // Condition changed.
            int lrval = buf[iBuf++];
            events[nevents++] = 0;     // Condition Event ID
            events[nevents++] = lrval+1; // Condition #
            events[nevents++] = clock5k; 
            //mexPrintf("Cond %d %d ", lrval+1, clock5k);
           
        }
        else if (byte <= 0x0F)
        {
            errorCode = 2;
            mexPrintf("End Byte = %d\n", byte);
            // For testing, dump the surounding bytes in the current SDCard sector
//             iBuf--;
//             for (int i= -200; i<+200; i++) {
//                 if (iBuf+i >= 0)
//                     mexPrintf("%d) 0x%02x\n", i, buf[iBuf+i]);
//             }
        }
        else if (byte <= 0x3F)
        {
            // low rate data channels are currently stored as events
            // Accelerometer X = Chan16, Y = Chan17, Z = Chan18, ADC1 = Chan 19, ADC2 = Chan20, ADC3/Temp = Chan21
            int lrval = 0;
            int lrchan = byte - 0x10;  // LR Chan ID 0..47
            
            lrval |= buf[iBuf++];
            if (byte >= 0x20)
            {
                lrval |= (buf[iBuf++] << 8);
                if (byte < 0x30)
                {
                    if (lrval & 0x8000)
                        lrval |= 0xFFFF0000; // Sign extend 16-bit values.
                }
                else
                { // Handle 32-bit values
                    lrval |= (buf[iBuf++] << 16);
                    lrval |= (buf[iBuf++] << 24);
                }
            }
            low_rate_value[lrchan] = lrval;
            //mexPrintf("Data %d %d\n", lrchan, lrval);
        }
        else if (byte <= 0x4F)
        {
            // 16-bit emg data.  TODO save EMG data to a file
            int emgchan = byte - 0x40;  // EMG Chan ID 0..15
            int emgval = buf[iBuf++];
            emgval |= (buf[iBuf++] << 8);
            //mexPrintf("EMG %d %d\n", emgchan, emgval);
        }
        else if ((byte & 0x80) == 0)
        {
            errorCode = 3;
            mexPrintf("End Byte = %d\n", byte);
            // For testing, dump the surounding bytes in the current SDCard sector
//             iBuf--;
//             for (int i= -200; i<+200; i++) {
//                 if (iBuf+i >= 0)
//                     mexPrintf("%d) 0x%02x\n", i, buf[iBuf+i]);
//             }
        }
        else
        { // Handle 5kHz sample packet
            int diginfo = buf[iBuf++]; // Low byte is discriminator flags.
            diginfo |= byte << 8;      // High byte is stim/cond flags
            diginfo &= 0x7FFF;  // Clear packet flag to reuse as CRC error bit.
            
//            mexPrintf("iBuf %d\n", iBuf);
            if (byte & 0x40)
            { // Packet count and session ID are included in packet
                unsigned int session;
                unsigned int packet;
                
                packet = buf[iBuf++];
                packet |= (buf[iBuf++] << 8);
//              if (niter++ < 1000) mexPrintf("packet id= %d\n", packet);

                session = buf[iBuf++];
                session |= (buf[iBuf++] << 8);
                               
                if (sessionID == -1)
                    sessionID = session;  // First occurrence of sesssion ID could possibly happen before the session id Param.
                
                if (high_packet == 1)
                { // Initialize packet counter first time it is encountered.
                    high_packet = 0;
                    last_packet = packet;
                    mexPrintf("Found Session ID = %d\n", session);
                }
                
                if (packet < last_packet)
                    high_packet += 0x10000; // Carry bit.
                
                last_packet = packet;
                                
                if (sessionID != session)
                { // Session ID mismatch.  We are done.
                    errorCode = 4; 
                    mexPrintf("End of Session found\n");
                }
                else
                {
                    // Write output buffers once per second of data.
                    
                    for (chan=0; chan<MAX_CHANNELS; chan++)
                    {
                        if (fid[chan] != -1)
                        {
                            write(fid[chan], samples[chan], nsamples[chan] << 1);
                            nsamples[chan] = 0;
                        }
                    }
                    
                    for (chan=0; chan<MAX_ACCEL_CHANNELS; chan++)
                    {
                        if (fidAccel[chan] != -1)
                        {
                            write(fidAccel[chan], accel[chan], naccel << 1);
                        }
                    }
                    naccel = 0;
                    
                    write(fidEvents, events, nevents << 2);
                    nevents = 0;
                   
                    write(fidDigitals, digitals, ndigitals << 1);
                    ndigitals = 0;                
                }
                
                // Update text on button
                if (clock5k >= updatetime)
                {
                    //               //mexPrintf("Seconds: %g\n", (double)updatetime / 2000.0);
                    //               set(hOpenButton, 'String', ['Open ' num2str(updatetime / 2000)]);
                    //               prhs[2 5] = hOpenbutton;
                    //               prhs[3 6] = 'String';
                    //               prhs[4 7] = 'hh:mm';
                    hours = clock5k / 18000000;
                    minutes = (clock5k - (hours * 18000000)) / 300000;
                    timestr[0] = '0' + (hours - (hours % 10)) / 10;
                    timestr[1] = '0' + (hours % 10);
                    timestr[3] = '0' + (minutes - (minutes % 10)) / 10;
                    timestr[4] = '0' + (minutes % 10);
                    prhs[7] = mxCreateString(timestr);
                    mexCallMATLAB(0, NULL, 3, (mxArray **)&(prhs[5]), "set");
                    mexCallMATLAB(0, NULL, 0, NULL, "drawnow");
                    mxDestroyArray((mxArray *)prhs[7]);
                    updatetime = updatetime + 300000; // Updated every  60 seconds of data. 5k sample rate * 60 secs
                }
                
                // Stop if we will reach our endTime on this tick.
                
                if (clock5k + 1 >= endClock)
                    break;
                
                // Save 1 second marker event.
                
                events[nevents++] = 4;            // Event ID
                events[nevents++] = high_packet | packet;  // packet index.
                events[nevents++] = clock5k;
                packet_count = 0;       // Clear packet counter for next 1 second.
                
                //mexPrintf("clock %d\n", clock5k);
            }

            if (byte & 0x20)
            { // Sample data is included in packet
                int i;
                short int b1;
                unsigned short crcchk;
                unsigned short crc;
       
                for (i=0; i<nchanOrder; i++)
                {
                	b1 = buf[iBuf++];
                	b1 |= buf[iBuf++] << 8;
                    chan = chanOrder[i];
                    samples[chan][nsamples[chan]++] = b1;
//                    if ((niter++ <= 1000)) mexPrintf("%d) Ch%d = 0x%04x, 0x%04x\n", iBuf-3, chan, /*diginfo*/niter, b1);
                }
                
                // Handle CCITT 16 bit cyclic redundancy check (crc)
                //                                 16  12  5
                // The CCITT CRC 16 polynomial is X + X + X + 1.
                // In binary, this is the bit pattern 1 0001 0000 0010 0001,
                // and in hex it is is 0x11021.  A 17 bit register is simulated
                // by testing the MSB before shifting the data. This allows
                // the polynomial to be represented by the 16 bit value 0x1021
                // This calculation would normally loop over bytes and bits but
                // those loops have been unrolled for efficiency.  Even so,
                // this calculation is still over 40% of the conversion time.
                // Start with the most significant bit of the least significant byte.
                
//                crc = crc16_ccitt(buf + crcstart, iBuf-crcstart);
// or                
//                crc = crc16i(buf + crcstart, iBuf-crcstart);
// or                
//                crc = 0xFFFF;
//                 for (i=crcstart; i<iBuf; i++)
//                 {
//                     b1 = buf[i];
//                     if ((crc & 0x8000) ^ ((b1 & 0x80)? 0x8000 : 0)) crc = (crc << 1) ^ 0x1021;  else crc <<= 1;
//                     if ((crc & 0x8000) ^ ((b1 & 0x40)? 0x8000 : 0)) crc = (crc << 1) ^ 0x1021;  else crc <<= 1;
//                     if ((crc & 0x8000) ^ ((b1 & 0x20)? 0x8000 : 0)) crc = (crc << 1) ^ 0x1021;  else crc <<= 1;
//                     if ((crc & 0x8000) ^ ((b1 & 0x10)? 0x8000 : 0)) crc = (crc << 1) ^ 0x1021;  else crc <<= 1;
//                     if ((crc & 0x8000) ^ ((b1 & 0x08)? 0x8000 : 0)) crc = (crc << 1) ^ 0x1021;  else crc <<= 1;
//                     if ((crc & 0x8000) ^ ((b1 & 0x04)? 0x8000 : 0)) crc = (crc << 1) ^ 0x1021;  else crc <<= 1;
//                     if ((crc & 0x8000) ^ ((b1 & 0x02)? 0x8000 : 0)) crc = (crc << 1) ^ 0x1021;  else crc <<= 1;
//                     if ((crc & 0x8000) ^ ((b1 & 0x01)? 0x8000 : 0)) crc = (crc << 1) ^ 0x1021;  else crc <<= 1;
//                 }
                
                
                // Read CRC check for this sample packet.
                crcchk = buf[iBuf++];
                crcchk |= (buf[iBuf++] << 8);              

                crc = crcchk;// Ignoring CRC checks for now.
                if (crc != crcchk)
                {
                    diginfo |= 0x8000; // badPacket = 1;
                    events[nevents++] = 5;   // CRC Error Event ID
                    events[nevents++] = (crc << 16) | crcchk;   
                    events[nevents++] = clock5k;
                    if (crcstart > 250)
                        errorCode = 5;
                    //mexPrintf("crc 0x%04x != chk 0x%04x 0x%04x 0x%04x\n", crc, crcchk, crcstart, iBuf);                   
                }
            }
            else
            { // include zero data for this packet.
                int i;
                for (i=0; i<nchanOrder; i++)
                {
                    chan = chanOrder[i];
                    samples[chan][nsamples[chan]++] = 0;
                }
                //mexPrintf("skip %d\n", clock5k);
            }
            
            // Update Accel output buffers with last accelerometer.
            
            if (clock5k % 50 == 0)
            { // Every 50 clock5K ticks = 100 Hz.
                int M; // Vector magnitude.
                int X = low_rate_value[16];
                int Y = low_rate_value[17];
                int Z = low_rate_value[18];
                accel[1][naccel] = (X * 1000) >> 14;
                accel[2][naccel] = (Y * 1000) >> 14;
                accel[3][naccel] = (Z * 1000) >> 14;
                accel[4][naccel] = low_rate_value[19] >> 6; // Get Temperature low data rate value, convert to Celcius
                // Reconstruct the Magnitude
                X >>= 4;  // Reduce X,Y,Z to 12 bits each so sum of squares stays in an int32.
                Y >>= 4;
                Z >>= 4;
                M = sqrt(X*X + Y*Y + Z*Z);
                accel[0][naccel] = (M * 1000) >> 10;
                //mexPrintf("Accel %g, %g, %g, %g, %g, %g, %g\n", floor((double)M),(double)X,(double)Y,(double)Z,(double)low_rate_value[19],(double)low_rate_value[20],(double)low_rate_value[21]);
                naccel++;
            }
            
            // Update Digital info and events;
            
            digitals[ndigitals++] = diginfo;
 
            if (diginfo != last_diginfo)
            { // Something changed.  Check events we haven't done yet.
                int evt;  // event #
                int bit;  // bit mask for event
                
                for (evt = 0; evt <= 7; evt++)
                { // Check digital events for rising edges
                    bit = 1 << evt;
                    if (((diginfo & bit) != 0) && ((last_diginfo & bit) == 0))
                    { // Event occurred.
                        events[nevents++] = 1;    // Discrim Event ID
                        events[nevents++] = evt + 1;  // Discrim# 1..8
                        events[nevents++] = clock5k;
                    }
                }
                
                // Check Stim events on rising edges of bit 0x1000
                
                if (((diginfo & 0x1000) != 0) && ((last_diginfo & 0x1000) == 0))
                { // Event occurred
                    events[nevents++] = 2;        // Stim Event ID
                    events[nevents++] = ((diginfo >> 8) & 0x0F) + 1; // Stim # 1..16
                    events[nevents++] = clock5k;
                }
                                      
                // Sample start and stop events.
                if ((diginfo & 0x2000) != (last_diginfo & 0x2000))
                {
                    events[nevents++] = 3;        // Sample Start/Stop Event ID
                    events[nevents++] = (diginfo & 0x2000) >> 13; // 1 = start, 0 = stop
                    events[nevents++] = clock5k;
                    //mexPrintf("start %d %d\n", clock5k, (diginfo & 0x2000) >> 13);
                }
            }
            last_diginfo = diginfo;
            
            // Update clock tick.
            if (++packet_count > 5000)
            {
                errorCode = 5;
                mexPrintf("Packet marker not found.  Assuming end of file.\n");
            }
            clock5k++;
        }
    }
    mexPrintf("Total Seconds = %g\n", (double)clock5k / 5000.0);
    
    // Close the SDCard and other open files.
    
    CloseHandle(hRawDisk);
    close(fidEvents);
    close(fidDigitals);
    
    for (chan=0; chan<nrates; chan++)  // Close all open channel files.  
        if (fid[chan] != -1)
            close(fid[chan]);

    for (chan=0; chan<MAX_ACCEL_CHANNELS; chan++)
        if (fidAccel[chan] != -1)
            close(fidAccel[chan]);

    // Update final data duration.
    
    hours = clock5k / 18000000;
    minutes = (clock5k - (hours * 18000000)) / 300000;
    secs = (clock5k - (hours * 18000000) - (minutes * 300000)) / 5000;
    timestr[0] = '0' + (hours - (hours % 10)) / 10;
    timestr[1] = '0' + (hours % 10);
    timestr[3] = '0' + (minutes - (minutes % 10)) / 10;
    timestr[4] = '0' + (minutes % 10);
    timestr[6] = '0' + (secs - (secs % 10)) / 10;
    timestr[7] = '0' + (secs % 10);
    prhs[7] = mxCreateString(timestr);
    mexCallMATLAB(0, NULL, 3, (mxArray **)&(prhs[5]), "set");
    mexCallMATLAB(0, NULL, 0, NULL, "drawnow");
    mxDestroyArray((mxArray *)prhs[7]);

    // Return the Session ID if it was found, otherwise return an error code.
    
    if (sessionID > -1)
        retvals[0] = sessionID;
    else
        retvals[0] = -errorCode;
    return;
}


