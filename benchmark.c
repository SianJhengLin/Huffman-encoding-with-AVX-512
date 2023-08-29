#include <stdio.h>      /* fprintf, fopen, ftello64 */
#include <stdlib.h>     /* malloc, free */
#include <string.h>     /* strcat */
#include <sys/types.h>  /* stat64 */
#include <sys/stat.h>   /* stat64 */
#include <assert.h>     /* assert */
#include <immintrin.h>
#include <time.h>
#include "bitstream.h"
#include "huf.h"
#include "fse.h"
#include "mem.h"
#include "xxhash.h"


#define NBLOOPS    4
#define TIMELOOP   (CLOCKS_PER_SEC * 2)

#define KB *(1U<<10)
#define MB *(1U<<20)
#define GB *(1U<<30)
#define KNUTH               2654435761U
#define MAX_MEM             ((sizeof(void*)==4) ? (2 GB - 64 MB) : (9ULL GB))
#define DEFAULT_CHUNKSIZE   (32 KB)

#define DISPLAY(...) fprintf(stderr, __VA_ARGS__)

#if !defined(S_ISREG)
#  define S_ISREG(x) (((x) & S_IFMT) == S_IFREG)
#endif

#if defined(_MSC_VER)   /* Visual */
#  pragma warning(disable : 4127)        /*disable: C4127: conditional expression is constant */
#  include <intrin.h>
#endif

typedef struct
{
	unsigned id;
	char*  origBuffer;
	size_t origSize;
	char*  compressedBuffer;
	size_t compressedSize;
	char*  destBuffer;
	size_t destSize;
} chunkParameters_t;
static clock_t BMK_clockSpan(clock_t start)
{
	return clock() - start;   /* works even if overflow; max duration ~ 30mn */
}

static size_t BMK_findMaxMem(U64 requiredMem)
{
	size_t step = (64 MB);
	BYTE* testmem = NULL;

	requiredMem = (((requiredMem >> 26) + 1) << 26);
	requiredMem += 2 * step;
	if (requiredMem > MAX_MEM) requiredMem = MAX_MEM;

	while (!testmem) {
		requiredMem -= step;
		if (requiredMem <= step) {
			requiredMem = step + 64;
			break;
		}
		testmem = (BYTE*)malloc((size_t)requiredMem);
	}

	free(testmem);
	return (size_t)(requiredMem - step);
}


static U64 BMK_GetFileSize(const char* infilename)
{
	int r;
#if defined(_MSC_VER)
	struct _stat64 statbuf;
	r = _stat64(infilename, &statbuf);
#else
	struct stat statbuf;
	r = stat(infilename, &statbuf);
#endif
	if (r || !S_ISREG(statbuf.st_mode)) return 0;   /* No good... */
	return (U64)statbuf.st_size;
}

void BMK_benchMem(chunkParameters_t* chunkP, int nbChunks,
	const char* inFileName, int benchedSize,
	U64* totalCompressedSize, double* totalCompressionTime, double* totalDecompressionTime,
	int nbSymbols, int memLog)
{
	int trial, chunkNb;
	size_t cSize = 0;
	double fastestC = 100000000., fastestD = 100000000.;
	double ratio = 0.;
	U32 crcCheck = 0;
	int nbDecodeLoops = ((100 MB) / (benchedSize + 1)) + 1;
	U32 const crcOrig = XXH32(chunkP[0].origBuffer, benchedSize, 0);
	size_t const nameLength = strlen(inFileName);

	/* Init */
	if (nameLength > 17) inFileName += nameLength - 17;   /* display last 17 characters */
	if (nbSymbols == 3) {   /* switch to special mode */
		//BMK_benchMem285(chunkP, nbChunks, inFileName, benchedSize, totalCompressedSize, totalCompressionTime, totalDecompressionTime, memLog);
		return;
	}

	DISPLAY("\r%79s\r", "");
	for (trial = 1; trial <= NBLOOPS; trial++) {
		int nbLoops = 0;
		clock_t clockStart, clockDuration;

		/* Compression */
		DISPLAY("%1i-%-15.15s : %9i ->\r", trial, inFileName, benchedSize);
		{ int i; for (i = 0; i < benchedSize; i++) chunkP[0].compressedBuffer[i] = (char)i; }    /* warmimg up memory */

		clockStart = clock();
		while (clock() == clockStart);
		clockStart = clock();

		while (BMK_clockSpan(clockStart) < TIMELOOP) {
			for (chunkNb = 0; chunkNb < nbChunks; chunkNb++) {
				size_t const cBSize = HUF_compress2(
					chunkP[chunkNb].compressedBuffer, FSE_compressBound(chunkP[chunkNb].origSize),
					chunkP[chunkNb].origBuffer, chunkP[chunkNb].origSize,nbSymbols,memLog);

				if (FSE_isError(cBSize)) {
					DISPLAY("!!! Error compressing block %i  !!!!  => %s   \n",
						chunkNb, FSE_getErrorName(cBSize));
					return;
				}
				chunkP[chunkNb].compressedSize = cBSize;
			}
			nbLoops++;
		}

		clockDuration = BMK_clockSpan(clockStart);
		clockDuration += !clockDuration;  /* to avoid division by zero */

		if ((double)clockDuration < fastestC * nbLoops * CLOCKS_PER_SEC)
			fastestC = (double)clockDuration / CLOCKS_PER_SEC / nbLoops;
		cSize = 0;
		for (chunkNb = 0; chunkNb < nbChunks; chunkNb++)
			cSize += chunkP[chunkNb].compressedSize ? chunkP[chunkNb].compressedSize : chunkP[chunkNb].origSize;
		ratio = (double)cSize / (double)benchedSize * 100.;

		DISPLAY("%1i-%-15.15s : %9i -> %9i (%5.2f%%),%7.1f MB/s\r",
			trial, inFileName, (int)benchedSize,
			(int)cSize, ratio,
			(double)benchedSize / (1 MB) / fastestC);


		//continue;
	   // if (loopNb == nbIterations) DISPLAY("\n"); continue;   /* skip decompression */
		/* Decompression */
		{ int i; for (i = 0; i < benchedSize; i++) chunkP[0].destBuffer[i] = 0; }     /* zeroing area, for CRC checking */
		clockStart = clock();
		while (clock() == clockStart);
		clockStart = clock();
		for (nbLoops = 0; nbLoops < nbDecodeLoops; nbLoops++) {
			for (chunkNb = 0; chunkNb < nbChunks; chunkNb++) {
				size_t regenSize;

				switch (chunkP[chunkNb].compressedSize)
				{
				case 0:   /* not compressed block; just memcpy() it */
					regenSize = chunkP[chunkNb].origSize;
					memcpy(chunkP[chunkNb].destBuffer, chunkP[chunkNb].origBuffer, regenSize);
					break;
				case 1:   /* single value byte; just memset() it */
					regenSize = chunkP[chunkNb].origSize;
					memset(chunkP[chunkNb].destBuffer, chunkP[chunkNb].origBuffer[0], chunkP[chunkNb].origSize);
					break;
				default:
					regenSize = HUF_decompress(chunkP[chunkNb].destBuffer, chunkP[chunkNb].origSize,
						chunkP[chunkNb].compressedBuffer, chunkP[chunkNb].compressedSize);
				}

				if (0) {  /* debugging => look for wrong bytes */
					const char* src = chunkP[chunkNb].origBuffer;
					const char* regen = chunkP[chunkNb].destBuffer;
					size_t origSize = chunkP[chunkNb].origSize;
					size_t n;
					for (n = 0; (n < origSize) && (src[n] == regen[n]); n++);
					if (n < origSize) {
						DISPLAY("\n!!! %15s : Invalid block %i !!! pos %u/%u\n",
							inFileName, chunkNb, (U32)n, (U32)origSize);
						break;
					}
				}
				if (regenSize != chunkP[chunkNb].origSize) {
					DISPLAY("!! Error decompressing block %i of cSize %u !! => (%s)  \r",
						chunkNb, (U32)chunkP[chunkNb].compressedSize, FSE_getErrorName(regenSize));
					return;
				}
			}
		}
		clockDuration = BMK_clockSpan(clockStart);

		if (clockDuration > 0) {
			if ((double)clockDuration < fastestD * nbDecodeLoops * CLOCKS_PER_SEC)
				fastestD = (double)clockDuration / CLOCKS_PER_SEC / nbDecodeLoops;
			assert(fastestD > 1. / 1000000000);   /* avoid overflow */
			nbDecodeLoops = (U32)(1. / fastestD) + 1;   /* aims for ~1sec */
		}
		else {
			assert(nbDecodeLoops < 20000000);  /* avoid overflow */
			nbDecodeLoops *= 100;
		}
		DISPLAY("%1i-%-15.15s : %9i -> %9i (%5.2f%%),%7.1f MB/s ,%7.1f MB/s\r",
			trial, inFileName, (int)benchedSize,
			(int)cSize, ratio,
			(double)benchedSize / (1 MB) / fastestC,
			(double)benchedSize / (1 MB) / fastestD);

		/* CRC Checking */
		crcCheck = XXH32(chunkP[0].destBuffer, benchedSize, 0);
		if (crcOrig != crcCheck) {
			const char* src = chunkP[0].origBuffer;
			const char* fin = chunkP[0].destBuffer;
			const char* const srcStart = src;
			while (*src == *fin)
				src++, fin++;
			DISPLAY("\n!!! %15s : Invalid Checksum !!! pos %i/%i\n",
				inFileName, (int)(src - srcStart), benchedSize);
			break;
		}
	}

	if (crcOrig == crcCheck) {
		if (ratio < 100.)
			DISPLAY("%-17.17s : %9i -> %9i (%5.2f%%),%7.1f MB/s ,%7.1f MB/s\n",
				inFileName, (int)benchedSize,
				(int)cSize, ratio,
				(double)benchedSize / (1 MB) / fastestC,
				(double)benchedSize / (1 MB) / fastestD);
		else
			DISPLAY("%-17.17s : %9i -> %9i (%5.1f%%),%7.1f MB/s ,%7.1f MB/s \n",
				inFileName, (int)benchedSize,
				(int)cSize, ratio,
				(double)benchedSize / (1 MB) / fastestC,
				(double)benchedSize / (1 MB) / fastestD);
	}
	else DISPLAY("\n");
	*totalCompressedSize += cSize;
	*totalCompressionTime += fastestC;
	*totalDecompressionTime += fastestD;
}

/*void BMK_benchMem(chunkParameters_t* chunkP, int nbChunks,
	const char* inFileName, int benchedSize,
	U64* totalCompressedSize, double* totalCompressionTime, double* totalDecompressionTime,
	int nbSymbols, int memLog, unsigned singleStream)
{
	int trial, chunkNb;
	size_t cSize = 0;
	double fastestC = 100000000., fastestD = 100000000.;
	double ratio = 0.;
	//U32 crcCheck = 0;
	//int nbDecodeLoops = ((100 MB) / (benchedSize + 1)) + 1;
	//U32 const crcOrig = XXH32(chunkP[0].origBuffer, benchedSize, 0);
	size_t const nameLength = strlen(inFileName);

	if (nameLength > 17) inFileName += nameLength - 17;   
	//if (nbSymbols == 3) { 
	//	BMK_benchMem285(chunkP, nbChunks, inFileName, benchedSize, totalCompressedSize, totalCompressionTime, totalDecompressionTime, memLog);
	//	return;
	//}
	DISPLAY("\r%79s\r", "");
	for (trial = 1; trial <= NBLOOPS; trial++) {
		int nbLoops = 0;
		clock_t clockStart, clockDuration;

		DISPLAY("%1i-%-15.15s : %9i ->\r", trial, inFileName, benchedSize);
		{ int i; for (i = 0; i < benchedSize; i++) chunkP[0].compressedBuffer[i] = (char)i; }   
		clockStart = clock();
		while (clock() == clockStart);
		clockStart = clock();
		while (BMK_clockSpan(clockStart) < TIMELOOP) {
			for (chunkNb = 0; chunkNb < nbChunks; chunkNb++) {
				size_t const cBSize = HUF_compress2(
					chunkP[chunkNb].compressedBuffer, FSE_compressBound(chunkP[chunkNb].origSize),
					chunkP[chunkNb].origBuffer, chunkP[chunkNb].origSize,nbSymbols,memLog);
				if (FSE_isError(cBSize)) {
					DISPLAY("!!! Error compressing block %i  !!!!  => %s   \n",
						chunkNb, FSE_getErrorName(cBSize));
					return;
				}
				chunkP[chunkNb].compressedSize = cBSize;
			}
			nbLoops++;
		}
		clockDuration = BMK_clockSpan(clockStart);
		clockDuration += !clockDuration;  
		if ((double)clockDuration < fastestC * nbLoops * CLOCKS_PER_SEC)
			fastestC = (double)clockDuration / CLOCKS_PER_SEC / nbLoops;
		cSize = 0;
		for (chunkNb = 0; chunkNb < nbChunks; chunkNb++)
			cSize += chunkP[chunkNb].compressedSize ? chunkP[chunkNb].compressedSize : chunkP[chunkNb].origSize;
		ratio = (double)cSize / (double)benchedSize * 100.;
		DISPLAY("%1i-%-15.15s : %9i -> %9i (%5.2f%%),%7.1f MB/s\r",
			trial, inFileName, (int)benchedSize,
			(int)cSize, ratio,
			(double)benchedSize / (1 MB) / fastestC);
		//continue;
	}
	DISPLAY("\n");
 	*totalCompressedSize += cSize;
 	*totalCompressionTime += fastestC;
 	*totalDecompressionTime += fastestD;
 }
*/
int main(int argc, char** argv)
{
	FILE *input_fp = NULL;
	argc = argc;
	char* filenamelist = NULL;
	char fileNamesTable[100][100] = { 0 };
	int nbFiles = 0;
	filenamelist = argv[1];
	if ((input_fp = fopen(filenamelist, "r+")) == NULL)
	{
		printf("Can't open file\n");
	}
	while (!feof(input_fp))
	{
		char *temp = fgets(fileNamesTable[nbFiles++], 100, input_fp);
		fileNamesTable[nbFiles - 1][strlen(fileNamesTable[nbFiles - 1]) - 1] = '\0';
		temp = temp;
	}


	size_t chunkSize = DEFAULT_CHUNKSIZE;
	int fileIdx = 0;
	U64 totalSourceSize_simd = 0;
	U64 totalCompressedSize_simd = 0;
	double totalc_simd = 0.;
	double totald_simd = 0.;

	while (fileIdx < nbFiles-1) {
		const char* const inFileName = fileNamesTable[fileIdx++];
		FILE* const inFile = fopen(inFileName, "rb");
		U64    inFileSize;
		size_t benchedSize;
		char* orig_buff;
		int nbChunks;
		int maxCompressedChunkSize;
		size_t readSize;
		char* compressedBuffer; int compressedBuffSize;
		char* destBuffer;
		chunkParameters_t* chunkP;

		/* Check file existence */

		if (inFile == NULL) { DISPLAY("Pb opening %s\n", inFileName); return 11; }

		/* Memory size evaluation */
		inFileSize = BMK_GetFileSize(inFileName);
		if (inFileSize == 0) { DISPLAY("file is empty\n"); fclose(inFile); return 11; }
		benchedSize = (size_t)BMK_findMaxMem(inFileSize * 3) / 3;
		if ((U64)benchedSize > inFileSize) benchedSize = (size_t)inFileSize;
		if (benchedSize < inFileSize)
			DISPLAY("Not enough memory for '%s' full size; testing %i MB only...\n",
				inFileName, (int)(benchedSize >> 20));

		/* Allocation */
		chunkP = (chunkParameters_t*)malloc(((benchedSize / chunkSize) + 1) * sizeof(chunkParameters_t));
		orig_buff = (char*)malloc((size_t)benchedSize);
		nbChunks = (int)(benchedSize / chunkSize) + 1;

		maxCompressedChunkSize = (int)FSE_compressBound(chunkSize);
		compressedBuffSize = nbChunks * maxCompressedChunkSize;
		compressedBuffer = (char*)malloc((size_t)compressedBuffSize);
		destBuffer = (char*)malloc((size_t)benchedSize);

		if (!orig_buff || !compressedBuffer || !destBuffer || !chunkP) {
			DISPLAY("\nError: not enough memory!\n");
			free(orig_buff);
			free(compressedBuffer);
			free(destBuffer);
			free(chunkP);
			fclose(inFile);
			return 12;
		}

		/* Init chunks */
		{   int i;
		size_t remaining = benchedSize;
		char* in = orig_buff;
		char* out = compressedBuffer;
		char* dst = destBuffer;
		for (i = 0; i < nbChunks; i++) {
			chunkP[i].id = i;
			chunkP[i].origBuffer = in; in += chunkSize;
			if (remaining > chunkSize) {
				chunkP[i].origSize = chunkSize;
				remaining -= chunkSize;
			}
			else {
				chunkP[i].origSize = (int)remaining;
				remaining = 0;
			}
			chunkP[i].compressedBuffer = out; out += maxCompressedChunkSize;
			chunkP[i].compressedSize = 0;
			chunkP[i].destBuffer = dst; dst += chunkSize;
		}   }

		/* Fill input buffer */
		DISPLAY("Loading %s...       \r", inFileName);
		readSize = fread(orig_buff, 1, benchedSize, inFile);
		fclose(inFile);

		if (readSize != benchedSize) {
			DISPLAY("\nError: problem reading file '%s' (%i read, should be %i) !!    \n",
				inFileName, (int)readSize, (int)benchedSize);
			free(orig_buff);
			free(compressedBuffer);
			free(destBuffer);
			free(chunkP);
			return 13;
		}

		/* Bench */

		BMK_benchMem(chunkP, nbChunks,
			inFileName, (int)benchedSize,
			&totalCompressedSize_simd, &totalc_simd, &totald_simd,
			255, 12);
		totalSourceSize_simd += benchedSize;

		free(orig_buff);
		free(compressedBuffer);
		free(destBuffer);
		free(chunkP);
	}

	if (nbFiles > 1)
		DISPLAY("%-17.17s :%10llu ->%10llu (%5.2f%%), %6.1f MB/s , %6.1f MB/s\n", "  TOTAL",
		(long long unsigned int)totalSourceSize_simd, (long long unsigned int)totalCompressedSize_simd,
			(double)totalCompressedSize_simd / (double)totalSourceSize_simd*100.,
			(double)totalSourceSize_simd / totalc_simd / CLOCKS_PER_SEC,
			(double)totalSourceSize_simd / totald_simd / CLOCKS_PER_SEC);
	return 0;
}
