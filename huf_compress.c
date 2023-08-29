/* ******************************************************************
   Huffman encoder, part of New Generation Entropy library
   Copyright (C) 2013-2016, Yann Collet.

   BSD 2-Clause License (http://www.opensource.org/licenses/bsd-license.php)

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:

	   * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
	   * Redistributions in binary form must reproduce the above
   copyright notice, this list of conditions and the following disclaimer
   in the documentation and/or other materials provided with the
   distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

	You can contact the author at :
	- FSE+HUF source repository : https://github.com/Cyan4973/FiniteStateEntropy
	- Public forum : https://groups.google.com/forum/#!forum/lz4c
****************************************************************** */

/* **************************************************************
*  Compiler specifics
****************************************************************/
#ifdef _MSC_VER    /* Visual Studio */
#  pragma warning(disable : 4127)        /* disable: C4127: conditional expression is constant */
#endif


/* **************************************************************
*  Includes
****************************************************************/
#include <string.h>     /* memcpy, memset */
#include <stdio.h>      /* printf (debug) */
#include <math.h>
#include "compiler.h"
#include "bitstream.h"
#include "hist.h"
#define FSE_STATIC_LINKING_ONLY   /* FSE_optimalTableLog_internal */
#include "fse.h"        /* header compression */
#define HUF_STATIC_LINKING_ONLY
#include "huf.h"
#include "error_private.h"


/* **************************************************************
*  Error Management
****************************************************************/
#define HUF_isError ERR_isError
#define HUF_STATIC_ASSERT(c) DEBUG_STATIC_ASSERT(c)   /* use only *after* variable declarations */
#define CHECK_V_F(e, f) size_t const e = f; if (ERR_isError(e)) return e
#define CHECK_F(f)   { CHECK_V_F(_var_err__, f); }


/* **************************************************************
*  Utils
****************************************************************/
unsigned HUF_optimalTableLog(unsigned maxTableLog, size_t srcSize, unsigned maxSymbolValue)
{
	return FSE_optimalTableLog_internal(maxTableLog, srcSize, maxSymbolValue, 1);
}


/* *******************************************************
*  HUF : Huffman block compression
*********************************************************/
/* HUF_compressWeights() :
 * Same as FSE_compress(), but dedicated to huff0's weights compression.
 * The use case needs much less stack memory.
 * Note : all elements within weightTable are supposed to be <= HUF_TABLELOG_MAX.
 */
#define MAX_FSE_TABLELOG_FOR_HUFF_HEADER 6
size_t HUF_compressWeights(void* dst, size_t dstSize, const void* weightTable, size_t wtSize)
{
	BYTE* const ostart = (BYTE*)dst;
	BYTE* op = ostart;
	BYTE* const oend = ostart + dstSize;

	U32 maxSymbolValue = HUF_TABLELOG_MAX;
	U32 tableLog = MAX_FSE_TABLELOG_FOR_HUFF_HEADER;

	FSE_CTable CTable[FSE_CTABLE_SIZE_U32(MAX_FSE_TABLELOG_FOR_HUFF_HEADER, HUF_TABLELOG_MAX)];
	BYTE scratchBuffer[1 << MAX_FSE_TABLELOG_FOR_HUFF_HEADER];

	U32 count[HUF_TABLELOG_MAX + 1];
	S16 norm[HUF_TABLELOG_MAX + 1];

	/* init conditions */
	if (wtSize <= 1) return 0;  /* Not compressible */

	/* Scan input and build symbol stats */
	{   unsigned const maxCount = HIST_count_simple(count, &maxSymbolValue, weightTable, wtSize);   /* never fails */
	if (maxCount == wtSize) return 1;   /* only a single symbol in src : rle */
	if (maxCount == 1) return 0;        /* each symbol present maximum once => not compressible */
	}

	tableLog = FSE_optimalTableLog(tableLog, wtSize, maxSymbolValue);
	CHECK_F(FSE_normalizeCount(norm, tableLog, count, wtSize, maxSymbolValue));

	/* Write table description header */
	{   CHECK_V_F(hSize, FSE_writeNCount(op, oend - op, norm, maxSymbolValue, tableLog));
	op += hSize;
	}

	/* Compress */
	CHECK_F(FSE_buildCTable_wksp(CTable, norm, maxSymbolValue, tableLog, scratchBuffer, sizeof(scratchBuffer)));
	{   CHECK_V_F(cSize, FSE_compress_usingCTable(op, oend - op, weightTable, wtSize, CTable));
	if (cSize == 0) return 0;   /* not enough space for compressed data */
	op += cSize;
	}

	return op - ostart;
}


struct HUF_CElt_s {
	U16  val;
	BYTE nbBits;
};   /* typedef'd to HUF_CElt within "huf.h" */

/*! HUF_writeCTable() :
	`CTable` : Huffman tree to save, using huf representation.
	@return : size of saved CTable */
size_t HUF_writeCTable(void* dst, size_t maxDstSize,
	const HUF_CElt* CTable, U32 maxSymbolValue, U32 huffLog)
{
	BYTE bitsToWeight[HUF_TABLELOG_MAX + 1];   /* precomputed conversion table */
	BYTE huffWeight[HUF_SYMBOLVALUE_MAX];
	BYTE* op = (BYTE*)dst;
	U32 n;

	/* check conditions */
	if (maxSymbolValue > HUF_SYMBOLVALUE_MAX) return ERROR(maxSymbolValue_tooLarge);

	/* convert to weight */
	bitsToWeight[0] = 0;
	for (n = 1; n < huffLog + 1; n++)
		bitsToWeight[n] = (BYTE)(huffLog + 1 - n);
	for (n = 0; n < maxSymbolValue; n++)
		huffWeight[n] = bitsToWeight[CTable[n].nbBits];

	/* attempt weights compression by FSE */
	{   CHECK_V_F(hSize, HUF_compressWeights(op + 1, maxDstSize - 1, huffWeight, maxSymbolValue));
	if ((hSize > 1) & (hSize < maxSymbolValue / 2)) {   /* FSE compressed */
		op[0] = (BYTE)hSize;
		return hSize + 1;
	}   }

	/* write raw values as 4-bits (max : 15) */
	if (maxSymbolValue > (256 - 128)) return ERROR(GENERIC);   /* should not happen : likely means source cannot be compressed */
	if (((maxSymbolValue + 1) / 2) + 1 > maxDstSize) return ERROR(dstSize_tooSmall);   /* not enough space within dst buffer */
	op[0] = (BYTE)(128 /*special case*/ + (maxSymbolValue - 1));
	huffWeight[maxSymbolValue] = 0;   /* to be sure it doesn't cause msan issue in final combination */
	for (n = 0; n < maxSymbolValue; n += 2)
		op[(n / 2) + 1] = (BYTE)((huffWeight[n] << 4) + huffWeight[n + 1]);
	return ((maxSymbolValue + 1) / 2) + 1;
}


size_t HUF_readCTable(HUF_CElt* CTable, U32* maxSymbolValuePtr, const void* src, size_t srcSize)
{
	BYTE huffWeight[HUF_SYMBOLVALUE_MAX + 1];   /* init not required, even though some static analyzer may complain */
	U32 rankVal[HUF_TABLELOG_ABSOLUTEMAX + 1];   /* large enough for values from 0 to 16 */
	U32 tableLog = 0;
	U32 nbSymbols = 0;

	/* get symbol weights */
	CHECK_V_F(readSize, HUF_readStats(huffWeight, HUF_SYMBOLVALUE_MAX + 1, rankVal, &nbSymbols, &tableLog, src, srcSize));

	/* check result */
	if (tableLog > HUF_TABLELOG_MAX) return ERROR(tableLog_tooLarge);
	if (nbSymbols > *maxSymbolValuePtr + 1) return ERROR(maxSymbolValue_tooSmall);

	/* Prepare base value per rank */
	{   U32 n, nextRankStart = 0;
	for (n = 1; n <= tableLog; n++) {
		U32 current = nextRankStart;
		nextRankStart += (rankVal[n] << (n - 1));
		rankVal[n] = current;
	}   }

	/* fill nbBits */
	{   U32 n; for (n = 0; n < nbSymbols; n++) {
		const U32 w = huffWeight[n];
		CTable[n].nbBits = (BYTE)(tableLog + 1 - w);
	}   }

	/* fill val */
	{   U16 nbPerRank[HUF_TABLELOG_MAX + 2] = { 0 };  /* support w=0=>n=tableLog+1 */
	U16 valPerRank[HUF_TABLELOG_MAX + 2] = { 0 };
	{ U32 n; for (n = 0; n < nbSymbols; n++) nbPerRank[CTable[n].nbBits]++; }
	/* determine stating value per rank */
	valPerRank[tableLog + 1] = 0;   /* for w==0 */
	{   U16 min = 0;
	U32 n; for (n = tableLog; n > 0; n--) {  /* start at n=tablelog <-> w=1 */
		valPerRank[n] = min;     /* get starting value within each rank */
		min += nbPerRank[n];
		min >>= 1;
	}   }
	/* assign value within rank, symbol order */
	{ U32 n; for (n = 0; n < nbSymbols; n++) CTable[n].val = valPerRank[CTable[n].nbBits]++; }
	}

	*maxSymbolValuePtr = nbSymbols - 1;
	return readSize;
}

U32 HUF_getNbBits(const void* symbolTable, U32 symbolValue)
{
	const HUF_CElt* table = (const HUF_CElt*)symbolTable;
	assert(symbolValue <= HUF_SYMBOLVALUE_MAX);
	return table[symbolValue].nbBits;
}


typedef struct nodeElt_s {
	U32 count;
	U16 parent;
	BYTE byte;
	BYTE nbBits;
} nodeElt;

static U32 HUF_setMaxHeight(nodeElt* huffNode, U32 lastNonNull, U32 maxNbBits)
{
	const U32 largestBits = huffNode[lastNonNull].nbBits;
	if (largestBits <= maxNbBits) return largestBits;   /* early exit : no elt > maxNbBits */

	/* there are several too large elements (at least >= 2) */
	{   int totalCost = 0;
	const U32 baseCost = 1 << (largestBits - maxNbBits);
	U32 n = lastNonNull;

	while (huffNode[n].nbBits > maxNbBits) {
		totalCost += baseCost - (1 << (largestBits - huffNode[n].nbBits));
		huffNode[n].nbBits = (BYTE)maxNbBits;
		n--;
	}  /* n stops at huffNode[n].nbBits <= maxNbBits */
	while (huffNode[n].nbBits == maxNbBits) n--;   /* n end at index of smallest symbol using < maxNbBits */

	/* renorm totalCost */
	totalCost >>= (largestBits - maxNbBits);  /* note : totalCost is necessarily a multiple of baseCost */

	/* repay normalized cost */
	{   U32 const noSymbol = 0xF0F0F0F0;
	U32 rankLast[HUF_TABLELOG_MAX + 2];
	int pos;

	/* Get pos of last (smallest) symbol per rank */
	memset(rankLast, 0xF0, sizeof(rankLast));
	{   U32 currentNbBits = maxNbBits;
	for (pos = n; pos >= 0; pos--) {
		if (huffNode[pos].nbBits >= currentNbBits) continue;
		currentNbBits = huffNode[pos].nbBits;   /* < maxNbBits */
		rankLast[maxNbBits - currentNbBits] = pos;
	}   }

	while (totalCost > 0) {
		U32 nBitsToDecrease = BIT_highbit32(totalCost) + 1;
		for (; nBitsToDecrease > 1; nBitsToDecrease--) {
			U32 highPos = rankLast[nBitsToDecrease];
			U32 lowPos = rankLast[nBitsToDecrease - 1];
			if (highPos == noSymbol) continue;
			if (lowPos == noSymbol) break;
			{   U32 const highTotal = huffNode[highPos].count;
			U32 const lowTotal = 2 * huffNode[lowPos].count;
			if (highTotal <= lowTotal) break;
			}
		}
		/* only triggered when no more rank 1 symbol left => find closest one (note : there is necessarily at least one !) */
		/* HUF_MAX_TABLELOG test just to please gcc 5+; but it should not be necessary */
		while ((nBitsToDecrease <= HUF_TABLELOG_MAX) && (rankLast[nBitsToDecrease] == noSymbol))
			nBitsToDecrease++;
		totalCost -= 1 << (nBitsToDecrease - 1);
		if (rankLast[nBitsToDecrease - 1] == noSymbol)
			rankLast[nBitsToDecrease - 1] = rankLast[nBitsToDecrease];   /* this rank is no longer empty */
		huffNode[rankLast[nBitsToDecrease]].nbBits++;
		if (rankLast[nBitsToDecrease] == 0)    /* special case, reached largest symbol */
			rankLast[nBitsToDecrease] = noSymbol;
		else {
			rankLast[nBitsToDecrease]--;
			if (huffNode[rankLast[nBitsToDecrease]].nbBits != maxNbBits - nBitsToDecrease)
				rankLast[nBitsToDecrease] = noSymbol;   /* this rank is now empty */
		}
	}   /* while (totalCost > 0) */

	while (totalCost < 0) {  /* Sometimes, cost correction overshoot */
		if (rankLast[1] == noSymbol) {  /* special case : no rank 1 symbol (using maxNbBits-1); let's create one from largest rank 0 (using maxNbBits) */
			while (huffNode[n].nbBits == maxNbBits) n--;
			huffNode[n + 1].nbBits--;
			rankLast[1] = n + 1;
			totalCost++;
			continue;
		}
		huffNode[rankLast[1] + 1].nbBits--;
		rankLast[1]++;
		totalCost++;
	}   }   }   /* there are several too large elements (at least >= 2) */

	return maxNbBits;
}


typedef struct {
	U32 base;
	U32 current;
} rankPos;

static void HUF_sort(nodeElt* huffNode, const U32* count, U32 maxSymbolValue)
{
	rankPos rank[32];
	U32 n;

	memset(rank, 0, sizeof(rank));
	for (n = 0; n <= maxSymbolValue; n++) {
		U32 r = BIT_highbit32(count[n] + 1);
		rank[r].base++;
	}
	for (n = 30; n > 0; n--) rank[n - 1].base += rank[n].base;
	for (n = 0; n < 32; n++) rank[n].current = rank[n].base;
	for (n = 0; n <= maxSymbolValue; n++) {
		U32 const c = count[n];
		U32 const r = BIT_highbit32(c + 1) + 1;
		U32 pos = rank[r].current++;
		while ((pos > rank[r].base) && (c > huffNode[pos - 1].count)) {
			huffNode[pos] = huffNode[pos - 1];
			pos--;
		}
		huffNode[pos].count = c;
		huffNode[pos].byte = (BYTE)n;
	}
}


/** HUF_buildCTable_wksp() :
 *  Same as HUF_buildCTable(), but using externally allocated scratch buffer.
 *  `workSpace` must be aligned on 4-bytes boundaries, and be at least as large as a table of HUF_CTABLE_WORKSPACE_SIZE_U32 unsigned.
 */
#define STARTNODE (HUF_SYMBOLVALUE_MAX+1)
typedef nodeElt huffNodeTable[HUF_CTABLE_WORKSPACE_SIZE_U32];
size_t HUF_buildCTable_wksp(HUF_CElt* tree, const U32* count, U32 maxSymbolValue, U32 maxNbBits, void* workSpace, size_t wkspSize)
{
	nodeElt* const huffNode0 = (nodeElt*)workSpace;
	nodeElt* const huffNode = huffNode0 + 1;
	U32 n, nonNullRank;
	int lowS, lowN;
	U16 nodeNb = STARTNODE;
	U32 nodeRoot;

	/* safety checks */
	if (((size_t)workSpace & 3) != 0) return ERROR(GENERIC);  /* must be aligned on 4-bytes boundaries */
	if (wkspSize < sizeof(huffNodeTable)) return ERROR(workSpace_tooSmall);
	if (maxNbBits == 0) maxNbBits = HUF_TABLELOG_DEFAULT;
	if (maxSymbolValue > HUF_SYMBOLVALUE_MAX) return ERROR(maxSymbolValue_tooLarge);
	memset(huffNode0, 0, sizeof(huffNodeTable));

	/* sort, decreasing order */
	HUF_sort(huffNode, count, maxSymbolValue);

	/* init for parents */
	nonNullRank = maxSymbolValue;
	while (huffNode[nonNullRank].count == 0) nonNullRank--;
	lowS = nonNullRank; nodeRoot = nodeNb + lowS - 1; lowN = nodeNb;
	huffNode[nodeNb].count = huffNode[lowS].count + huffNode[lowS - 1].count;
	huffNode[lowS].parent = huffNode[lowS - 1].parent = nodeNb;
	nodeNb++; lowS -= 2;
	for (n = nodeNb; n <= nodeRoot; n++) huffNode[n].count = (U32)(1U << 30);
	huffNode0[0].count = (U32)(1U << 31);  /* fake entry, strong barrier */

	/* create parents */
	while (nodeNb <= nodeRoot) {
		U32 n1 = (huffNode[lowS].count < huffNode[lowN].count) ? lowS-- : lowN++;
		U32 n2 = (huffNode[lowS].count < huffNode[lowN].count) ? lowS-- : lowN++;
		huffNode[nodeNb].count = huffNode[n1].count + huffNode[n2].count;
		huffNode[n1].parent = huffNode[n2].parent = nodeNb;
		nodeNb++;
	}

	/* distribute weights (unlimited tree height) */
	huffNode[nodeRoot].nbBits = 0;
	for (n = nodeRoot - 1; n >= STARTNODE; n--)
		huffNode[n].nbBits = huffNode[huffNode[n].parent].nbBits + 1;
	for (n = 0; n <= nonNullRank; n++)
		huffNode[n].nbBits = huffNode[huffNode[n].parent].nbBits + 1;

	/* enforce maxTableLog */
	maxNbBits = HUF_setMaxHeight(huffNode, nonNullRank, maxNbBits);

	/* fill result into tree (val, nbBits) */
	{   U16 nbPerRank[HUF_TABLELOG_MAX + 1] = { 0 };
	U16 valPerRank[HUF_TABLELOG_MAX + 1] = { 0 };
	if (maxNbBits > HUF_TABLELOG_MAX) return ERROR(GENERIC);   /* check fit into table */
	for (n = 0; n <= nonNullRank; n++)
		nbPerRank[huffNode[n].nbBits]++;
	/* determine stating value per rank */
	{   U16 min = 0;
	for (n = maxNbBits; n > 0; n--) {
		valPerRank[n] = min;      /* get starting value within each rank */
		min += nbPerRank[n];
		min >>= 1;
	}   }
	for (n = 0; n <= maxSymbolValue; n++)
		tree[huffNode[n].byte].nbBits = huffNode[n].nbBits;   /* push nbBits per symbol, symbol order */
	for (n = 0; n <= maxSymbolValue; n++)
		tree[n].val = valPerRank[tree[n].nbBits]++;   /* assign value within rank, symbol order */
	}

	return maxNbBits;
}

/** HUF_buildCTable() :
 * @return : maxNbBits
 *  Note : count is used before tree is written, so they can safely overlap
 */
size_t HUF_buildCTable(HUF_CElt* tree, const U32* count, U32 maxSymbolValue, U32 maxNbBits)
{
	huffNodeTable nodeTable;
	return HUF_buildCTable_wksp(tree, count, maxSymbolValue, maxNbBits, nodeTable, sizeof(nodeTable));
}

static size_t HUF_estimateCompressedSize(HUF_CElt* CTable, const unsigned* count, unsigned maxSymbolValue)
{
	size_t nbBits = 0;
	int s;
	for (s = 0; s <= (int)maxSymbolValue; ++s) {
		nbBits += CTable[s].nbBits * count[s];
	}
	return nbBits >> 3;
}

static int HUF_validateCTable(const HUF_CElt* CTable, const unsigned* count, unsigned maxSymbolValue) {
	int bad = 0;
	int s;
	for (s = 0; s <= (int)maxSymbolValue; ++s) {
		bad |= (count[s] != 0) & (CTable[s].nbBits == 0);
	}
	return !bad;
}

size_t HUF_compressBound(size_t size) { return HUF_COMPRESSBOUND(size); }

FORCE_INLINE_TEMPLATE void
HUF_encodeSymbol(BIT_CStream_t* bitCPtr, U32 symbol, const HUF_CElt* CTable)
{
	BIT_addBitsFast(bitCPtr, CTable[symbol].val, CTable[symbol].nbBits);
}

#define HUF_FLUSHBITS(s)  BIT_flushBits(s)

#define HUF_FLUSHBITS_1(stream) \
    if (sizeof((stream)->bitContainer)*8 < HUF_TABLELOG_MAX*2+7) HUF_FLUSHBITS(stream)

#define HUF_FLUSHBITS_2(stream) \
    if (sizeof((stream)->bitContainer)*8 < HUF_TABLELOG_MAX*4+7) HUF_FLUSHBITS(stream)

size_t
HUF_compress1X_usingCTable_internal_body(void* dst, size_t dstSize,
	const void* src, size_t srcSize,
	const HUF_CElt* CTable)
{
	const BYTE* ip = (const BYTE*)src;
	BYTE* const ostart = (BYTE*)dst;
	BYTE* const oend = ostart + dstSize;
	BYTE* op = ostart;
	size_t n;
	BIT_CStream_t bitC;

	/* init */
	if (dstSize < 8) return 0;   /* not enough space to compress */
	{ size_t const initErr = BIT_initCStream(&bitC, op, oend - op);
	if (HUF_isError(initErr)) return 0; }

	n = srcSize & ~3;  /* join to mod 4 */
	switch (srcSize & 3)
	{
	case 3: HUF_encodeSymbol(&bitC, ip[n + 2], CTable);
		HUF_FLUSHBITS_2(&bitC);
		/* fall-through */
	case 2: HUF_encodeSymbol(&bitC, ip[n + 1], CTable);
		HUF_FLUSHBITS_1(&bitC);
		/* fall-through */
	case 1: HUF_encodeSymbol(&bitC, ip[n + 0], CTable);
		HUF_FLUSHBITS(&bitC);
		/* fall-through */
	case 0: /* fall-through */
	default: break;
	}

	for (; n > 0; n -= 4) {  /* note : n&3==0 at this stage */
		HUF_encodeSymbol(&bitC, ip[n - 1], CTable);
		HUF_FLUSHBITS_1(&bitC);
		HUF_encodeSymbol(&bitC, ip[n - 2], CTable);
		HUF_FLUSHBITS_2(&bitC);
		HUF_encodeSymbol(&bitC, ip[n - 3], CTable);
		HUF_FLUSHBITS_1(&bitC);
		HUF_encodeSymbol(&bitC, ip[n - 4], CTable);
		HUF_FLUSHBITS(&bitC);
	}

	return BIT_closeCStream(&bitC);
}

#if DYNAMIC_BMI2

static TARGET_ATTRIBUTE("bmi2") size_t
HUF_compress1X_usingCTable_internal_bmi2(void* dst, size_t dstSize,
	const void* src, size_t srcSize,
	const HUF_CElt* CTable)
{
	return HUF_compress1X_usingCTable_internal_body(dst, dstSize, src, srcSize, CTable);
}

static size_t
HUF_compress1X_usingCTable_internal_default(void* dst, size_t dstSize,
	const void* src, size_t srcSize,
	const HUF_CElt* CTable)
{
	return HUF_compress1X_usingCTable_internal_body(dst, dstSize, src, srcSize, CTable);
}

static size_t
HUF_compress1X_usingCTable_internal(void* dst, size_t dstSize,
	const void* src, size_t srcSize,
	const HUF_CElt* CTable, const int bmi2)
{
	if (bmi2) {
		return HUF_compress1X_usingCTable_internal_bmi2(dst, dstSize, src, srcSize, CTable);
	}
	return HUF_compress1X_usingCTable_internal_default(dst, dstSize, src, srcSize, CTable);
}

#else

static size_t
HUF_compress1X_usingCTable_internal(void* dst, size_t dstSize,
	const void* src, size_t srcSize,
	const HUF_CElt* CTable, const int bmi2)
{
	(void)bmi2;
	return HUF_compress1X_usingCTable_internal_body(dst, dstSize, src, srcSize, CTable);
}

#endif

size_t HUF_compress1X_usingCTable(void* dst, size_t dstSize, const void* src, size_t srcSize, const HUF_CElt* CTable)
{
	return HUF_compress1X_usingCTable_internal(dst, dstSize, src, srcSize, CTable, /* bmi2 */ 0);
}


/** Huffman_SIMD
*
*simd优化函数，单次循环处理128bytes
*
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////

size_t HUF_compressSIMD_internal_body(void* dst, size_t dstSize, const void* src, size_t srcSize, __m512i table[], const HUF_CElt* CTable)
{
	__m512i  flag_pad0, flag_pad1, flag_shuffle;
	U16 uflag_pad0[] = { 1, 1, 0, 0, 3, 3, 2, 2, 5, 5, 4, 4, 7, 7, 6, 6, 9, 9, 8, 8, 11, 11, 10, 10, 13, 13, 12, 12, 15, 15, 14, 14 };
	U16 uflag_pad1[] = { 17, 17, 16, 16, 19, 19, 18, 18, 21, 21, 20, 20, 23, 23, 22, 22, 25, 25, 24, 24, 27, 27, 26, 26, 29, 29, 28, 28, 31, 31, 30, 30 };
	BYTE uflag_shuffle[] = { 1,0,2,0,5,0,6,0,9,0,10,0,13,0,14,0,17,0,18,0,21,0,22,0,25,0,26,0,29,0,30,0,33,0,34,0,37,0,38,0,41,0,42,0,45,0,46,0,49,0,50,0,53,0,54,0,57,0,58,0,61,0,62,0 };
	flag_pad0 = _mm512_loadu_si512(uflag_pad0);
	flag_pad1 = _mm512_loadu_si512(uflag_pad1);
	flag_shuffle = _mm512_loadu_si512(uflag_shuffle);
	__m512i zero = _mm512_setzero_epi32();
	__m512i one = _mm512_set1_epi32(1);
	__m512i num64 = _mm512_set1_epi64(64ll);
	__m512i idx_length_32[4], idx_length_64[4], result[4];
	__m512i  ls_64[2], temp[4], idx_length_16_odd[4], idx_length_16_even[4];
	__mmask64 buf_64_length[8] = { 0 }, buf_64[8] = { 0 };
	__mmask32 pre_length = 0, k1[4], k2[4];

	const BYTE* ip = (const BYTE*)src;
	BYTE* const ostart = (BYTE*)dst;
	BYTE* const oend = ostart + dstSize;
	BYTE* op = ostart;
	BIT_CStream_t bitC;
	int i=0;


	if (dstSize < 8) return 0;
	{ size_t const initErr = BIT_initCStream(&bitC, op, oend - op);
	if (HUF_isError(initErr)) return 0; }
	size_t n = srcSize & ~3;
	switch (srcSize & 3)
	{
	case 3: HUF_encodeSymbol(&bitC, ip[n + 2], CTable);
		HUF_FLUSHBITS_2(&bitC);
	case 2: HUF_encodeSymbol(&bitC, ip[n + 1], CTable);
		HUF_FLUSHBITS_1(&bitC);

	case 1: HUF_encodeSymbol(&bitC, ip[n + 0], CTable);
		HUF_FLUSHBITS(&bitC);

	case 0:
	default: break;
	}
	pre_length = bitC.bitPos;
	ip += n;
	__mmask64 * ptr = (__mmask64*)bitC.ptr;
	for (; n >= 128; n -= 128)
	{
		//512bits
		ip -= 64;
		result[0] = _mm512_loadu_si512(ip);
		result[1] = _mm512_permutexvar_epi16(flag_pad1, result[0]);
		result[1] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[1], flag_shuffle);
		result[0] = _mm512_permutexvar_epi16(flag_pad0, result[0]);
		result[0] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[0], flag_shuffle);
		ip -= 64;
		result[2] = _mm512_loadu_si512(ip);
		result[3] = _mm512_permutexvar_epi16(flag_pad1, result[2]);
		result[3] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[3], flag_shuffle);
		result[2] = _mm512_permutexvar_epi16(flag_pad0, result[2]);
		result[2] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[2], flag_shuffle);
		//设置mask,mask前5位标记子表中的位置，第6位标记每一组表中第几个子表，第7-8位标记第几组表。
		//查表
		temp[0] = _mm512_slli_epi16(result[0], 8);
		temp[1] = _mm512_slli_epi16(result[1], 8);
		k1[0] = _mm512_movepi16_mask(temp[0]);
		k1[1] = _mm512_movepi16_mask(temp[1]);
		k2[0] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[0], 1));
		k2[1] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[1], 1));
		temp[2] = _mm512_slli_epi16(result[2], 8);
		temp[3] = _mm512_slli_epi16(result[3], 8);
		k1[2] = _mm512_movepi16_mask(temp[2]);
		k1[3] = _mm512_movepi16_mask(temp[3]);
		k2[2] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[2], 1));
		k2[3] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[3], 1));

		result[0] = _mm512_mask2_permutex2var_epi16(table[0], result[0], ~k1[0] & ~k2[0], table[1]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[0], result[1], ~k1[1] & ~k2[1], table[1]);
		result[0] = _mm512_mask2_permutex2var_epi16(table[2], result[0], ~k1[0] & k2[0], table[3]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[2], result[1], ~k1[1] & k2[1], table[3]);
		result[0] = _mm512_mask2_permutex2var_epi16(table[4], result[0], k1[0] & ~k2[0], table[5]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[4], result[1], k1[1] & ~k2[1], table[5]);
		result[0] = _mm512_mask2_permutex2var_epi16(table[6], result[0], k1[0] & k2[0], table[7]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[6], result[1], k1[1] & k2[1], table[7]);

		result[2] = _mm512_mask2_permutex2var_epi16(table[0], result[2], ~k1[2] & ~k2[2], table[1]);
		result[3] = _mm512_mask2_permutex2var_epi16(table[0], result[3], ~k1[3] & ~k2[3], table[1]);
		result[2] = _mm512_mask2_permutex2var_epi16(table[2], result[2], ~k1[2] & k2[2], table[3]);
		result[3] = _mm512_mask2_permutex2var_epi16(table[2], result[3], ~k1[3] & k2[3], table[3]);
		result[2] = _mm512_mask2_permutex2var_epi16(table[4], result[2], k1[2] & ~k2[2], table[5]);
		result[3] = _mm512_mask2_permutex2var_epi16(table[4], result[3], k1[3] & ~k2[3], table[5]);
		result[2] = _mm512_mask2_permutex2var_epi16(table[6], result[2], k1[2] & k2[2], table[7]);
		result[3] = _mm512_mask2_permutex2var_epi16(table[6], result[3], k1[3] & k2[3], table[7]);
		//获取移位长度
		//移位	
		idx_length_16_odd[0] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[0]), one);
		idx_length_16_odd[1] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[1]), one);
		idx_length_16_even[0] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[0], 16)), one);
		idx_length_16_even[1] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[1], 16)), one);
		idx_length_32[0] = _mm512_add_epi32(idx_length_16_even[0], idx_length_16_odd[0]);
		idx_length_32[1] = _mm512_add_epi32(idx_length_16_even[1], idx_length_16_odd[1]);
		result[0] = _mm512_sllv_epi16(result[0], idx_length_16_even[0]);//左移
		result[1] = _mm512_sllv_epi16(result[1], idx_length_16_even[1]);//左移
		idx_length_16_odd[0] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[0]);
		idx_length_16_odd[1] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[1]);
		idx_length_16_even[0] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[0]);
		idx_length_16_even[1] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[1]);
		result[0] = _mm512_sllv_epi32(result[0], idx_length_16_odd[0]);//左移
		result[1] = _mm512_sllv_epi32(result[1], idx_length_16_odd[1]);//左移
		result[0] = _mm512_srlv_epi32(result[0], idx_length_16_even[0]);//右移
		result[1] = _mm512_srlv_epi32(result[1], idx_length_16_even[1]);//右移
		ls_64[0] = _mm512_srli_epi64(idx_length_32[0], 32);
		ls_64[1] = _mm512_srli_epi64(idx_length_32[1], 32);
		result[0] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[0], ls_64[0]), ls_64[0]);//左移
		result[1] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[1], ls_64[1]), ls_64[1]);//左移
		ls_64[0] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[0]);
		ls_64[1] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[1]);
		result[0] = _mm512_srlv_epi64(result[0], ls_64[0]);
		result[1] = _mm512_srlv_epi64(result[1], ls_64[1]);
		idx_length_64[0] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[0], _mm512_slli_epi64(idx_length_32[0], 32)), 32));
		idx_length_64[1] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[1], _mm512_slli_epi64(idx_length_32[1], 32)), 32));

		idx_length_16_odd[2] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[2]), one);
		idx_length_16_odd[3] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[3]), one);
		idx_length_16_even[2] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[2], 16)), one);
		idx_length_16_even[3] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[3], 16)), one);
		idx_length_32[2] = _mm512_add_epi32(idx_length_16_even[2], idx_length_16_odd[2]);
		idx_length_32[3] = _mm512_add_epi32(idx_length_16_even[3], idx_length_16_odd[3]);
		result[2] = _mm512_sllv_epi16(result[2], idx_length_16_even[2]);//左移
		result[3] = _mm512_sllv_epi16(result[3], idx_length_16_even[3]);//左移
		idx_length_16_odd[2] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[2]);
		idx_length_16_odd[3] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[3]);
		idx_length_16_even[2] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[2]);
		idx_length_16_even[3] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[3]);
		result[2] = _mm512_sllv_epi32(result[2], idx_length_16_odd[2]);//左移
		result[3] = _mm512_sllv_epi32(result[3], idx_length_16_odd[3]);//左移
		result[2] = _mm512_srlv_epi32(result[2], idx_length_16_even[2]);//右移
		result[3] = _mm512_srlv_epi32(result[3], idx_length_16_even[3]);//右移
		ls_64[2] = _mm512_srli_epi64(idx_length_32[2], 32);
		ls_64[3] = _mm512_srli_epi64(idx_length_32[3], 32);
		result[2] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[2], ls_64[2]), ls_64[2]);//左移
		result[3] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[3], ls_64[3]), ls_64[3]);//左移
		ls_64[2] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[2]);;
		ls_64[3] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[3]);;
		result[2] = _mm512_srlv_epi64(result[2], ls_64[2]);
		result[3] = _mm512_srlv_epi64(result[3], ls_64[3]);
		idx_length_64[2] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[2], _mm512_slli_epi64(idx_length_32[2], 32)), 32));
		idx_length_64[3] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[3], _mm512_slli_epi64(idx_length_32[3], 32)), 32));
		//写出数据	
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[1]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[1]);

		for (i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];
			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[0]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[0]);

		for (i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];

			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}

		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[3]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[3]);

		for (i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];
			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[2]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[2]);

		for (i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];

			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
	}
	if(n>=64)
	{
		//512bits
		ip -= 64;
		n-=64;
		result[0] = _mm512_loadu_si512(ip);
		result[1] = _mm512_permutexvar_epi16(flag_pad1, result[0]);
		result[1] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[1], flag_shuffle);
		result[0] = _mm512_permutexvar_epi16(flag_pad0, result[0]);
		result[0] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[0], flag_shuffle);
		//设置mask,mask前5位标记子表中的位置，第6位标记每一组表中第几个子表，第7-8位标记第几组表。
		//查表
		temp[0] = _mm512_slli_epi16(result[0], 8);
		temp[1] = _mm512_slli_epi16(result[1], 8);
		k1[0] = _mm512_movepi16_mask(temp[0]);
		k1[1] = _mm512_movepi16_mask(temp[1]);
		k2[0] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[0], 1));
		k2[1] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[1], 1));

		result[0] = _mm512_mask2_permutex2var_epi16(table[0], result[0], ~k1[0] & ~k2[0], table[1]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[0], result[1], ~k1[1] & ~k2[1], table[1]);
		result[0] = _mm512_mask2_permutex2var_epi16(table[2], result[0], ~k1[0] & k2[0], table[3]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[2], result[1], ~k1[1] & k2[1], table[3]);
		result[0] = _mm512_mask2_permutex2var_epi16(table[4], result[0], k1[0] & ~k2[0], table[5]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[4], result[1], k1[1] & ~k2[1], table[5]);
		result[0] = _mm512_mask2_permutex2var_epi16(table[6], result[0], k1[0] & k2[0], table[7]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[6], result[1], k1[1] & k2[1], table[7]);
		//获取移位长度
		//移位	
		idx_length_16_odd[0] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[0]), one);
		idx_length_16_odd[1] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[1]), one);
		idx_length_16_even[0] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[0], 16)), one);
		idx_length_16_even[1] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[1], 16)), one);
		idx_length_32[0] = _mm512_add_epi32(idx_length_16_even[0], idx_length_16_odd[0]);
		idx_length_32[1] = _mm512_add_epi32(idx_length_16_even[1], idx_length_16_odd[1]);
		result[0] = _mm512_sllv_epi16(result[0], idx_length_16_even[0]);//左移
		result[1] = _mm512_sllv_epi16(result[1], idx_length_16_even[1]);//左移
		idx_length_16_odd[0] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[0]);
		idx_length_16_odd[1] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[1]);
		idx_length_16_even[0] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[0]);
		idx_length_16_even[1] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[1]);
		result[0] = _mm512_sllv_epi32(result[0], idx_length_16_odd[0]);//左移
		result[1] = _mm512_sllv_epi32(result[1], idx_length_16_odd[1]);//左移
		result[0] = _mm512_srlv_epi32(result[0], idx_length_16_even[0]);//右移
		result[1] = _mm512_srlv_epi32(result[1], idx_length_16_even[1]);//右移
		ls_64[0] = _mm512_srli_epi64(idx_length_32[0], 32);
		ls_64[1] = _mm512_srli_epi64(idx_length_32[1], 32);
		result[0] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[0], ls_64[0]), ls_64[0]);//左移
		result[1] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[1], ls_64[1]), ls_64[1]);//左移

		ls_64[0] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[0]);;
		ls_64[1] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[1]);;
		result[0] = _mm512_srlv_epi64(result[0], ls_64[0]);
		result[1] = _mm512_srlv_epi64(result[1], ls_64[1]);
		idx_length_64[0] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[0], _mm512_slli_epi64(idx_length_32[0], 32)), 32));
		idx_length_64[1] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[1], _mm512_slli_epi64(idx_length_32[1], 32)), 32));
		//写出数据	
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[1]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[1]);

		for (i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];
			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[0]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[0]);

		for (i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];

			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
	}
	if (n > 32&&n<64)
	{//处理长度不足64bytes的数据
		BYTE temp_in[64] = { 0 };
		ip = src;
		memcpy(temp_in + 64 - n, ip, n);
		result[0] = _mm512_loadu_si512(temp_in);
		result[1] = _mm512_permutexvar_epi16(flag_pad1, result[0]);
		result[1] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[1], flag_shuffle);
		result[0] = _mm512_permutexvar_epi16(flag_pad0, result[0]);
		result[0] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[0], flag_shuffle);
		//设置mask,mask前5位标记子表中的位置，第6位标记每一组表中第几个子表，第7-8位标记第几组表。
		//查表
		temp[0] = _mm512_slli_epi16(result[0], 8);
		temp[1] = _mm512_slli_epi16(result[1], 8);
		k1[0] = _mm512_movepi16_mask(temp[0]);
		k1[1] = _mm512_movepi16_mask(temp[1]);
		k2[0] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[0], 1));
		k2[1] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[1], 1));

		result[0] = _mm512_mask2_permutex2var_epi16(table[0], result[0], ~k1[0] & ~k2[0], table[1]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[0], result[1], ~k1[1] & ~k2[1], table[1]);
		result[0] = _mm512_mask2_permutex2var_epi16(table[2], result[0], ~k1[0] & k2[0], table[3]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[2], result[1], ~k1[1] & k2[1], table[3]);
		result[0] = _mm512_mask2_permutex2var_epi16(table[4], result[0], k1[0] & ~k2[0], table[5]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[4], result[1], k1[1] & ~k2[1], table[5]);
		result[0] = _mm512_mask2_permutex2var_epi16(table[6], result[0], k1[0] & k2[0], table[7]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[6], result[1], k1[1] & k2[1], table[7]);

		__mmask32 mask_zero = (__mmask32)((1llu << (n - 32)) - 1) << (64 - n);
		result[0] = _mm512_mask_blend_epi16(mask_zero, zero, result[0]);
		//获取移位长度
		//移位	
		idx_length_16_odd[0] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[0]), one);
		idx_length_16_odd[1] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[1]), one);
		idx_length_16_even[0] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[0], 16)), one);
		idx_length_16_even[1] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[1], 16)), one);
		idx_length_32[0] = _mm512_add_epi32(idx_length_16_even[0], idx_length_16_odd[0]);
		idx_length_32[1] = _mm512_add_epi32(idx_length_16_even[1], idx_length_16_odd[1]);
		result[0] = _mm512_sllv_epi16(result[0], idx_length_16_even[0]);//左移
		result[1] = _mm512_sllv_epi16(result[1], idx_length_16_even[1]);//左移
		idx_length_16_odd[0] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[0]);
		idx_length_16_odd[1] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[1]);
		idx_length_16_even[0] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[0]);
		idx_length_16_even[1] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[1]);
		result[0] = _mm512_sllv_epi32(result[0], idx_length_16_odd[0]);//左移
		result[1] = _mm512_sllv_epi32(result[1], idx_length_16_odd[1]);//左移
		result[0] = _mm512_srlv_epi32(result[0], idx_length_16_even[0]);//右移
		result[1] = _mm512_srlv_epi32(result[1], idx_length_16_even[1]);//右移
		ls_64[0] = _mm512_srli_epi64(idx_length_32[0], 32);
		ls_64[1] = _mm512_srli_epi64(idx_length_32[1], 32);
		result[0] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[0], ls_64[0]), ls_64[0]);//左移
		result[1] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[1], ls_64[1]), ls_64[1]);//左移

		ls_64[0] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[0]);;
		ls_64[1] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[1]);;
		result[0] = _mm512_srlv_epi64(result[0], ls_64[0]);
		result[1] = _mm512_srlv_epi64(result[1], ls_64[1]);
		idx_length_64[0] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[0], _mm512_slli_epi64(idx_length_32[0], 32)), 32));
		idx_length_64[1] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[1], _mm512_slli_epi64(idx_length_32[1], 32)), 32));
		//写出数据	
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[1]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[1]);

		for (i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];
			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[0]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[0]);

		for (i = 7;(i>=0) && buf_64[i]; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];
			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
	}
	if (n <= 32 && n)
	{//处理长度不足32bytes的数据

		BYTE temp_in[64] = { 0 };
		ip = src;
		memcpy(temp_in + 64 - n, ip, n);
		result[0] = _mm512_loadu_si512(temp_in);
		result[1] = _mm512_permutexvar_epi16(flag_pad1, result[0]);
		result[1] = _mm512_maskz_shuffle_epi8(0x5555555555555555llu, result[1], flag_shuffle);
		temp[1] = _mm512_slli_epi16(result[1], 8);
		k1[1] = _mm512_movepi16_mask(temp[1]);
		k2[1] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[1], 1));
		result[1] = _mm512_mask2_permutex2var_epi16(table[0], result[1], ~k1[1] & ~k2[1], table[1]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[2], result[1], ~k1[1] & k2[1], table[3]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[4], result[1], k1[1] & ~k2[1], table[5]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[6], result[1], k1[1] & k2[1], table[7]);
		__mmask32 mask_zero = (__mmask32)((1llu << n) - 1) << (32 - n);
		result[1] = _mm512_mask_blend_epi16(mask_zero, zero, result[1]);

		idx_length_16_odd[1] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[1]), one);
		idx_length_16_even[1] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[1], 16)), one);
		idx_length_32[1] = _mm512_add_epi32(idx_length_16_even[1], idx_length_16_odd[1]);
		result[1] = _mm512_sllv_epi16(result[1], idx_length_16_even[1]);//左移
		idx_length_16_odd[1] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[1]);
		idx_length_16_even[1] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[1]);
		result[1] = _mm512_sllv_epi32(result[1], idx_length_16_odd[1]);//左移
		result[1] = _mm512_srlv_epi32(result[1], idx_length_16_even[1]);//右移
		ls_64[1] = _mm512_srli_epi64(idx_length_32[1], 32);
		result[1] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[1], ls_64[1]), ls_64[1]);//左移
		ls_64[1] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[1]);
		result[1] = _mm512_srlv_epi64(result[1], ls_64[1]);
		idx_length_64[1] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[1], _mm512_slli_epi64(idx_length_32[1], 32)), 32));
		//写出数据	
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[1]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[1]);

		for (i = 7; (i>=0) && buf_64[i]; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			//printf("%d,%llu\n",i,ptr[0]);
			pre_length += buf_64_length[i];
			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
	}
	*ptr |= 1llu << pre_length;
	pre_length++;
	size_t const nbBytes = pre_length >> 3;
	return ((BYTE*)ptr - ostart) + nbBytes + ((pre_length & 7) > 0);
}

/* 256byte */
//单次循环处理256byte，
/*
size_t HUF_compressSIMD_internal_body(void* dst, size_t dstSize, const void* src, size_t srcSize, __m512i table[], const HUF_CElt* CTable)
{
	__m512i  flag_pad0, flag_pad1, flag_shuffle;
	const __m512i zero = _mm512_setzero_epi32();
	const __m512i one = _mm512_set1_epi32(1ll);
	const __m512i num64 = _mm512_set1_epi64(64ll);
	__m512i idx_length_32[8], idx_length_64[8], result[8];
	__m512i  ls_64[8], temp[8], idx_length_16_odd[8], idx_length_16_even[8];
	__mmask64 buf_64_length[8] = { 0 }, buf_64[8] = { 0 };
	__mmask64 pre_length = 0;
	__mmask32 k1[8], k2[8];
	U16 uflag_pad1[] = { 1, 1, 0, 0, 3, 3, 2, 2, 5, 5, 4, 4, 7, 7, 6, 6, 9, 9, 8, 8, 11, 11, 10, 10, 13, 13, 12, 12, 15, 15, 14, 14 };
	U16 uflag_pad2[] = { 17, 17, 16, 16, 19, 19, 18, 18, 21, 21, 20, 20, 23, 23, 22, 22, 25, 25, 24, 24, 27, 27, 26, 26, 29, 29, 28, 28, 31, 31, 30, 30 };
	BYTE uflag_shuffle[] = { 1,0,2,0,5,0,6,0,9,0,10,0,13,0,14,0,17,0,18,0,21,0,22,0,25,0,26,0,29,0,30,0,33,0,34,0,37,0,38,0,41,0,42,0,45,0,46,0,49,0,50,0,53,0,54,0,57,0,58,0,61,0,62,0 };
	flag_pad0 = _mm512_loadu_si512(uflag_pad1);
	flag_pad1 = _mm512_loadu_si512(uflag_pad2);
	flag_shuffle = _mm512_loadu_si512(uflag_shuffle);

	const BYTE* ip = (const BYTE*)src;
	BYTE* const ostart = (BYTE*)dst;
	BYTE* const oend = ostart + dstSize;
	BYTE* op = ostart;

	BIT_CStream_t bitC;

	if (dstSize < 8) return 0;
	{ size_t const initErr = BIT_initCStream(&bitC, op, oend - op);
	if (HUF_isError(initErr)) return 0; }
	size_t n = srcSize & ~3;
	switch (srcSize & 3)
	{
	case 3: HUF_encodeSymbol(&bitC, ip[n + 2], CTable);
		HUF_FLUSHBITS_2(&bitC);
	case 2: HUF_encodeSymbol(&bitC, ip[n + 1], CTable);
		HUF_FLUSHBITS_1(&bitC);

	case 1: HUF_encodeSymbol(&bitC, ip[n + 0], CTable);
		HUF_FLUSHBITS(&bitC);

	case 0:
	default: break;
	}
	pre_length = bitC.bitPos;
	ip += n;
	__mmask64 * ptr = (__mmask64*)bitC.ptr;
	for (; n >= 256; n -= 256)
	{
		ip -= 64;
		result[0] = _mm512_loadu_si512(ip);
		ip -= 64;
		result[2] = _mm512_loadu_si512(ip);
		ip -= 64;
		result[4] = _mm512_loadu_si512(ip);
		ip -= 64;
		result[6] = _mm512_loadu_si512(ip);

		result[1] = _mm512_permutexvar_epi16(flag_pad1, result[0]);
		result[3] = _mm512_permutexvar_epi16(flag_pad1, result[2]);
		result[5] = _mm512_permutexvar_epi16(flag_pad1, result[4]);
		result[7] = _mm512_permutexvar_epi16(flag_pad1, result[6]);
		result[1] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[1], flag_shuffle);
		result[3] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[3], flag_shuffle);
		result[5] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[5], flag_shuffle);
		result[7] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[7], flag_shuffle);

		result[0] = _mm512_permutexvar_epi16(flag_pad0, result[0]);
		result[2] = _mm512_permutexvar_epi16(flag_pad0, result[2]);
		result[4] = _mm512_permutexvar_epi16(flag_pad0, result[4]);
		result[6] = _mm512_permutexvar_epi16(flag_pad0, result[6]);
		result[0] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[0], flag_shuffle);
		result[2] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[2], flag_shuffle);
		result[4] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[4], flag_shuffle);
		result[6] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[6], flag_shuffle);

		temp[0] = _mm512_slli_epi16(result[0], 8);
		temp[1] = _mm512_slli_epi16(result[1], 8);
		temp[2] = _mm512_slli_epi16(result[2], 8);
		temp[3] = _mm512_slli_epi16(result[3], 8);		
		temp[4] = _mm512_slli_epi16(result[4], 8);
		temp[5] = _mm512_slli_epi16(result[5], 8);
		temp[6] = _mm512_slli_epi16(result[6], 8);
		temp[7] = _mm512_slli_epi16(result[7], 8);

		k1[0] = _mm512_movepi16_mask(temp[0]);
		k1[1] = _mm512_movepi16_mask(temp[1]);
		k1[2] = _mm512_movepi16_mask(temp[2]);
		k1[3] = _mm512_movepi16_mask(temp[3]);
		k1[4] = _mm512_movepi16_mask(temp[4]);
		k1[5] = _mm512_movepi16_mask(temp[5]);
		k1[6] = _mm512_movepi16_mask(temp[6]);
		k1[7] = _mm512_movepi16_mask(temp[7]);

		k2[0] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[0], 1));
		k2[1] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[1], 1));
		k2[2] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[2], 1));
		k2[3] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[3], 1));
		k2[4] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[4], 1));
		k2[5] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[5], 1));
		k2[6] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[6], 1));
		k2[7] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[7], 1));

		result[0] = _mm512_mask2_permutex2var_epi16(table[0], result[0], ~k1[0] & ~k2[0], table[1]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[0], result[1], ~k1[1] & ~k2[1], table[1]);
		result[2] = _mm512_mask2_permutex2var_epi16(table[0], result[2], ~k1[2] & ~k2[2], table[1]);
		result[3] = _mm512_mask2_permutex2var_epi16(table[0], result[3], ~k1[3] & ~k2[3], table[1]);
		result[4] = _mm512_mask2_permutex2var_epi16(table[0], result[4], ~k1[4] & ~k2[4], table[1]);
		result[5] = _mm512_mask2_permutex2var_epi16(table[0], result[5], ~k1[5] & ~k2[5], table[1]);
		result[6] = _mm512_mask2_permutex2var_epi16(table[0], result[6], ~k1[6] & ~k2[6], table[1]);
		result[7] = _mm512_mask2_permutex2var_epi16(table[0], result[7], ~k1[7] & ~k2[7], table[1]);

		result[0] = _mm512_mask2_permutex2var_epi16(table[2], result[0], ~k1[0] & k2[0], table[3]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[2], result[1], ~k1[1] & k2[1], table[3]);
		result[2] = _mm512_mask2_permutex2var_epi16(table[2], result[2], ~k1[2] & k2[2], table[3]);
		result[3] = _mm512_mask2_permutex2var_epi16(table[2], result[3], ~k1[3] & k2[3], table[3]);
		result[4] = _mm512_mask2_permutex2var_epi16(table[2], result[4], ~k1[4] & k2[4], table[3]);
		result[5] = _mm512_mask2_permutex2var_epi16(table[2], result[5], ~k1[5] & k2[5], table[3]);
		result[6] = _mm512_mask2_permutex2var_epi16(table[2], result[6], ~k1[6] & k2[6], table[3]);
		result[7] = _mm512_mask2_permutex2var_epi16(table[2], result[7], ~k1[7] & k2[7], table[3]);

		result[0] = _mm512_mask2_permutex2var_epi16(table[4], result[0], k1[0] & ~k2[0], table[5]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[4], result[1], k1[1] & ~k2[1], table[5]);
		result[2] = _mm512_mask2_permutex2var_epi16(table[4], result[2], k1[2] & ~k2[2], table[5]);
		result[3] = _mm512_mask2_permutex2var_epi16(table[4], result[3], k1[3] & ~k2[3], table[5]);
		result[4] = _mm512_mask2_permutex2var_epi16(table[4], result[4], k1[4] & ~k2[4], table[5]);
		result[5] = _mm512_mask2_permutex2var_epi16(table[4], result[5], k1[5] & ~k2[5], table[5]);
		result[6] = _mm512_mask2_permutex2var_epi16(table[4], result[6], k1[6] & ~k2[6], table[5]);
		result[7] = _mm512_mask2_permutex2var_epi16(table[4], result[7], k1[7] & ~k2[7], table[5]);

		result[0] = _mm512_mask2_permutex2var_epi16(table[6], result[0], k1[0] & k2[0], table[7]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[6], result[1], k1[1] & k2[1], table[7]);
		result[2] = _mm512_mask2_permutex2var_epi16(table[6], result[2], k1[2] & k2[2], table[7]);
		result[3] = _mm512_mask2_permutex2var_epi16(table[6], result[3], k1[3] & k2[3], table[7]);
		result[4] = _mm512_mask2_permutex2var_epi16(table[6], result[4], k1[4] & k2[4], table[7]);
		result[5] = _mm512_mask2_permutex2var_epi16(table[6], result[5], k1[5] & k2[5], table[7]);
		result[6] = _mm512_mask2_permutex2var_epi16(table[6], result[6], k1[6] & k2[6], table[7]);
		result[7] = _mm512_mask2_permutex2var_epi16(table[6], result[7], k1[7] & k2[7], table[7]);
		//Shift and merge
		idx_length_16_odd[0] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[0]), one);
		idx_length_16_odd[1] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[1]), one);
		idx_length_16_odd[2] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[2]), one);
		idx_length_16_odd[3] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[3]), one);
		idx_length_16_odd[4] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[4]), one);
		idx_length_16_odd[5] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[5]), one);
		idx_length_16_odd[6] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[6]), one);
		idx_length_16_odd[7] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[7]), one);

		idx_length_16_even[0] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[0], 16)), one);
		idx_length_16_even[1] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[1], 16)), one);
		idx_length_16_even[2] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[2], 16)), one);
		idx_length_16_even[3] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[3], 16)), one);
		idx_length_16_even[4] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[4], 16)), one);
		idx_length_16_even[5] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[5], 16)), one);
		idx_length_16_even[6] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[6], 16)), one);
		idx_length_16_even[7] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[7], 16)), one);

		idx_length_32[0] = _mm512_add_epi32(idx_length_16_even[0], idx_length_16_odd[0]);
		idx_length_32[1] = _mm512_add_epi32(idx_length_16_even[1], idx_length_16_odd[1]);
		idx_length_32[2] = _mm512_add_epi32(idx_length_16_even[2], idx_length_16_odd[2]);
		idx_length_32[3] = _mm512_add_epi32(idx_length_16_even[3], idx_length_16_odd[3]);
		idx_length_32[4] = _mm512_add_epi32(idx_length_16_even[4], idx_length_16_odd[4]);
		idx_length_32[5] = _mm512_add_epi32(idx_length_16_even[5], idx_length_16_odd[5]);
		idx_length_32[6] = _mm512_add_epi32(idx_length_16_even[6], idx_length_16_odd[6]);
		idx_length_32[7] = _mm512_add_epi32(idx_length_16_even[7], idx_length_16_odd[7]);

		result[0] = _mm512_sllv_epi16(result[0], idx_length_16_even[0]);
		result[1] = _mm512_sllv_epi16(result[1], idx_length_16_even[1]);
		result[2] = _mm512_sllv_epi16(result[2], idx_length_16_even[2]);
		result[3] = _mm512_sllv_epi16(result[3], idx_length_16_even[3]);
		result[4] = _mm512_sllv_epi16(result[4], idx_length_16_even[4]);
		result[5] = _mm512_sllv_epi16(result[5], idx_length_16_even[5]);
		result[6] = _mm512_sllv_epi16(result[6], idx_length_16_even[6]);
		result[7] = _mm512_sllv_epi16(result[7], idx_length_16_even[7]);

		idx_length_16_odd[0] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[0]);
		idx_length_16_odd[1] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[1]);
		idx_length_16_odd[2] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[2]);
		idx_length_16_odd[3] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[3]);
		idx_length_16_odd[4] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[4]);
		idx_length_16_odd[5] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[5]);
		idx_length_16_odd[6] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[6]);
		idx_length_16_odd[7] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[7]);

		idx_length_16_even[0] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[0]);
		idx_length_16_even[1] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[1]);
		idx_length_16_even[2] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[2]);
		idx_length_16_even[3] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[3]);
		idx_length_16_even[4] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[4]);
		idx_length_16_even[5] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[5]);
		idx_length_16_even[6] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[6]);
		idx_length_16_even[7] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[7]);

		result[0] = _mm512_sllv_epi32(result[0], idx_length_16_odd[0]);
		result[1] = _mm512_sllv_epi32(result[1], idx_length_16_odd[1]);
		result[2] = _mm512_sllv_epi32(result[2], idx_length_16_odd[2]);
		result[3] = _mm512_sllv_epi32(result[3], idx_length_16_odd[3]);
		result[4] = _mm512_sllv_epi32(result[4], idx_length_16_odd[4]);
		result[5] = _mm512_sllv_epi32(result[5], idx_length_16_odd[5]);
		result[6] = _mm512_sllv_epi32(result[6], idx_length_16_odd[6]);
		result[7] = _mm512_sllv_epi32(result[7], idx_length_16_odd[7]);

		result[0] = _mm512_srlv_epi32(result[0], idx_length_16_even[0]);
		result[1] = _mm512_srlv_epi32(result[1], idx_length_16_even[1]);
		result[2] = _mm512_srlv_epi32(result[2], idx_length_16_even[2]);
		result[3] = _mm512_srlv_epi32(result[3], idx_length_16_even[3]);
		result[4] = _mm512_srlv_epi32(result[4], idx_length_16_even[4]);
		result[5] = _mm512_srlv_epi32(result[5], idx_length_16_even[5]);
		result[6] = _mm512_srlv_epi32(result[6], idx_length_16_even[6]);
		result[7] = _mm512_srlv_epi32(result[7], idx_length_16_even[7]);

		ls_64[0] = _mm512_srli_epi64(idx_length_32[0], 32);
		ls_64[1] = _mm512_srli_epi64(idx_length_32[1], 32);
		ls_64[2] = _mm512_srli_epi64(idx_length_32[2], 32);
		ls_64[3] = _mm512_srli_epi64(idx_length_32[3], 32);
		ls_64[4] = _mm512_srli_epi64(idx_length_32[4], 32);
		ls_64[5] = _mm512_srli_epi64(idx_length_32[5], 32);
		ls_64[6] = _mm512_srli_epi64(idx_length_32[6], 32);
		ls_64[7] = _mm512_srli_epi64(idx_length_32[7], 32);

		result[0] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[0], ls_64[0]), ls_64[0]);
		result[1] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[1], ls_64[1]), ls_64[1]);
		result[2] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[2], ls_64[2]), ls_64[2]);
		result[3] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[3], ls_64[3]), ls_64[3]);
		result[4] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[4], ls_64[4]), ls_64[4]);
		result[5] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[5], ls_64[5]), ls_64[5]);
		result[6] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[6], ls_64[6]), ls_64[6]);
		result[7] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[7], ls_64[7]), ls_64[7]);

		ls_64[0] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[0]);
		ls_64[1] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[1]);
		ls_64[2] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[2]);
		ls_64[3] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[3]);
		ls_64[4] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[4]);
		ls_64[5] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[5]);
		ls_64[6] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[6]);
		ls_64[7] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[7]);

		result[0] = _mm512_srlv_epi64(result[0], ls_64[0]);
		result[1] = _mm512_srlv_epi64(result[1], ls_64[1]);
		result[2] = _mm512_srlv_epi64(result[2], ls_64[2]);
		result[3] = _mm512_srlv_epi64(result[3], ls_64[3]);
		result[4] = _mm512_srlv_epi64(result[4], ls_64[4]);
		result[5] = _mm512_srlv_epi64(result[5], ls_64[5]);
		result[6] = _mm512_srlv_epi64(result[6], ls_64[6]);
		result[7] = _mm512_srlv_epi64(result[7], ls_64[7]);

		idx_length_64[0] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[0], _mm512_slli_epi64(idx_length_32[0], 32)), 32));
		idx_length_64[1] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[1], _mm512_slli_epi64(idx_length_32[1], 32)), 32));
		idx_length_64[2] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[2], _mm512_slli_epi64(idx_length_32[2], 32)), 32));
		idx_length_64[3] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[3], _mm512_slli_epi64(idx_length_32[3], 32)), 32));
		idx_length_64[4] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[4], _mm512_slli_epi64(idx_length_32[4], 32)), 32));
		idx_length_64[5] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[5], _mm512_slli_epi64(idx_length_32[5], 32)), 32));
		idx_length_64[6] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[6], _mm512_slli_epi64(idx_length_32[6], 32)), 32));
		idx_length_64[7] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[7], _mm512_slli_epi64(idx_length_32[7], 32)), 32));

		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[1]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[1]);
		for (int i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];
			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[0]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[0]);
		for (int i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];

			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[3]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[3]);
		for (int i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];
			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[2]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[2]);
		for (int i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];

			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[5]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[5]);
		for (int i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];
			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[4]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[4]);
		for (int i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];

			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[7]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[7]);
		for (int i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];
			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[6]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[6]);
		for (int i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];

			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
	}
	//n<256 and n>128
	if(n >= 128)
	{
		n -= 128;
		ip -= 64;
		result[0] = _mm512_loadu_si512(ip);
		ip -= 64;
		result[2] = _mm512_loadu_si512(ip);

		result[1] = _mm512_permutexvar_epi16(flag_pad1, result[0]);
		result[3] = _mm512_permutexvar_epi16(flag_pad1, result[2]);
		result[1] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[1], flag_shuffle);
		result[3] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[3], flag_shuffle);

		result[0] = _mm512_permutexvar_epi16(flag_pad0, result[0]);
		result[2] = _mm512_permutexvar_epi16(flag_pad0, result[2]);
		result[0] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[0], flag_shuffle);
		result[2] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[2], flag_shuffle);

		temp[0] = _mm512_slli_epi16(result[0], 8);
		temp[1] = _mm512_slli_epi16(result[1], 8);
		temp[2] = _mm512_slli_epi16(result[2], 8);
		temp[3] = _mm512_slli_epi16(result[3], 8);

		k1[0] = _mm512_movepi16_mask(temp[0]);
		k1[1] = _mm512_movepi16_mask(temp[1]);
		k1[2] = _mm512_movepi16_mask(temp[2]);
		k1[3] = _mm512_movepi16_mask(temp[3]);

		k2[0] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[0], 1));
		k2[1] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[1], 1));
		k2[2] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[2], 1));
		k2[3] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[3], 1));

		result[0] = _mm512_mask2_permutex2var_epi16(table[0], result[0], ~k1[0] & ~k2[0], table[1]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[0], result[1], ~k1[1] & ~k2[1], table[1]);
		result[2] = _mm512_mask2_permutex2var_epi16(table[0], result[2], ~k1[2] & ~k2[2], table[1]);
		result[3] = _mm512_mask2_permutex2var_epi16(table[0], result[3], ~k1[3] & ~k2[3], table[1]);

		result[0] = _mm512_mask2_permutex2var_epi16(table[2], result[0], ~k1[0] & k2[0], table[3]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[2], result[1], ~k1[1] & k2[1], table[3]);
		result[2] = _mm512_mask2_permutex2var_epi16(table[2], result[2], ~k1[2] & k2[2], table[3]);
		result[3] = _mm512_mask2_permutex2var_epi16(table[2], result[3], ~k1[3] & k2[3], table[3]);

		result[0] = _mm512_mask2_permutex2var_epi16(table[4], result[0], k1[0] & ~k2[0], table[5]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[4], result[1], k1[1] & ~k2[1], table[5]);
		result[2] = _mm512_mask2_permutex2var_epi16(table[4], result[2], k1[2] & ~k2[2], table[5]);
		result[3] = _mm512_mask2_permutex2var_epi16(table[4], result[3], k1[3] & ~k2[3], table[5]);

		result[0] = _mm512_mask2_permutex2var_epi16(table[6], result[0], k1[0] & k2[0], table[7]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[6], result[1], k1[1] & k2[1], table[7]);
		result[2] = _mm512_mask2_permutex2var_epi16(table[6], result[2], k1[2] & k2[2], table[7]);
		result[3] = _mm512_mask2_permutex2var_epi16(table[6], result[3], k1[3] & k2[3], table[7]);
		//Shift and merge
		idx_length_16_odd[0] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[0]), one);
		idx_length_16_odd[1] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[1]), one);
		idx_length_16_odd[2] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[2]), one);
		idx_length_16_odd[3] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[3]), one);

		idx_length_16_even[0] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[0], 16)), one);
		idx_length_16_even[1] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[1], 16)), one);
		idx_length_16_even[2] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[2], 16)), one);
		idx_length_16_even[3] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[3], 16)), one);

		idx_length_32[0] = _mm512_add_epi32(idx_length_16_even[0], idx_length_16_odd[0]);
		idx_length_32[1] = _mm512_add_epi32(idx_length_16_even[1], idx_length_16_odd[1]);
		idx_length_32[2] = _mm512_add_epi32(idx_length_16_even[2], idx_length_16_odd[2]);
		idx_length_32[3] = _mm512_add_epi32(idx_length_16_even[3], idx_length_16_odd[3]);

		result[0] = _mm512_sllv_epi16(result[0], idx_length_16_even[0]);
		result[1] = _mm512_sllv_epi16(result[1], idx_length_16_even[1]);
		result[2] = _mm512_sllv_epi16(result[2], idx_length_16_even[2]);
		result[3] = _mm512_sllv_epi16(result[3], idx_length_16_even[3]);

		idx_length_16_odd[0] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[0]);
		idx_length_16_odd[1] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[1]);
		idx_length_16_odd[2] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[2]);
		idx_length_16_odd[3] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[3]);

		idx_length_16_even[0] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[0]);
		idx_length_16_even[1] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[1]);
		idx_length_16_even[2] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[2]);
		idx_length_16_even[3] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[3]);

		result[0] = _mm512_sllv_epi32(result[0], idx_length_16_odd[0]);
		result[1] = _mm512_sllv_epi32(result[1], idx_length_16_odd[1]);
		result[2] = _mm512_sllv_epi32(result[2], idx_length_16_odd[2]);
		result[3] = _mm512_sllv_epi32(result[3], idx_length_16_odd[3]);

		result[0] = _mm512_srlv_epi32(result[0], idx_length_16_even[0]);
		result[1] = _mm512_srlv_epi32(result[1], idx_length_16_even[1]);
		result[2] = _mm512_srlv_epi32(result[2], idx_length_16_even[2]);
		result[3] = _mm512_srlv_epi32(result[3], idx_length_16_even[3]);

		ls_64[0] = _mm512_srli_epi64(idx_length_32[0], 32);
		ls_64[1] = _mm512_srli_epi64(idx_length_32[1], 32);
		ls_64[2] = _mm512_srli_epi64(idx_length_32[2], 32);
		ls_64[3] = _mm512_srli_epi64(idx_length_32[3], 32);

		result[0] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[0], ls_64[0]), ls_64[0]);
		result[1] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[1], ls_64[1]), ls_64[1]);
		result[2] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[2], ls_64[2]), ls_64[2]);
		result[3] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[3], ls_64[3]), ls_64[3]);

		ls_64[0] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[0]);
		ls_64[1] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[1]);
		ls_64[2] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[2]);
		ls_64[3] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[3]);

		result[0] = _mm512_srlv_epi64(result[0], ls_64[0]);
		result[1] = _mm512_srlv_epi64(result[1], ls_64[1]);
		result[2] = _mm512_srlv_epi64(result[2], ls_64[2]);
		result[3] = _mm512_srlv_epi64(result[3], ls_64[3]);

		idx_length_64[0] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[0], _mm512_slli_epi64(idx_length_32[0], 32)), 32));
		idx_length_64[1] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[1], _mm512_slli_epi64(idx_length_32[1], 32)), 32));
		idx_length_64[2] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[2], _mm512_slli_epi64(idx_length_32[2], 32)), 32));
		idx_length_64[3] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[3], _mm512_slli_epi64(idx_length_32[3], 32)), 32));

		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[1]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[1]);
		for (int i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];
			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[0]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[0]);
		for (int i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];

			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[3]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[3]);
		for (int i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];
			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[2]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[2]);
		for (int i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];

			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
	}
	if (n >= 64)
	{
		ip -= 64;
		n -= 64;
		result[0] = _mm512_loadu_si512(ip);
		result[1] = _mm512_permutexvar_epi16(flag_pad1, result[0]);
		result[1] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[1], flag_shuffle);
		result[0] = _mm512_permutexvar_epi16(flag_pad0, result[0]);
		result[0] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[0], flag_shuffle);

		temp[0] = _mm512_slli_epi16(result[0], 8);
		temp[1] = _mm512_slli_epi16(result[1], 8);
		k1[0] = _mm512_movepi16_mask(temp[0]);
		k1[1] = _mm512_movepi16_mask(temp[1]);
		k2[0] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[0], 1));
		k2[1] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[1], 1));

		result[0] = _mm512_mask2_permutex2var_epi16(table[0], result[0], ~k1[0] & ~k2[0], table[1]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[0], result[1], ~k1[1] & ~k2[1], table[1]);
		result[0] = _mm512_mask2_permutex2var_epi16(table[2], result[0], ~k1[0] & k2[0], table[3]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[2], result[1], ~k1[1] & k2[1], table[3]);
		result[0] = _mm512_mask2_permutex2var_epi16(table[4], result[0], k1[0] & ~k2[0], table[5]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[4], result[1], k1[1] & ~k2[1], table[5]);
		result[0] = _mm512_mask2_permutex2var_epi16(table[6], result[0], k1[0] & k2[0], table[7]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[6], result[1], k1[1] & k2[1], table[7]);
		//获取移位长度
		//移位	
		idx_length_16_odd[0] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[0]), one);
		idx_length_16_odd[1] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[1]), one);
		idx_length_16_even[0] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[0], 16)), one);
		idx_length_16_even[1] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[1], 16)), one);
		idx_length_32[0] = _mm512_add_epi32(idx_length_16_even[0], idx_length_16_odd[0]);
		idx_length_32[1] = _mm512_add_epi32(idx_length_16_even[1], idx_length_16_odd[1]);
		result[0] = _mm512_sllv_epi16(result[0], idx_length_16_even[0]);//左移
		result[1] = _mm512_sllv_epi16(result[1], idx_length_16_even[1]);//左移
		idx_length_16_odd[0] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[0]);
		idx_length_16_odd[1] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[1]);
		idx_length_16_even[0] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[0]);
		idx_length_16_even[1] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[1]);
		result[0] = _mm512_sllv_epi32(result[0], idx_length_16_odd[0]);//左移
		result[1] = _mm512_sllv_epi32(result[1], idx_length_16_odd[1]);//左移
		result[0] = _mm512_srlv_epi32(result[0], idx_length_16_even[0]);//右移
		result[1] = _mm512_srlv_epi32(result[1], idx_length_16_even[1]);//右移
		ls_64[0] = _mm512_srli_epi64(idx_length_32[0], 32);
		ls_64[1] = _mm512_srli_epi64(idx_length_32[1], 32);
		result[0] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[0], ls_64[0]), ls_64[0]);//左移
		result[1] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[1], ls_64[1]), ls_64[1]);//左移

		ls_64[0] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[0]);;
		ls_64[1] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[1]);;
		result[0] = _mm512_srlv_epi64(result[0], ls_64[0]);
		result[1] = _mm512_srlv_epi64(result[1], ls_64[1]);
		idx_length_64[0] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[0], _mm512_slli_epi64(idx_length_32[0], 32)), 32));
		idx_length_64[1] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[1], _mm512_slli_epi64(idx_length_32[1], 32)), 32));
		//写出数据	
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[1]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[1]);

		for (int i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];
			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[0]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[0]);

		for (int i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];

			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
	}
	if (n > 32 && n < 64)
	{
		BYTE temp_in[64] = { 0 };
		ip = src;
		memcpy(temp_in + 64 - n, ip, n);
		result[0] = _mm512_loadu_si512(temp_in);
		result[1] = _mm512_permutexvar_epi16(flag_pad1, result[0]);
		result[1] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[1], flag_shuffle);
		result[0] = _mm512_permutexvar_epi16(flag_pad0, result[0]);
		result[0] = _mm512_maskz_shuffle_epi8(0x5555555555555555, result[0], flag_shuffle);
		//设置mask,mask前5位标记子表中的位置，第6位标记每一组表中第几个子表，第7-8位标记第几组表。
		//查表
		temp[0] = _mm512_slli_epi16(result[0], 8);
		temp[1] = _mm512_slli_epi16(result[1], 8);
		k1[0] = _mm512_movepi16_mask(temp[0]);
		k1[1] = _mm512_movepi16_mask(temp[1]);
		k2[0] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[0], 1));
		k2[1] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[1], 1));

		result[0] = _mm512_mask2_permutex2var_epi16(table[0], result[0], ~k1[0] & ~k2[0], table[1]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[0], result[1], ~k1[1] & ~k2[1], table[1]);
		result[0] = _mm512_mask2_permutex2var_epi16(table[2], result[0], ~k1[0] & k2[0], table[3]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[2], result[1], ~k1[1] & k2[1], table[3]);
		result[0] = _mm512_mask2_permutex2var_epi16(table[4], result[0], k1[0] & ~k2[0], table[5]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[4], result[1], k1[1] & ~k2[1], table[5]);
		result[0] = _mm512_mask2_permutex2var_epi16(table[6], result[0], k1[0] & k2[0], table[7]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[6], result[1], k1[1] & k2[1], table[7]);

		__mmask32 mask_zero = (__mmask32)((1llu << (n - 32)) - 1) << (64 - n);
		result[0] = _mm512_mask_blend_epi16(mask_zero, zero, result[0]);
		//获取移位长度
		//移位	
		idx_length_16_odd[0] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[0]), one);
		idx_length_16_odd[1] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[1]), one);
		idx_length_16_even[0] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[0], 16)), one);
		idx_length_16_even[1] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[1], 16)), one);
		idx_length_32[0] = _mm512_add_epi32(idx_length_16_even[0], idx_length_16_odd[0]);
		idx_length_32[1] = _mm512_add_epi32(idx_length_16_even[1], idx_length_16_odd[1]);
		result[0] = _mm512_sllv_epi16(result[0], idx_length_16_even[0]);//左移
		result[1] = _mm512_sllv_epi16(result[1], idx_length_16_even[1]);//左移
		idx_length_16_odd[0] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[0]);
		idx_length_16_odd[1] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[1]);
		idx_length_16_even[0] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[0]);
		idx_length_16_even[1] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[1]);
		result[0] = _mm512_sllv_epi32(result[0], idx_length_16_odd[0]);//左移
		result[1] = _mm512_sllv_epi32(result[1], idx_length_16_odd[1]);//左移
		result[0] = _mm512_srlv_epi32(result[0], idx_length_16_even[0]);//右移
		result[1] = _mm512_srlv_epi32(result[1], idx_length_16_even[1]);//右移
		ls_64[0] = _mm512_srli_epi64(idx_length_32[0], 32);
		ls_64[1] = _mm512_srli_epi64(idx_length_32[1], 32);
		result[0] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[0], ls_64[0]), ls_64[0]);//左移
		result[1] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[1], ls_64[1]), ls_64[1]);//左移

		ls_64[0] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[0]);;
		ls_64[1] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[1]);;
		result[0] = _mm512_srlv_epi64(result[0], ls_64[0]);
		result[1] = _mm512_srlv_epi64(result[1], ls_64[1]);
		idx_length_64[0] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[0], _mm512_slli_epi64(idx_length_32[0], 32)), 32));
		idx_length_64[1] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[1], _mm512_slli_epi64(idx_length_32[1], 32)), 32));
		//写出数据	
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[1]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[1]);

		for (int i = 7; i >= 0; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];
			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[0]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[0]);

		for (int i = 7; (i >= 0) && buf_64[i]; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];
			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
	}
	if (n <= 32 && n)
	{

		BYTE temp_in[64] = { 0 };
		ip = src;
		memcpy(temp_in + 64 - n, ip, n);
		result[0] = _mm512_loadu_si512(temp_in);
		result[1] = _mm512_permutexvar_epi16(flag_pad1, result[0]);
		result[1] = _mm512_maskz_shuffle_epi8(0x5555555555555555llu, result[1], flag_shuffle);
		temp[1] = _mm512_slli_epi16(result[1], 8);
		k1[1] = _mm512_movepi16_mask(temp[1]);
		k2[1] = _mm512_movepi16_mask(_mm512_slli_epi16(temp[1], 1));
		result[1] = _mm512_mask2_permutex2var_epi16(table[0], result[1], ~k1[1] & ~k2[1], table[1]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[2], result[1], ~k1[1] & k2[1], table[3]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[4], result[1], k1[1] & ~k2[1], table[5]);
		result[1] = _mm512_mask2_permutex2var_epi16(table[6], result[1], k1[1] & k2[1], table[7]);
		__mmask32 mask_zero = (__mmask32)((1llu << n) - 1) << (32 - n);
		result[1] = _mm512_mask_blend_epi16(mask_zero, zero, result[1]);

		idx_length_16_odd[1] = _mm512_add_epi32(_mm512_lzcnt_epi32(result[1]), one);
		idx_length_16_even[1] = _mm512_add_epi32(_mm512_lzcnt_epi32(_mm512_slli_epi32(result[1], 16)), one);
		idx_length_32[1] = _mm512_add_epi32(idx_length_16_even[1], idx_length_16_odd[1]);
		result[1] = _mm512_sllv_epi16(result[1], idx_length_16_even[1]);//左移
		idx_length_16_odd[1] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_16_odd[1]);
		idx_length_16_even[1] = _mm512_mask_blend_epi32(0xAAAA, zero, idx_length_16_even[1]);
		result[1] = _mm512_sllv_epi32(result[1], idx_length_16_odd[1]);//左移
		result[1] = _mm512_srlv_epi32(result[1], idx_length_16_even[1]);//右移
		ls_64[1] = _mm512_srli_epi64(idx_length_32[1], 32);
		result[1] = _mm512_srlv_epi64(_mm512_sllv_epi64(result[1], ls_64[1]), ls_64[1]);//左移
		ls_64[1] = _mm512_mask_blend_epi32(0x5555, zero, idx_length_32[1]);
		result[1] = _mm512_srlv_epi64(result[1], ls_64[1]);
		idx_length_64[1] = _mm512_sub_epi64(num64, _mm512_srli_epi64(_mm512_add_epi32(idx_length_32[1], _mm512_slli_epi64(idx_length_32[1], 32)), 32));
		//写出数据	
		_mm512_mask_storeu_epi64(buf_64_length, 0xFF, idx_length_64[1]);
		_mm512_mask_storeu_epi64(buf_64, 0xFF, result[1]);
		for (int i = 7; (i >= 0) && buf_64[i]; i--)
		{
			*ptr |= buf_64[i] << pre_length;
			pre_length += buf_64_length[i];
			if (pre_length >= 64)
			{
				ptr++;
				pre_length -= 64;
				*ptr = buf_64[i] >> (buf_64_length[i] - pre_length);
			}
		}
	}
	*ptr |= 1llu << pre_length;
	pre_length++;
	size_t const nbBytes = pre_length >> 3;
	return ((BYTE*)ptr - ostart) + nbBytes + ((pre_length & 7) > 0);
}*/

size_t HUF_compressSIMD_internal(void* dst, size_t dstSize,
	const void* src, size_t srcSize,
	const HUF_CElt* CTable, int bmi2)
{
	size_t const segmentSize = (srcSize + 3) / 4;   /* first 3 segments */
	const BYTE* ip = (const BYTE*)src;
	const BYTE* const iend = ip + srcSize;
	BYTE* const ostart = (BYTE*)dst;
	BYTE* const oend = ostart + dstSize;
	BYTE* op = ostart;
	memset((BYTE*)dst, 0, dstSize);

	__m512i table_simd[8];
	U16 tab[256] = { 0 };
	for (int i = 0; i < 256; i++)
	{
		tab[i] = CTable[i].nbBits ? (1 << CTable[i].nbBits) ^ CTable[i].val : 0;
	}
	for (int i = 0; i < 8; i++)
		table_simd[i] = _mm512_loadu_si512(&tab[i * 32]);

	if (dstSize < 6 + 1 + 1 + 1 + 8) return 0;   /* minimum space to compress successfully */
	if (srcSize < 12) return 0;   /* no saving possible : too small input */
	op += 6;   /* jumpTable */
	{   CHECK_V_F(cSize, (op, oend - op, ip, segmentSize, table_simd, CTable));
	if (cSize == 0) return 0;
	assert(cSize <= 65535);
	MEM_writeLE16(ostart, (U16)cSize);
	op += cSize;
	}
	ip += segmentSize;
	{   CHECK_V_F(cSize, HUF_compressSIMD_internal_body(op, oend - op, ip, segmentSize, table_simd, CTable));
	if (cSize == 0) return 0;
	assert(cSize <= 65535);
	MEM_writeLE16(ostart + 2, (U16)cSize);
	op += cSize;
	}
	ip += segmentSize;
	{   CHECK_V_F(cSize, HUF_compressSIMD_internal_body(op, oend - op, ip, segmentSize, table_simd, CTable));
	if (cSize == 0) return 0;
	assert(cSize <= 65535);
	MEM_writeLE16(ostart + 4, (U16)cSize);
	op += cSize;
	}
	ip += segmentSize;
	{   CHECK_V_F(cSize, HUF_compressSIMD_internal_body(op, oend - op, ip, iend - ip, table_simd, CTable));
	if (cSize == 0) return 0;
	op += cSize;
	}
	return op - ostart;

}
/////////////////////////////////////////////////////////////////////////////////
static size_t
HUF_compress4X_usingCTable_internal(void* dst, size_t dstSize,
	const void* src, size_t srcSize,
	const HUF_CElt* CTable, int bmi2)
{
	size_t const segmentSize = (srcSize + 3) / 4;   /* first 3 segments */
	const BYTE* ip = (const BYTE*)src;
	const BYTE* const iend = ip + srcSize;
	BYTE* const ostart = (BYTE*)dst;
	BYTE* const oend = ostart + dstSize;
	BYTE* op = ostart;


	if (dstSize < 6 + 1 + 1 + 1 + 8) return 0;   /* minimum space to compress successfully */
	if (srcSize < 12) return 0;   /* no saving possible : too small input */
	op += 6;   /* jumpTable */

	{   CHECK_V_F(cSize, HUF_compress1X_usingCTable_internal(op, oend - op, ip, segmentSize, CTable, bmi2));
	if (cSize == 0) return 0;
	assert(cSize <= 65535);
	MEM_writeLE16(ostart, (U16)cSize);
	op += cSize;
	}

	ip += segmentSize;
	{CHECK_V_F(cSize, HUF_compress1X_usingCTable_internal(op, oend - op, ip, segmentSize, CTable, bmi2));
	if (cSize == 0) return 0;
	assert(cSize <= 65535);
	MEM_writeLE16(ostart + 2, (U16)cSize);
	op += cSize;
	}
	ip += segmentSize;
	{   CHECK_V_F(cSize, HUF_compress1X_usingCTable_internal(op, oend - op, ip, segmentSize, CTable, bmi2));
	if (cSize == 0) return 0;
	assert(cSize <= 65535);
	MEM_writeLE16(ostart + 4, (U16)cSize);
	op += cSize;
	}
	ip += segmentSize;
	{   CHECK_V_F(cSize, HUF_compress1X_usingCTable_internal(op, oend - op, ip, iend - ip, CTable, bmi2));
	if (cSize == 0) return 0;
	op += cSize;
	}
	return op - ostart;
}

size_t HUF_compress4X_usingCTable(void* dst, size_t dstSize, const void* src, size_t srcSize, const HUF_CElt* CTable)
{
	return HUF_compress4X_usingCTable_internal(dst, dstSize, src, srcSize, CTable, /* bmi2 */ 0);
}


/*
*
*压缩函数入口，singleStream为1，调用SIMD优化函数进行压缩，否则使用原有函数压缩
*/
static size_t HUF_compressCTable_internal(
	BYTE* const ostart, BYTE* op, BYTE* const oend,
	const void* src, size_t srcSize,
	unsigned singleStream, const HUF_CElt* CTable, const int bmi2)
{
	singleStream = 1;
	/*size_t const cSize = singleStream ?
						 HUF_compress1X_usingCTable_internal(op, oend - op, src, srcSize, CTable, bmi2) :
						 HUF_compress4X_usingCTable_internal(op, oend - op, src, srcSize, CTable,
						 bmi2);*/
	size_t const cSize = singleStream ?
		HUF_compressSIMD_internal(op, oend - op, src, srcSize, CTable, bmi2) :
		HUF_compress4X_usingCTable_internal(op, oend - op, src, srcSize, CTable, bmi2);
	if (HUF_isError(cSize)) { return cSize; }
	if (cSize == 0) { return 0; }   /* uncompressible */
	op += cSize;
	/* check compressibility */
	if ((size_t)(op - ostart) >= srcSize - 1) { return 0; }
	return op - ostart;
}



typedef struct {
	U32 count[HUF_SYMBOLVALUE_MAX + 1];
	HUF_CElt CTable[HUF_SYMBOLVALUE_MAX + 1];
	huffNodeTable nodeTable;
} HUF_compress_tables_t;

/* HUF_compress_internal() :
 * `workSpace` must a table of at least HUF_WORKSPACE_SIZE_U32 unsigned */
static size_t HUF_compress_internal(
	void* dst, size_t dstSize,
	const void* src, size_t srcSize,
	unsigned maxSymbolValue, unsigned huffLog,
	unsigned singleStream,
	void* workSpace, size_t wkspSize,
	HUF_CElt* oldHufTable, HUF_repeat* repeat, int preferRepeat,
	const int bmi2)
{
	HUF_compress_tables_t* const table = (HUF_compress_tables_t*)workSpace;
	BYTE* const ostart = (BYTE*)dst;
	BYTE* const oend = ostart + dstSize;
	BYTE* op = ostart;

	/* checks & inits */
	if (((size_t)workSpace & 3) != 0) return ERROR(GENERIC);  /* must be aligned on 4-bytes boundaries */
	if (wkspSize < sizeof(*table)) return ERROR(workSpace_tooSmall);
	if (!srcSize) return 0;  /* Uncompressed */
	if (!dstSize) return 0;  /* cannot fit anything within dst budget */
	if (srcSize > HUF_BLOCKSIZE_MAX) return ERROR(srcSize_wrong);   /* current block size limit */
	if (huffLog > HUF_TABLELOG_MAX) return ERROR(tableLog_tooLarge);
	if (maxSymbolValue > HUF_SYMBOLVALUE_MAX) return ERROR(maxSymbolValue_tooLarge);
	if (!maxSymbolValue) maxSymbolValue = HUF_SYMBOLVALUE_MAX;
	if (!huffLog) huffLog = HUF_TABLELOG_DEFAULT;

	/* Heuristic : If old table is valid, use it for small inputs */
	if (preferRepeat && repeat && *repeat == HUF_repeat_valid) {
		return HUF_compressCTable_internal(ostart, op, oend,
			src, srcSize,
			singleStream, oldHufTable, bmi2);
	}

	/* Scan input and build symbol stats */
	{   CHECK_V_F(largest, HIST_count_wksp(table->count, &maxSymbolValue, (const BYTE*)src, srcSize, table->count));
	if (largest == srcSize) { *ostart = ((const BYTE*)src)[0]; return 1; }   /* single symbol, rle */
	if (largest <= (srcSize >> 7) + 1) return 0;   /* heuristic : probably not compressible enough */
	}

	/* Check validity of previous table */
	if (repeat
		&& *repeat == HUF_repeat_check
		&& !HUF_validateCTable(oldHufTable, table->count, maxSymbolValue)) {
		*repeat = HUF_repeat_none;
	}
	/* Heuristic : use existing table for small inputs */
	if (preferRepeat && repeat && *repeat != HUF_repeat_none) {
		return HUF_compressCTable_internal(ostart, op, oend,
			src, srcSize,
			singleStream, oldHufTable, bmi2);
	}

	/* Build Huffman Tree */
	huffLog = HUF_optimalTableLog(huffLog, srcSize, maxSymbolValue);
	{   CHECK_V_F(maxBits, HUF_buildCTable_wksp(table->CTable, table->count,
		maxSymbolValue, huffLog,
		table->nodeTable, sizeof(table->nodeTable)));
	huffLog = (U32)maxBits;
	/* Zero unused symbols in CTable, so we can check it for validity */
	memset(table->CTable + (maxSymbolValue + 1), 0,
		sizeof(table->CTable) - ((maxSymbolValue + 1) * sizeof(HUF_CElt)));
	}

	/* Write table description header */
	{   CHECK_V_F(hSize, HUF_writeCTable(op, dstSize, table->CTable, maxSymbolValue, huffLog));
	/* Check if using previous huffman table is beneficial */
	if (repeat && *repeat != HUF_repeat_none) {
		size_t const oldSize = HUF_estimateCompressedSize(oldHufTable, table->count, maxSymbolValue);
		size_t const newSize = HUF_estimateCompressedSize(table->CTable, table->count, maxSymbolValue);
		if (oldSize <= hSize + newSize || hSize + 12 >= srcSize) {
			return HUF_compressCTable_internal(ostart, op, oend,
				src, srcSize,
				singleStream, oldHufTable, bmi2);
		}
	}

	/* Use the new huffman table */
	if (hSize + 12ul >= srcSize) { return 0; }
	op += hSize;
	if (repeat) { *repeat = HUF_repeat_none; }
	if (oldHufTable)
		memcpy(oldHufTable, table->CTable, sizeof(table->CTable));  /* Save new table */
	}

	/************************************************************************/

	/************************************************************************/

	return HUF_compressCTable_internal(ostart, op, oend, src, srcSize, singleStream, table->CTable, bmi2);


}








size_t HUF_compress1X_wksp(void* dst, size_t dstSize,
	const void* src, size_t srcSize,
	unsigned maxSymbolValue, unsigned huffLog,
	void* workSpace, size_t wkspSize)
{
	return HUF_compress_internal(dst, dstSize, src, srcSize,
		maxSymbolValue, huffLog, 1 /*single stream*/,
		workSpace, wkspSize,
		NULL, NULL, 0, 0 /*bmi2*/);
}

size_t HUF_compress1X_repeat(void* dst, size_t dstSize,
	const void* src, size_t srcSize,
	unsigned maxSymbolValue, unsigned huffLog,
	void* workSpace, size_t wkspSize,
	HUF_CElt* hufTable, HUF_repeat* repeat, int preferRepeat, int bmi2)
{
	return HUF_compress_internal(dst, dstSize, src, srcSize,
		maxSymbolValue, huffLog, 1 /*single stream*/,
		workSpace, wkspSize, hufTable,
		repeat, preferRepeat, bmi2);
}

size_t HUF_compress1X(void* dst, size_t dstSize,
	const void* src, size_t srcSize,
	unsigned maxSymbolValue, unsigned huffLog)
{
	unsigned workSpace[HUF_WORKSPACE_SIZE_U32];
	return HUF_compress1X_wksp(dst, dstSize, src, srcSize, maxSymbolValue, huffLog, workSpace, sizeof(workSpace));
}

/* HUF_compress4X_repeat():
 * compress input using 4 streams.
 * provide workspace to generate compression tables */
size_t HUF_compress4X_wksp(void* dst, size_t dstSize,
	const void* src, size_t srcSize,
	unsigned maxSymbolValue, unsigned huffLog,
	void* workSpace, size_t wkspSize)
{
	return HUF_compress_internal(dst, dstSize, src, srcSize,
		maxSymbolValue, huffLog, 1 /*4 streams*/,
		workSpace, wkspSize,
		NULL, NULL, 0, 0 /*bmi2*/);
}

/* HUF_compress4X_repeat():
 * compress input using 4 streams.
 * re-use an existing huffman compression table */
size_t HUF_compress4X_repeat(void* dst, size_t dstSize,
	const void* src, size_t srcSize,
	unsigned maxSymbolValue, unsigned huffLog,
	void* workSpace, size_t wkspSize,
	HUF_CElt* hufTable, HUF_repeat* repeat, int preferRepeat, int bmi2)
{
	return HUF_compress_internal(dst, dstSize, src, srcSize,
		maxSymbolValue, huffLog, 1 /* 4 streams */,
		workSpace, wkspSize,
		hufTable, repeat, preferRepeat, bmi2);
}

size_t HUF_compress2(void* dst, size_t dstSize,
	const void* src, size_t srcSize,
	unsigned maxSymbolValue, unsigned huffLog)
{
	unsigned workSpace[HUF_WORKSPACE_SIZE_U32];
	return HUF_compress4X_wksp(dst, dstSize, src, srcSize, maxSymbolValue, huffLog, workSpace, sizeof(workSpace));
}

size_t HUF_compress(void* dst, size_t maxDstSize, const void* src, size_t srcSize)
{
	return HUF_compress2(dst, maxDstSize, src, srcSize, 255, HUF_TABLELOG_DEFAULT);
}
