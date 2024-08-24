#ifndef lint
static const char RCSid[] = "$Id: abitmapio.cpp,v 2.2 2024/08/24 23:25:24 greg Exp $";
#endif
/*
 *  abitmapio.cpp
 *  panlib
 *
 *  BitMap class file i/o using BMP format.
 *
 *  Created by Greg Ward on 6/30/16.
 *  Copyright 2016 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "abitmap.h"
#include "bmpfile.h"
#include "dmessage.h"

static const char	BitMap1Dmagic[] = "1-D BitMap";
static const char	BitMap2Dmagic[] = "2-D BitMap";

// Function to write a bitmap to a BMP file
bool
WriteBitMap(const ABitMap &bm, const char *fname)
{
	if (!bm.Length() || fname == NULL || !*fname)
		return false;
	if (bm.Length() >= 1L<<31) {
		DMESG(DMCparameter, "BitMap too long to write as one scanline");
		return false;
	}
	BMPHeader *	hdr = BMPmappedHeader(bm.Length(), 1, 16, 2);
	if (hdr == NULL) {
		DMESG(DMCparameter, "Cannot create BMP header");
		return false;
	}
	strncpy(BMPinfo(hdr), BitMap1Dmagic, 16);

	BMPWriter *	wtr = BMPopenOutputFile(fname, hdr);
	if (wtr == NULL) {
		DMESGF(DMCresource, "Cannot open BMP output file '%s'", fname);
		free(hdr);
		return false;
	}
	DMESGF(DMCtrace, "Writing 1-D bitmap to '%s'", fname);

	memset(wtr->scanline, 0, (bm.Length()+7)>>3);
	for (uint32 i = 0; bm.Find(&i); i++)
		wtr->scanline[i>>3] |= 128>>(i&7);

	int	rval = BMPwriteScanline(wtr);

	BMPcloseOutput(wtr);

	if (rval != BIR_OK) {
		DMESG(DMCdata, BMPerrorMessage(rval));
		return false;
	}
	return true;
}

// Function to write a 2-D bitmap to a BMP file
bool
WriteBitMap2(const ABitMap2 &bm2, const char *fname)
{
	if (!bm2.Width() | !bm2.Height() || fname == NULL || !*fname)
		return false;

	BMPHeader *	hdr = BMPmappedHeader(bm2.Width(), bm2.Height(), 16, 2);
	if (hdr == NULL) {
		DMESG(DMCparameter, "Cannot create BMP header");
		return false;
	}
	strncpy(BMPinfo(hdr), BitMap2Dmagic, 16);
	hdr->yIsDown = 1;		// sane scanline ordering

	BMPWriter *	wtr = BMPopenOutputFile(fname, hdr);
	if (wtr == NULL) {
		DMESGF(DMCresource, "Cannot open BMP output file '%s'", fname);
		free(hdr);
		return false;
	}
	DMESGF(DMCtrace, "Writing 2-D bitmap to '%s'", fname);
	int	rval = BIR_OK;
	int	my = 0;
	for (int y = 0; y < bm2.Height(); y++) {
		memset(wtr->scanline, 0, (bm2.Width()+7)>>3);
		while (y < my) {	// write empty scanlines
			if ((rval = BMPwriteScanline(wtr)) != BIR_OK)
				break;
			++y;
		}
		if ((rval != BIR_OK) | (y >= bm2.Height()))
			break;
		for (int x = 0; bm2.Find(&x, &my) && my == y; x++)
			wtr->scanline[x>>3] |= 128>>(x&7);

		if ((rval = BMPwriteScanline(wtr)) != BIR_OK)
			break;
	}
	BMPcloseOutput(wtr);

	if (rval != BIR_OK) {
		DMESG(DMCdata, BMPerrorMessage(rval));
		return false;
	}
	return true;
}

// Function to read a bitmap from a BMP file
bool
ReadBitMap(ABitMap *bmp, const char *fname)
{
	if (bmp == NULL || fname == NULL || !*fname)
		return false;
						// open call reads 1st scanline
	BMPReader *	rdr = BMPopenInputFile(fname);

	if (rdr == NULL) {
		DMESGF(DMCresource, "Cannot open BMP input file '%s'", fname);
		return false;
	}
	if (rdr->hdr->nColors > 2) {
		DMESGF(DMCdata, "BMP input file '%s' is not a bitmap", fname);
		BMPcloseInput(rdr);
		return false;
	}
	if (rdr->hdr->height != 1) {
		DMESGF(DMCdata, "BMP input file '%s' has more than one scan line", fname);
		BMPcloseInput(rdr);
		return false;
	}
	if (rdr->hdr->infoSiz <= 0 ||
			strncmp(BMPinfo(rdr->hdr), BitMap1Dmagic, rdr->hdr->infoSiz))
		DMESGF(DMCwarning, "BMP file '%s' not flagged as 1-D bitmap", fname);

	DASSERT(rdr->yscan == 0);
	DMESGF(DMCtrace, "Reading 1-D bitmap from '%s'", fname);
	bmp->NewBitMap(rdr->hdr->width);
	for (uint32 i = rdr->hdr->width; i--; )
		if (rdr->scanline[i>>3] & 128>>(i&7))
			bmp->Set(i);
	if (rdr->hdr->palette[0].g >= 128)
		bmp->Invert();		// against normal convention!
	BMPcloseInput(rdr);
	return true;
}

// Function to read a 2-D bitmap from a BMP file
bool
ReadBitMap2(ABitMap2 *bm2p, const char *fname)
{
	if (bm2p == NULL || fname == NULL || !*fname)
		return false;

	BMPReader *	rdr = BMPopenInputFile(fname);

	if (rdr == NULL) {
		DMESGF(DMCresource, "Cannot open BMP input file '%s'", fname);
		return false;
	}
	if (rdr->hdr->nColors > 2) {
		DMESGF(DMCdata, "BMP input file '%s' is not a bitmap", fname);
		BMPcloseInput(rdr);
		return false;
	}
	if (rdr->hdr->infoSiz <= 0 ||
			strncmp(BMPinfo(rdr->hdr), BitMap2Dmagic, rdr->hdr->infoSiz))
		DMESGF(DMCwarning, "BMP file '%s' not flagged as 2-D bitmap", fname);

	int	y, rval = BIR_OK;
	bm2p->NewBitMap(rdr->hdr->width, rdr->hdr->height);
	DMESGF(DMCtrace, "Reading 2-D bitmap from '%s'", fname);
	for (y = 0; y < rdr->hdr->height; y++) {
		const int	my = rdr->hdr->yIsDown ? y
						: rdr->hdr->height-1 - y;
		while (rdr->yscan < y && (rval = BMPreadScanline(rdr)) == BIR_OK)
			;
		if (rval != BIR_OK)
			break;
		for (uint32 x = rdr->hdr->width; x--; )
			if (rdr->scanline[x>>3] & 128>>(x&7))
				bm2p->Set(x, my);
	}
	if (rval == BIR_OK && rdr->hdr->palette[0].g >= 128)
		bm2p->Invert();		// against normal convention!
	BMPcloseInput(rdr);

	if (rval != BIR_OK) {
		bm2p->NewBitMap(0,0);
		sprintf(dmessage_buf, "%s at y==%d", BMPerrorMessage(rval), y);
		DMESG(DMCdata, dmessage_buf);
		return false;
	}
	return true;
}
