/* RCSid $Id: abitmap.h,v 2.2 2024/08/24 23:25:24 greg Exp $ */
/*
 *  abitmap.h
 *  panlib
 *
 *  General bitmap class (mostly inline)
 *
 *  Created by gward on Tue May 15 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */
#ifndef _ABITMAP_H_
#define _ABITMAP_H_

#ifndef _TIFF_
#include "tiff.h"		/* needed for uint32 type */
#endif

#define	ABMend		0xffffffff	// terminal return

// Inline bitmap class
class ABitMap {
private:
	uint32 *	bmap;		// bitmap storage
	uint32		len;		// bitmap size
public:
			ABitMap() { len=0; bmap=0; }
			ABitMap(uint32 n, bool clrset=false) {
				bmap=0; len=0; NewBitMap(n,clrset);
			}
			ABitMap(const ABitMap &orig) {
				bmap=0; len=0; *this=orig;
			}
			~ABitMap() {
				delete [] bmap;
			}
			// Access to raw word data
	uint32 *	base() {
				return bmap;
			}
			// Read-only word access
	const uint32 *	base() const {
				return bmap;
			}
			// Number of words from #bits
	static int32	bmlen(uint32 nbits) {
				return (nbits+0x1f)>>5;
			}
			// Total number of words
	int32		bmlen() const {
				return bmlen(len);
			}
			// Reallocate and clear bitmap
	bool		NewBitMap(uint32 n, bool clrset=false);
			// Clear bitmap to all 0's or 1's
	void		ClearBitMap(bool clrset=false);
			// Steal another bitmap's contents
	bool		Take(ABitMap *srcp) {
				if (srcp == this) return false;
				delete [] bmap; bmap=0; len=0;
				if (!srcp) return false;
				bmap=srcp->bmap; len=srcp->len;
				srcp->bmap=0; srcp->len=0;
				return len;
			}
			// Return number of bits in bitmap
	uint32		Length() const {
				return len;
			}
			// Get original length if RLE (or 0)
	uint32		RLength() const;
			// Compress into a run-length encoded bitmap
	bool		GetRLE(ABitMap *rlep) const;
			// Reconstitute bits from RLE encoding
	bool		SetFromRLE(const ABitMap &rle);
			// Extract bitmap section (true if some overlap)
	bool		GetBits(ABitMap *dp, uint32 i) const;
			// Extract bitmap section (fill if past end)
	ABitMap		GetBits(uint32 i, uint32 n, bool fill=false) const {
				ABitMap	bm(n, fill);
				GetBits(&bm, i);
				return bm;
			}
			// Overlay bitmap section (ignores bits past end)
	bool		AssignBits(uint32 i, const ABitMap &src);
			// Apply operation to bitmap section
	bool		OpBits(uint32 i, char op, const ABitMap &src);
			// Clear bitmap section
	void		ClearBits(uint32 i, uint32 n, bool clrset=false);
			// Clear bits set in second bitmap
	bool		ClearBitsFrom(const ABitMap &src);
			// Count number of bits set in bitmap
	uint32		SumTotal(bool bit2cnt = true) const;
			// Return the next bit position matching val (or ABMend)
	uint32		Find(uint32 i = 0, bool val = true) const;
	uint32		Find(int i, bool val = true) const {
				return Find((uint32)i, val);
			}
	uint32		Find(long i, bool val = true) const {
				return Find((uint32)i, val);
			}
			// Different interface to find the next bit matching val
	bool		Find(uint32 *ip, bool val = true) const {
				if (!ip) return false;
				if (*ip >= len) return false;
				*ip = Find(*ip, val);
				return (*ip < len);
			}
			// Access word for i'th bit
	uint32 &	Word(uint32 i) {
				static uint32	dummy;
				if (i >= len) return dummy;
				return bmap[i>>5];
			}
			// Get word value for i'th bit (0 if out of range)
	uint32		WordV(uint32 i) const {
				if (i >= len) return 0;
				return bmap[i>>5];
			}
			// Index corresponding bit in word
	static uint32	Bit(uint32 i) {
				return 1 << (i & 0x1f);
			}
			// Return value of i'th bit in bitmap
	bool		Check(uint32 i) const {
				return (WordV(i) & Bit(i));
			}
			// Set i'th bit
	void		Set(uint32 i) {
				Word(i) |= Bit(i);
			}
			// Reset i'th bit
	void		Reset(uint32 i) {
				Word(i) &= ~Bit(i);
			}
			// Set i'th bit explicitly on or off
	void		Set(uint32 i, bool switchon) {
				uint32		b = Bit(i);
				uint32 &	w = Word(i);
				if (switchon) w |= b;
				else w &= ~b;
			}
			// Toggle i'th bit
	void		Toggle(uint32 i) {
				Word(i) ^= Bit(i);
			}
			// Set i'th bit if it's clear (or fail)
	bool		TestAndSet(uint32 i) {
				if (i >= len) return false;
				uint32		b = Bit(i);
				uint32 &	w = Word(i);
				uint32		wasset = w & b;
				w |= b;
				return !wasset;
			}
			// Clear i'th bit if it's set (or fail)
	bool		TestAndReset(uint32 i) {
				if (i >= len) return false;
				uint32		b = Bit(i);
				uint32 &	w = Word(i);
				uint32		wasset = w & b;
				w &= ~b;
				return wasset;
			}
			// Set i'th bit on/off (fail if no change)
	bool		TestAndSet(uint32 i, bool switchon) {
				return switchon ? TestAndSet(i) : TestAndReset(i);
			}
			// Invert the entire bitmap
	void		Invert();
			// Downward shift operator, zero fill
	ABitMap &	operator>>=(uint32 nbits);
			// Upward shift operator, zero fill
	ABitMap &	operator<<=(uint32 nbits);
			// Copy operator
	ABitMap &	operator=(const ABitMap &src);
			// Bitwise OR-copy operator
	ABitMap &	operator|=(const ABitMap &src);
			// Bitwise AND-assign operator
	ABitMap &	operator&=(const ABitMap &src);
			// Bitwise XOR-assign operator
	ABitMap &	operator^=(const ABitMap &src);
			// Subtraction operator, synonym for ClearBitsFrom()
	ABitMap &	operator-=(const ABitMap &src) {
				ClearBitsFrom(src);
				return *this;
			}
			// Compare two bitmaps for equality
	bool		operator==(const ABitMap &that) const;
};

inline bool
operator!=(const ABitMap &bm1, const ABitMap &bm2)
{
	return !(bm1 == bm2);
}

inline ABitMap
operator>>(ABitMap bmLeft, uint32 nbr)
{
	return bmLeft >>= nbr;
}

inline ABitMap
operator<<(ABitMap bmLeft, uint32 nbl)
{
	return bmLeft <<= nbl;
}

inline ABitMap
operator|(ABitMap bmLeft, const ABitMap &bmRight)
{
	return bmLeft |= bmRight;
}

inline ABitMap
operator&(ABitMap bmLeft, const ABitMap &bmRight)
{
	return bmLeft &= bmRight;
}

inline ABitMap
operator^(ABitMap bmLeft, const ABitMap &bmRight)
{
	return bmLeft ^= bmRight;
}

inline ABitMap
operator-(ABitMap bmLeft, const ABitMap &bmRight)
{
	return bmLeft -= bmRight;
}

inline ABitMap
operator~(ABitMap bmUnary)
{
	bmUnary.Invert();
	return bmUnary;
}

// 2-dimensional bitmap class
class ABitMap2 : protected ABitMap {
private:
	int		width, height;		// bitmap dimensions
protected:
	uint32		bmi(int x, int y) const {
				if (OffBitMap(x, y)) return ABMend;
				return (uint32)y*width + x;
			}
public:
			ABitMap2() { width=height=0; }
			ABitMap2(int w, int h, bool clrset=false) :
					ABitMap((w>0)&(h>0)&&(w <= ABMend/h)
						? (uint32)w*h : (uint32)0, clrset) {
				if (Length()) {
					width=w; height=h;
				} else { width=height=0; }
			}
			ABitMap2(const ABitMap2 &orig) {
				*this=orig;
			}
			// Access to raw word data
	uint32 *	base() {
				return ABitMap::base();
			}
			// Read-only word access
	const uint32 *	base() const {
				return ABitMap::base();
			}
			// Total number of words
	int32		bmlen() const {
				return ABitMap::bmlen();
			}
			// Return bitmap width
	int		Width() const {
				return width;
			}
			// Return bitmap height
	int		Height() const {
				return height;
			}
			// Is the indicated bit off our bitmap?
	bool		OffBitMap(int x, int y) const {
				return ((x < 0) | (x >= width) |
					(y < 0) | (y >= height));
			}
			// Count number of bits set in bitmap
	uint32		SumTotal(bool bit2cnt = true) const {
				return ABitMap::SumTotal(bit2cnt);
			}
			// Reallocate and clear bitmap
	bool		NewBitMap(int w, int h, bool clrset=false) {
				if ((w <= 0) | (h <= 0) || w > ABMend/h)
					w = h = 0;
				width=w; height=h;
				return ABitMap::NewBitMap((uint32)w*h, clrset);
			}
			// Clear bitmap to all 0's or 1's
	void		ClearBitMap(bool clrset=false) {
				ABitMap::ClearBitMap(clrset);
			}
			// Compress with run-length encoding into a 1-D bitmap
	bool		GetRLE(ABitMap *rlep) const {
				return ABitMap::GetRLE(rlep);
			}
			// Reconstitute bits from RLE encoding (size must match)
	bool		SetFromRLE(int w, int h, const ABitMap &rle);
			// Steal another bitmap's contents
	bool		Take(ABitMap2 *srcp) {
				if (srcp == this) return false;
				width=height=0;
				if (!ABitMap::Take(srcp)) return false;
				width=srcp->width; height=srcp->height;
				srcp->width=srcp->height=0;
				return true;
			}
			// Extract bitmap section (true if some overlap)
	bool		GetRect(ABitMap2 *dp, int sx, int sy) const;
			// Extract bitmap section (fill outside overlap)
	ABitMap2	GetRect(int sx, int sy, int w, int h, bool fill=false) const {
				ABitMap2	bm2(w, h, fill);
				GetRect(&bm2, sx, sy);
				return bm2;
			}
			// Assign bitmap section (ignores anything past edges)
	bool		AssignRect(int dx, int dy, const ABitMap2 &src);
			// Apply operation to bitmap section
	bool		OpRect(int dx, int dy, char op, const ABitMap2 &src);
			// Clear a rectangle
	void		ClearRect(int x, int y, int w, int h, bool clrset=false);
			// Find the next bit matching val (scanline order)
	bool		Find(int *xp, int *yp, bool val=true) const {
				if (!xp | !yp) return false;
				if (width <= 0) return false;
				if ((*xp < 0) | (*yp < 0)) *xp = *yp = 0;
				else if (*xp >= width) { *xp=0; ++(*yp); }
				uint32	i = ABitMap::Find(bmi(*xp,*yp), val);
				if (i == ABMend) { *yp = height; return false; }
				*yp = int(i / width);
				*xp = int(i - *yp*width);
				return true;
			}
			// Get bounds of assigned region
	bool		GetBoundRect(int xymin[2], int wh[2], bool val=true) const;
			// Return value of bit in bitmap
	bool		Check(int x, int y) const {
				return ABitMap::Check(bmi(x,y));
			}
			// Set bit
	void		Set(int x, int y) {
				ABitMap::Set(bmi(x,y));
			}
			// Reset bit
	void		Reset(int x, int y) {
				ABitMap::Reset(bmi(x,y));
			}
			// Set bit explicitly on or off
	void		Set(int x, int y, bool switchon) {
				ABitMap::Set(bmi(x,y), switchon);
			}
			// Toggle bit
	void		Toggle(int x, int y) {
				ABitMap::Toggle(bmi(x,y));
			}
			// Set bit if it's clear (or fail)
	bool		TestAndSet(int x, int y) {
				return ABitMap::TestAndSet(bmi(x,y));
			}
			// Clear bit if it's set (or fail)
	bool		TestAndReset(int x, int y) {
				return ABitMap::TestAndReset(bmi(x,y));
			}
			// Set bit on/off (fail if no change)
	bool		TestAndSet(int x, int y, bool switchon) {
				return ABitMap::TestAndSet(bmi(x,y), switchon);
			}
			// Invert the entire bitmap
	void		Invert() {
				ABitMap::Invert();
			}
			// Shift bitmap, filling uncovered area as indicated
	void		Shift(int dx, int dy, int fill=0);
			// Dilate (or erode) selection by given radius
	void		Expand(double rad, bool val=true);
			// Clear bits set in second map
	bool		ClearBitsFrom(const ABitMap2 &src) {
				if (width != src.width)
					return false;
				return ABitMap::ClearBitsFrom(src);
			}
			// Copy operator
	ABitMap2 &	operator=(const ABitMap2 &src) {
				ABitMap::operator=(src);
				width=src.width; height=src.height;
				return *this;
			}
			// Bitwise OR-copy operator
	ABitMap2 &	operator|=(const ABitMap2 &src) {
				ABitMap::operator|=(src);
				width=src.width; height=src.height;
				return *this;
			}
			// Bitwise AND-assign operator
	ABitMap2 &	operator&=(const ABitMap2 &src) {
				if ((width!=src.width)|(height!=src.height))
					return *this;
				ABitMap::operator&=(src);
				return *this;
			}
			// Bitwise XOR-assign operator
	ABitMap2 &	operator^=(const ABitMap2 &src) {
				if ((width!=src.width)|(height!=src.height))
					return *this;
				ABitMap::operator^=(src);
				return *this;
			}
			// Subtraction operator, synonym for ClearBitsFrom()
	ABitMap2 &	operator-=(const ABitMap2 &src) {
				ClearBitsFrom(src);
				return *this;
			}
			// Compare two bitmaps for equality
	bool		operator==(const ABitMap2 &that) const {
				if (width != that.width)
					return false;
				return ABitMap::operator==(that);
			}
};

inline bool
operator!=(const ABitMap2 &bm1, const ABitMap2 &bm2)
{
	return !(bm1 == bm2);
}

inline ABitMap2
operator|(ABitMap2 bmLeft, const ABitMap2 &bmRight)
{
	return bmLeft |= bmRight;
}

inline ABitMap2
operator&(ABitMap2 bmLeft, const ABitMap2 &bmRight)
{
	return bmLeft &= bmRight;
}

inline ABitMap2
operator^(ABitMap2 bmLeft, const ABitMap2 &bmRight)
{
	return bmLeft ^= bmRight;
}

inline ABitMap2
operator-(ABitMap2 bmLeft, const ABitMap2 &bmRight)
{
	return bmLeft -= bmRight;
}

inline ABitMap2
operator~(ABitMap2 bmUnary)
{
	bmUnary.Invert();
	return bmUnary;
}

// Function to write a bitmap to a BMP file
extern bool	WriteBitMap(const ABitMap &bm, const char *fname);

// Function to write a 2-D bitmap to a BMP file
extern bool	WriteBitMap2(const ABitMap2 &bm2, const char *fname);

// Function to read a bitmap from a BMP file
extern bool	ReadBitMap(ABitMap *bmp, const char *fname);

// Function to read a 2-D bitmap from a BMP file
extern bool	ReadBitMap2(ABitMap2 *bm2p, const char *fname);

#endif	// ! _ABITMAP_H_
