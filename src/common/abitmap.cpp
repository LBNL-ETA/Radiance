#ifndef lint
static const char RCSid[] = "$Id: abitmap.cpp,v 2.2 2024/08/24 23:25:24 greg Exp $";
#endif
/*
 *  abitmap.cpp
 *  panlib
 *
 *  General bitmap class implementation
 *
 *  Created by gward on Wed Oct 31 2001.
 *  Copyright (c) 2022 Anyhere Software. All rights reserved.
 *
 */

#include <string.h>
#include <math.h>
#include "abitmap.h"

// Private workhorse bit copy function; handles overlaps but no length checks
static void
moveBits(uint32 *bdst, uint32 idst, const uint32 *bsrc, uint32 isrc, uint32 n)
{
	if (!n)
		return;
	bdst += idst >> 5; idst &= 0x1f;
	bsrc += isrc >> 5; isrc &= 0x1f;
	if ((bdst == bsrc) & (idst == isrc))
		return;
	uint32	sword[2];
	if (n <= 32) {				// short string case
		sword[0] = bsrc[0];
		sword[1] = bsrc[isrc+n > 32];
		bsrc = sword;
		while (n--) {
			*bdst &= ~(1 << idst);
			*bdst |= ((*bsrc >> isrc) & 1) << idst;
			if (++idst > 0x1f) { idst=0; ++bdst; }
			if (++isrc > 0x1f) { isrc=0; ++bsrc; }
		}
		return;
	}
	const bool	reverse = (bdst==bsrc) ? (idst > isrc) : (bdst > bsrc);
	const int	boff = isrc-idst & 0x1f;
	const bool	woff = (isrc > idst);
	const bool	partstart = (idst != 0);
	const bool	partend = ((idst+n & 0x1f) != 0);
	const int	lastword = (idst+n) >> 5;
	uint32		mask;
				// messy starting-word stuff
#define DO_FIRSTPART if (partstart) { \
		sword[0] = bsrc[0]; sword[1] = bsrc[1]; \
		bdst[0] &= mask = (1<<idst)-1; \
		if (boff) { \
			bdst[0] |= sword[woff]<<(32-boff) & ~mask; \
			if (woff) bdst[0] |= sword[0]>>boff & ~mask; \
		} else \
			bdst[0] |= sword[0] & ~mask; \
	} else
				// messy ending-word stuff
#define DO_LASTPART if (partend) { \
		mask = ~0 << (idst+n & 0x1f); \
		bool	beyond = (1<<(32-boff) & ~mask); \
		sword[0] = bsrc[lastword-1]; \
		if (!boff | woff | beyond) sword[1] = bsrc[lastword]; \
		bdst[lastword] &= mask; \
		if (boff) { \
			bdst[lastword] |= (sword[woff]>>boff) & ~mask; \
			if (beyond) \
				bdst[lastword] |= (woff ? bsrc[lastword+1] : sword[1]) \
							<< (32-boff) & ~mask; \
		} else \
			bdst[lastword] |= sword[1] & ~mask; \
	} else
	if (reverse)
		DO_LASTPART;
	else
		DO_FIRSTPART;
	if (boff) {				// middle part (unaligned case)
		int	i;
		if (reverse) {
			for (i = lastword; --i > 0; )
				bdst[i] = bsrc[i+woff]<<(32-boff) | bsrc[i+woff-1]>>boff;
			if (!partstart) {
				bdst[0] = bsrc[woff]<<(32-boff);
				if (woff) bdst[0] |= bsrc[0]>>boff;
			}
		} else {
			if (!partstart) {
				bdst[0] = bsrc[woff]<<(32-boff);
				if (woff) bdst[0] |= bsrc[0]>>boff;
			}
			for (i = 0; ++i < lastword; )
				bdst[i] = bsrc[i+woff]<<(32-boff) | bsrc[i+woff-1]>>boff;
		}
	} else {				// middle (aligned word case)
		memmove(bdst+partstart, bsrc+partstart,
				(lastword-partstart)*sizeof(uint32));
	}
	if (reverse)
		DO_FIRSTPART;
	else
		DO_LASTPART;
#undef DO_FIRSTPART
#undef DO_LASTPART
}

// Create and clear a new bitmap
bool
ABitMap::NewBitMap(uint32 n, bool clrset)
{
	int32	onwords = bmlen();
	len = n;
	int32	nwords = bmlen();
	if (nwords != onwords) {
		delete [] bmap;
		if (nwords) bmap = new uint32 [nwords];
		else bmap = 0;
	}
	if (!nwords)
		return false;
	ClearBitMap(clrset);
	return true;
}

// Clear bitmap to given value
void
ABitMap::ClearBitMap(bool clrset)
{
	memset(bmap, clrset * 0xff, sizeof(uint32)*bmlen());
}

// Invert the entire bitmap
void
ABitMap::Invert()
{
	uint32 *	wp = bmap + bmlen();
	while (wp-- > bmap)
		*wp = ~*wp;
}

// Extract bitmap section
bool
ABitMap::GetBits(ABitMap *dp, uint32 i) const
{
	if (!dp | (dp == this))
		return false;
	if (i >= len)
		return false;
	if (!dp->len && !dp->NewBitMap(len-i))
		return false;
	if (!i & (dp->len == len)) {
		*dp = *this;
		return true;
	}
	uint32	n = dp->len;
	if (n > len-i)
		n = len-i;
	moveBits(dp->bmap, 0, bmap, i, n);
	return true;
}

// Overlay bitmap section (ignores bits past end)
bool
ABitMap::AssignBits(uint32 i, const ABitMap &src)
{
	if (!src.len)
		return true;
	if (!i & (src.len == len)) {
		*this = src;
		return true;
	}
	if (i >= len)
		return false;
	moveBits(bmap, i, src.bmap, 0, (len-i < src.len) ? len-i : src.len);
	return true;
}

// Apply operation to bitmap section
bool
ABitMap::OpBits(uint32 i, char op, const ABitMap &src)
{
	if (!src.len | !len)
		return false;
	if (op == '=')
		return AssignBits(i, src);
	ABitMap	bits(src.len);
	if (!GetBits(&bits, i))
		return false;
	switch (op) {
	case '|':
		bits |= src;
		break;
	case '&':
		bits &= src;
		break;
	case '^':
		bits ^= src;
		break;
	case '-':
	case '>':
		bits -= src;
		break;
	case '<':
		bits.Invert();
		bits &= src;
		break;
	default:
		return false;
	}
	return AssignBits(i, bits);
}

// Clear bitmap section
void
ABitMap::ClearBits(uint32 i, uint32 n, bool clrset)
{
	if (i >= len)
		return;
	if (n >= len - i) {
		if (!i) {
			ClearBitMap(clrset);
			return;
		}
		n = len - i;
	} else if (!n)
		return;

	const uint32 * const	sectEnd = bmap + ((i+n)>>5);
	uint32 *		wp = bmap + (i>>5);
	if (wp == sectEnd) {		// single word clear?
		const uint32	bits = (~0 << (i & 0x1f) &
					(1 << (i+n & 0x1f)) - 1);
		if (clrset)
			*wp |= bits;
		else
			*wp &= ~bits;
		return;
	}
	const uint32		clrWord = clrset * ~0;
	if (i & 0x1f) {			// partial first word?
		if (clrset)
			*wp++ |= ~0 << (i & 0x1f);
		else
			*wp++ &= (1 << (i & 0x1f)) - 1;
	}
	while (wp < sectEnd)		// central words
		*wp++ = clrWord;
	if (i+n & 0x1f) {		// partial last word?
		if (clrset)
			*wp |= (1 << (i+n & 0x1f)) - 1;
		else
			*wp &= ~0 << (i+n & 0x1f);
	}
}

// Total all bits set in bitmap
uint32
ABitMap::SumTotal(bool bit2cnt) const
{
	static char	bitCount[256];
	int		i;
	if (!bitCount[255]) {
		for (i = 256; --i; )
			for (int bits = i; bits; bits >>= 1)
				bitCount[i] += (bits & 1);
	}
	uint32			count = 0;
	const unsigned char *	cp;
	cp = (const unsigned char *)(bmap + bmlen());
	if (len & 0x1f) {		// partial last word
		uint32		lastBits = WordV(len-1);
		lastBits &= (1<<(len&0x1f)) - 1;
		for (i = sizeof(uint32); i--; lastBits >>= 8)
			count += bitCount[lastBits & 0xff];
		cp -= sizeof(uint32);
	}
	while (cp > (const unsigned char *)bmap)
		count += bitCount[*--cp];
	if (bit2cnt)
		return count;
	return len - count;
}

// Return the next bit position matching val (or ABMend if no match)
uint32
ABitMap::Find(uint32 i, bool val) const
{
	const uint32	clrWord = !val * ~0;
	const uint32 *	wp = bmap + (i>>5);
	uint32		b = Bit(i);

	wp -= !(i & 0x1f);

	while (i < len) {
		if (!(i & 0x1f)) {
			while (*++wp == clrWord)
				if ((i += 0x20) >= len)
					return ABMend;
			b = 1;
		}
		if (((*wp & b) != 0) == val)
			return i;
		b <<= 1;
		++i;
	}
	return ABMend;
}

// Shift bits downward in bitmap, zero fill
ABitMap &
ABitMap::operator>>=(uint32 nbits)
{
	if (!nbits)
		return *this;
	if (nbits >= len) {
		ClearBitMap();
		return *this;
	}
	moveBits(bmap, 0, bmap, nbits, len-nbits);
	ClearBits(len-nbits, nbits);
	return *this;
}

// Shift bits upwards in bitmap, zero fill
ABitMap &
ABitMap::operator<<=(uint32 nbits)
{
	if (!nbits)
		return *this;
	if (nbits >= len) {
		ClearBitMap();
		return *this;
	}
	moveBits(bmap, nbits, bmap, 0, len-nbits);
	ClearBits(0, nbits);
	return *this;
}

// Bitmap copy operator
ABitMap &
ABitMap::operator=(const ABitMap &src)
{
	if (this == &src)
		return *this;
	int32	nwords = src.bmlen();
	if (nwords != bmlen()) {
		delete [] bmap;
		if (nwords) bmap = new uint32 [nwords];
		else bmap = 0;
	}
	len = src.len;
	memcpy(bmap, src.bmap, nwords*sizeof(uint32));
	return *this;
}

// Bitmap OR-copy operator (reverts to copy for different size bitmaps)
ABitMap &
ABitMap::operator|=(const ABitMap &src)
{
	if (this == &src)
		return *this;
	if (len != src.len)
		return *this = src;
	const int32	nwords = bmlen();
	const uint32 *	owp = src.bmap + nwords;
	uint32 *	wp = bmap + nwords;
	while (wp > bmap)
		*--wp |= *--owp;
	return *this;
}

// Bitmap AND-assign operator (no effect for different size bitmaps)
ABitMap &
ABitMap::operator&=(const ABitMap &src)
{
	if (this == &src)
		return *this;
	if (len != src.len)
		return *this;
	const int32	nwords = bmlen();
	const uint32 *	owp = src.bmap + nwords;
	uint32 *	wp = bmap + nwords;
	while (wp > bmap)
		*--wp &= *--owp;
	return *this;
}

// Bitmap XOR-assign operator (no effect for different size bitmaps)
ABitMap &
ABitMap::operator^=(const ABitMap &src)
{
	if (this == &src) {
		ClearBitMap();
		return *this;
	}
	if (len != src.len)
		return *this;
	const int32	nwords = bmlen();
	const uint32 *	owp = src.bmap + nwords;
	uint32 *	wp = bmap + nwords;
	while (wp > bmap)
		*--wp ^= *--owp;
	return *this;
}

// Clear bits set in second bitmap
bool
ABitMap::ClearBitsFrom(const ABitMap &src)
{
	if (this == &src) {
		ClearBitMap();
		return true;
	}
	if (src.len != len)
		return false;
	uint32 *	wp = bmap + bmlen();
	const uint32 *	sp = src.bmap + src.bmlen();
	while (wp > bmap)
		*--wp &= ~*--sp;
	return true;
}

// Compare two bitmaps for equality
bool
ABitMap::operator==(const ABitMap &that) const
{
	if (this == &that)
		return true;
	if (len != that.len)
		return false;
	if (!len)
		return true;
	const int	nwords = len >> 5;
	const uint32 *	owp = that.bmap + nwords;
	const uint32 *	wp = bmap + nwords;
	if (len & 0x1f && (*--wp ^ *--owp) & (1<<(len&0x1f))-1)
		return false;
	while (wp > bmap)
		if (*--wp != *--owp)
			return false;
	return true;
}

/*************** Run-length compander section ***************
 *  First bit is:
 * 	0	=> non-run
 *	1	=> run
 *  Next N "0" bits indicate counter length-3, "1"-terminated:
 *	001	=> e.g., 5-bit counter follows
 *  Next N+3 bits constitute counter M, with implied "1" in leftmost bit
 *	00001	=> e.g., 33 run or non-run length
 *  Next bit is "0" for a run of M 0's, or "1" for a run of 1's,
 *  OR M bits of non-run data.
 *
 *  Example 1:	"00101101"	=> "0 1 000 00101101"
 *  Example 2:	49 1's		=> "1 001 10001 1"
 *  Example 3:	7326 0's	=> "1 0000000001 110010011110 0"
 *
 *  Note that any run or non-run must span 8 bits or more in source, and
 *  will take up at least 6 encoded bits for a run and 13 for a non-run.
 *  Encoding a bitmap < 72 bits long is never a win, and < 8 is
 *  not possible.  A bitmap <= 64 bits in length will be passed as is.
 *  Decoding is trivial compared to logical constraints during encode.
 */

#define RLEmagic	0x1700A9A5	// magic number (32-bits)
#define	MinCntB		3		// minimum counter length
#define MinCnt		(1<<MinCntB)	// minimum encodable bit sequence

// Get original bitmap length if RLE (or 0)
uint32
ABitMap::RLength() const
{
	if (bmlen() < 3 || bmap[0] != RLEmagic) return 0;
	return bmap[1];
}

// Compress into a run-length encoded bitmap
bool
ABitMap::GetRLE(ABitMap *rlep) const
{
	if (!rlep)
		return false;
	if (RLength())			// already encoded?
		return false;
	if (len <= 64) {		// don't bother?
		*rlep = *this;
		return len;
	}
					// create draft bitmap
	ABitMap	tmap(len + len/50 + 128);
	tmap.bmap[0] = RLEmagic;	// mark as RLE
	tmap.bmap[1] = len;		// record original length
	uint32	i=0, o=64;		// encode bits
	while (i < len) {
		uint32	cnt, cbex;
		uint32	start = i++;	// check for usable run
		if (!Find(&i, !Check(start)))
			i = len;
		else if (i > len - MinCnt)
			i = len - MinCnt;

		if (i >= start + MinCnt) {
			tmap.Set(o++);	// encode a run
			cnt = (i - start)>>MinCntB;
			cbex = MinCntB;	// counter width
			while (cnt > 1) {
				o++;
				cbex++;
				cnt >>= 1;
			}
			tmap.Set(o++);	// 1 terminator, then count
			cnt = i - start;
			while (cbex-- > 0) {
				if (cnt & 1<<cbex) tmap.Set(o);
				++o;
			}		// followed by repeat bit
			if (Check(start))
				tmap.Set(o);
			++o;
			continue;	// move on to next
		}
		i = start + MinCnt;	// else encode non-run
		if (i + MinCnt > len)
			i = len;	// non-run to end
		while (i < len) {	// stop at next useful run
			bool	cand = Check(i);
			uint32	candStart = i;
			while (++i < len && Check(i) == cand)
				if (i - candStart > 2*(2+MinCntB)+1) {
					i = candStart;
					candStart = ABMend;
					break;
				}
			if (candStart == ABMend)
				break;	// found run worth stopping for
		}
		o++;			// encode our non-run
		cnt = (i - start)>>MinCntB;
		if (!cnt)
			goto calamity;	// should never occur!
		cbex = MinCntB;		// counter width
		while (cnt > 1) {
			o++;
			cbex++;
			cnt >>= 1;
		}
		tmap.Set(o++);		// 1 terminator, then count
		cnt = i - start;
		while (cbex-- > 0) {
			if (cnt & 1<<cbex) tmap.Set(o);
			++o;
		}			// finally, copy bit data
		if (o + cnt > tmap.len)
			goto calamity;	// over-ran temp buffer!
		moveBits(tmap.bmap, o, bmap, start, cnt);
		o += cnt;		// onwards...
	}
					// copy to right-sized array
	if (rlep->NewBitMap(o) && tmap.GetBits(rlep, 0))
		return true;
calamity:		// XXX should really speak up if we get here!
	return false;
}

// Reconstitute bits from RLE encoding
bool
ABitMap::SetFromRLE(const ABitMap &rle)
{
	if (rle.len <= 64) {		// never was encoded?
		*this = rle;
		return len;
	}
	if (&rle == this) {		// cuidado!
		ABitMap		tmap;
		return tmap.Take(this) && SetFromRLE(tmap);
	}
	if (!NewBitMap(rle.RLength()))	// start from 0's
		return false;
	uint32	i=64, o=0;		// decode bits
	while (i < rle.len) {
		bool	isrun = rle.Check(i++);
		int	cntlen = MinCntB;
		while (!rle.Check(i++)) cntlen++;
		if (cntlen > 31)
			return false;
		uint32	rlen = 1;	// get output count
		while (cntlen--)
			rlen = (rlen<<1) + rle.Check(i++);
		if (!isrun) {		// copy bits
			if ((i+rlen > rle.len) | (o+rlen > len))
				return false;
			moveBits(bmap, o, rle.bmap, i, rlen);
			i += rlen;
		} else if (rle.Check(i++))
			ClearBits(o, rlen, true);

		if ((o += rlen) >= len)	// advance output index
			break;
	}
	return (i == rle.len) & (o == len);
}

/*************** 2-D bitmap section ***************/

// Reconstitute bits from RLE encoding (size must match)
bool
ABitMap2::SetFromRLE(int w, int h, const ABitMap &rle)
{
	uint32	orig_len = rle.RLength();

	if ((w <= 0) | (h <= 0) | (orig_len != (uint32)w*h))
		return false;

	if (!ABitMap::SetFromRLE(rle))
		return false;

	width = w; height = h;
	return true;
}

// Extract bitmap section, true if some overlap
bool
ABitMap2::GetRect(ABitMap2 *dp, int sx, int sy) const
{
	if (!dp | (dp == this))
		return false;
	if ((sx >= width) | (sy >= height))
		return false;
	if (dp->width <= 0 && !dp->NewBitMap(width-sx, height-sy))
		return false;
	if (!sx & !sy & (dp->width == width) & (dp->height == height)) {
		*dp = *this;
		return true;
	}
	int		dx=0, dy=0;
	if (sx < 0) {
		if ((dx = -sx) >= dp->width) return false;
		sx = 0;
	}
	if (sy < 0) {
		if ((dy = -sy) >= dp->height) return false;
		sy = 0;
	}
	const int	rowwidth = (dp->width-dx > width-sx) ? width-sx : dp->width-dx;
	int		rowcount = (dp->height-dy > height-sy) ? height-sy : dp->height-dy;
	if ((rowwidth == width) & (width == dp->width))
		moveBits(dp->base(), dp->width*dy, base(), width*sy, rowwidth*rowcount);
	else
		while (rowcount--)
			moveBits(dp->base(), dp->width*dy++ + dx,
					base(), width*sy++ + sx, rowwidth);
	return true;	
}

// Assign bitmap section (ignores anything past edges)
bool
ABitMap2::AssignRect(int dx, int dy, const ABitMap2 &src)
{
	if (src.width <= 0)
		return true;
	if ((dx >= width) | (dy >= height))
		return false;
	if (!dx & !dy && (src.width == width) & (src.height == height)) {
		*this = src;
		return true;
	}
	int		sx=0, sy=0;
	int		w=src.width, h=src.height;
	if (dx < 0) {
		if ((sx = -dx) >= w) return false;
		dx = 0; w -= sx;
	}
	if (dy < 0) {
		if ((sy = -dy) >= h) return false;
		dy = 0; h -= sy;
	}
	if (dx+w > width) w = width - dx;
	if (dy+h > height) h = height - dy;
	if ((w <= 0) | (h <= 0))
		return false;
	if ((w == width) & (width == src.width))
		moveBits(base(), width*dy, src.base(), src.width*sy, w*h);
	else
		while (h--)
			moveBits(base(), width*dy++ + dx,
					src.base(), src.width*sy++ + sx, w);
	return true;
}

// Apply operation to bitmap section
bool
ABitMap2::OpRect(int dx, int dy, char op, const ABitMap2 &src)
{
	if ((src.width <= 0) | (width <= 0))
		return false;
	if (op == '=')
		return AssignRect(dx, dy, src);
	ABitMap2	rbits(src.width, src.height);
	if (!GetRect(&rbits, dx, dy))
		return false;
	switch (op) {
	case '|':
		rbits |= src;
		break;
	case '&':
		rbits &= src;
		break;
	case '^':
		rbits ^= src;
		break;
	case '-':
	case '>':
		rbits -= src;
		break;
	case '<':
		rbits.Invert();
		rbits &= src;
		break;
	default:
		return false;
	}
	return AssignRect(dx, dy, rbits);
}

// Clear a rectangle
void
ABitMap2::ClearRect(int x, int y, int w, int h, bool clrset)
{
	if (x < 0) { w += x; x = 0; }
	if (y < 0) { h += y; y = 0; }
	if (w > width - x)
		w = width - x;
	if (w <= 0)
		return;
	if (h > height - y)
		h = height - y;
	if (h <= 0)
		return;
	if (w == width)			// contiguous case
		ClearBits(y*width, h*width, clrset);
	else
		while (h--)		// discontiguous
			ClearBits(bmi(x,y++), w, clrset);
}

// Get bounds of assigned region
bool
ABitMap2::GetBoundRect(int xymin[2], int wh[2], bool val) const
{
	int	x=0, y=0;

	if (!xymin | !wh)
		return false;
	if (!Find(&x, &y, val)) {
		xymin[0] = xymin[1] = wh[0] = wh[1] = 0;
		return false;
	}
	xymin[0] = x; xymin[1] = y;
	int	xymax[2] = {x, y};
	do {
		if (x < xymin[0]) xymin[0] = x;
		xymax[1] = y;

		if (Check(width-1, y) == val) {
			xymax[0] = width-1;
			x = 0; ++y;
			continue;
		}
		if (++x < xymax[0]) x = xymax[0];
		Find(&x, &y, !val);	// never false after above Check()
		if (x-1 > xymax[0]) xymax[0] = x-1;
	} while (Find(&x, &y, val));

	wh[0] = xymax[0] - xymin[0] + 1;
	wh[1] = xymax[1] - xymin[1] + 1;
	return true;
}

// Shift bitmap image right-left and down-up, filling as indicated
void
ABitMap2::Shift(int dx, int dy, int fill)
{
	int	adx = (dx > 0) ? dx : -dx;
	int	ady = (dy > 0) ? dy : -dy;
	if ((adx >= width) | (ady >= height)) {
		if (fill >= 0)
			ClearBitMap(fill);
		return;
	}
	int	bitshift = dy*width + dx;
	if (bitshift > 0)		// shift underlying bits
		operator<<=(bitshift);
	else if (bitshift < 0)
		operator>>=(-bitshift);
	else
		return;
	if (fill < 0)			// no fill -- we're done
		return;
	int	hmarg[4];		// new horizontal margin
	hmarg[0] = (dx < 0) ? width-adx : 0;
	hmarg[1] = (dy > 0) ? dy : 0;
	hmarg[2] = adx;
	hmarg[3] = height - ady;
	if (fill) {			// fill corner with 1's
		ClearRect(hmarg[0], hmarg[1], hmarg[2], hmarg[3], true);
		ClearRect(0, (dy < 0) ? height-ady : 0, width, ady, true);
	} else {			// fill side with 0's
		ClearRect(hmarg[0], hmarg[1]-1, hmarg[2], hmarg[3]+2);
	}
}

static inline int
iSqrt(double x)
{
	if (x <= 0) return 0;
	return int(sqrt(x) + .5);
}

// Dilate (or erode) selection by given radius
void
ABitMap2::Expand(double rad, bool val)
{
	if (rad < 0) {
		rad = -rad;
		val = !val;
	}
	if ((width <= 0) | (rad < 1))
		return;
					// check what we have here
	int		xyorg[2], wh[2];
	if (!GetBoundRect(xyorg, wh, val))
		return;			// empty bitmap!
	const int	wh_orig[2] = {width, height};
	const int	irad = int(rad+.5);
					// optimize if >= 3/4 empty
	if ((width >= 128) & (height >= 128) &&
				((wh[0] += irad<<1) <= width>>1) &
				((wh[1] += irad<<1) <= height>>1)) {
		ABitMap2	subRgn(wh[0], wh[1], !val);
		GetRect(&subRgn, xyorg[0] -= irad, xyorg[1] -= irad);
		Take(&subRgn);		// work with subregion
	}
					// copy original pattern
	ABitMap2	origMap = *this;
					// stamp out larger one
	for (int y = -irad; y <= irad; y++) {
		const int	spanRad = iSqrt(rad*rad - y*y);
		for (int x = -spanRad; x <= spanRad; x++) {
			if (!x & !y) continue;
			ABitMap2	stamp = origMap;
			stamp.Shift(x, y, !val);
			if (val)
				*this |= stamp;
			else
				*this &= stamp;
		}
	}
	if ((width == wh_orig[0]) & (height == wh_orig[1]))
		return;
					// restore original dimensions
	origMap.NewBitMap(wh_orig[0], wh_orig[1], !val) &&
		origMap.AssignRect(xyorg[0], xyorg[1], *this) &&
		Take(&origMap);
}
