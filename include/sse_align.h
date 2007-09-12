/*******************************************************************************
 * $Id: sse_align.h,v 1.1 2007-09-12 19:33:13 bjoo Exp $
 * 
 *
 * Define alignment macro ALIGN. This can vary depending on compiler 
 *
 * Author: Balint Joo 
 * Date: 11/13/2002
 *
 *******************************************************************************/

/* Include guard... */
#ifndef __INCLUDED_SSE_ALIGN_H__
#define __INCLUDED_SSE_ALIGN_H__

#include <sse_config.h>

#ifndef ALIGN
#if ((defined SSE)||(defined SSE2))


#if defined P4
/* 64 bit alignment on P4 targets */
#define ALIGN_BITS (128)
#else 
/* 32 bit alignment on non P3 targets */
//#define ALIGN_BITS (32)
#define ALIGN_BITS (128)
#endif /* If defined P4 */

/* Now get the right alignment for your particular compiler */

/* Gnu compilers */
#if __GNUC__

/* GCC-2.96 -- buggy -- wants alignment in bits */
#if __GNUC__ == 2 && __GNUC_MINOR__ == 96
#warning "GCC 2.96 detected aligning in terms of bits"
#define ALIGN_ATTRIB  ALIGN_BITS
#else 
/* Other GNU targets -- want alignment in bytes ( divide by 8 == shiftright by 8 */
#define ALIGN_ATTRIB  (ALIGN_BITS >> 3)
#endif

#define ALIGN __attribute__ ((aligned (ALIGN_ATTRIB)))
#else

/* Other compiler targets */
#error "Non GNU Compilers not supported"

#endif /* ifdef __GNUC__ */
 
#else /* if defined SSE || defined SSE2 */

#define ALIGN 

#endif /* if defined SSE || defined SSE2 */

#endif /* ifndef ALIGN */

/* End of include guard */
#endif 


