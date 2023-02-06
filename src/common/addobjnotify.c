#ifndef lint
static const char	RCSid[] = "$Id: addobjnotify.c,v 3.4 2023/02/06 22:40:21 greg Exp $";
#endif
/*
 * Dummy declaration of addobjnotify[]
 */

#include "copyright.h"

#include <stdio.h>
#include "tiff.h"
#include "fvect.h"
#include "object.h"

void  (*addobjnotify[1])(OBJECT) = {0};
