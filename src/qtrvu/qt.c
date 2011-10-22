#ifndef lint
static const char	RCSid[] = "$Id: qt.c,v 1.1 2011/10/22 22:38:10 greg Exp $";
#endif
/*
 *  qt.c - driver for qt
 */
#include "math.h"
#include "copyright.h"
#include  "standard.h"
#include  "platform.h"
#include  "color.h"
#include  "driver.h"
#include  "rpaint.h"
#include "color.h"

static char lastPrompt[1024];
static const char* currentCommand =0;
static int abort_render = 0;
static int progress = 0;
static int last_total_progress = 0;
#define GAMMA		2.2		/* default exponent correction */

static dr_closef_t qt_close;
static dr_clearf_t qt_clear;
static dr_paintrf_t qt_paintr;
static dr_getcurf_t qt_getcur;
static dr_comoutf_t qt_comout;
static dr_cominf_t qt_comin;
static dr_flushf_t qt_flush;

static struct driver  qt_driver = {
  qt_close, qt_clear, qt_paintr, qt_getcur,
  qt_comout, qt_comin, qt_flush, 1.0
};

/* functions from qt_rvu_main.cxx */
extern qt_rvu_init(char*, char*, int*, int* );
extern qt_resize(int, int);
extern qt_rvu_paint(int r,int g,int b,int xmin,
                    int ymin,int xmax,int ymax);
extern void qt_flush_display();
extern int qt_rvu_run();
extern void qt_window_comout(const char*);
extern void qt_set_progress(int);
extern int qt_open_text_dialog(char* , const char*);

extern struct driver *
qt_init(		/* initialize driver */
  char  *name,
  char  *id
)
{
  lastPrompt[0] = 0;
  make_gmap(GAMMA);
  qt_rvu_init(name, id, &qt_driver.xsiz, &qt_driver.ysiz);
  return(&qt_driver);
}


static void
qt_close(void)			/* close our display */
{
  fprintf(stderr, "Qt close\n");
}

static void
qt_clear(			/* clear our display */
  int  xres,
  int  yres
)
{
  qt_resize(xres, yres);
}


static void
qt_paintr(		/* fill a rectangle */
  COLOR  col,
  int  xmin,
  int  ymin,
  int  xmax,
  int  ymax
)
{
  uby8 rgb[3];
  map_color(rgb, col);
  qt_rvu_paint(rgb[RED], rgb[GRN], rgb[BLU],
               xmin, ymin, xmax, ymax);
}

#define SCALE 1.2

static void
qt_flush(void)			/* flush output */
{
  int p;
  if(last_total_progress)
    {
    progress++;
    p = pdepth*10 +
      (int)floor(10 *
                 atan( SCALE * progress / last_total_progress) * 2.0/3.141593);
    qt_set_progress(p);
    }
  qt_flush_display();
}


static int comincalls = 0;
static void
qt_comin(		/* read in command line qt gui */
  char  *inp,
  char  *prompt
)
{
  if (prompt && strlen(prompt))
    {
    qt_comout(prompt);
    }
  if(comincalls < 2)
    {
    comincalls++;
    }
  /* if prompt is not set, then use the lastPrompt which
     comes from the last call to comout */
  if(!prompt)
    {
    prompt = lastPrompt;
    }
  /* first time call to qt_comin start qt event loop */
  if(comincalls == 1)
    {
    qt_rvu_run();
    /* when qt is done we are done so exit */
    exit(0);
    }
  /* This code will be called from inside the running
     qt_rvu_run loop in two cases:
     1. when a user enters a new command to the text widget. In
     this case the currentCommand is not null.
     2. when a command calls comin directly, that is the else
     case. */
  if(currentCommand)
    {
    /* A user has entered a new command */
    qt_comout((char*)currentCommand);
    /* copy the user command into inp, set
       currentCommand to null so as not to use it again,
       and then return */
    strcpy(inp, currentCommand);
    currentCommand = 0;
    return;
    }
  else
    {
    /* in this case comin has been called from inside a
       command call like view.
       OK tricky bit here.  If the prompt is empty, it
       is just trying to get the user to hit enter so that
       the most recent comout does not scrool away. Since the
       qt interface has a scroll bar we don't have to do that,
       and do not want to prompt the user. Also, if the prompt
       is set to "done: ", then we don't want to prompt the user.
       In other cases we assume that we are inside a command
       that is asking for more input and we use a dialog box
       to get it. */
    if(strlen(prompt) > 0 && strcmp(prompt, "done: ") != 0)
      {
      /* if the user entered new text, then echo that
         to the output window. */
      if(qt_open_text_dialog(inp, prompt))
        {
        qt_window_comout(inp);
        }
      lastPrompt[0] = 0;
      }
    }
}

static void
qt_comout(		/* write out string to stdout */
  char	*outp
)
{
  if (!outp)
    {
    return;
    }
  if (strlen(outp) > 0)
    {
    if(outp[strlen(outp)-2] == ':')
      {
      strcpy(lastPrompt, outp);
      }
    qt_window_comout(outp);
    }
}

extern void qt_rvu_get_position(int *x, int *y);

static int
qt_getcur(		/* get cursor position */
  int  *xp,
  int  *yp
)
{
  qt_rvu_get_position(xp, yp);
  return 0;
}

void qt_set_abort(int value)
{
  dev->inpready = value;
}


/* process one command from the GUI */
void qt_process_command(const char* com)
{
  char  buf[512];
  buf[0] = 0;
  /* set the currentCommand to the command */
  currentCommand = com;
  /* Call command with no prompt, this will in
     turn call qt_comin which will use the currentCommand
     pointer as the command to process */
  command("");
  /* after processing a command try to do a render */
  for ( ; ; )
    {
    if (hresolu <= 1<<pdepth && vresolu <= 1<<pdepth)
      {
      qt_set_progress(100);
      qt_comout("done");
      return;
      }
    errno = 0;
    if (hresolu <= psample<<pdepth && vresolu <= psample<<pdepth)
      {
      sprintf(buf, "%d sampling...\n", 1<<pdepth);
      qt_comout(buf);
      last_total_progress = progress ? progress:1;
      progress = 0;
      rsample();
      }
    else
      {
      sprintf(buf, "%d refining...\n", 1<<pdepth);
      qt_comout(buf);
      last_total_progress = progress ? progress:1;
      progress = 0;
      refine(&ptrunk, pdepth+1);
      }
    if (dev->inpready)
      {
      qt_comout("abort");
      dev->inpready = 0;
      pdepth = 10;
      return;
      }
    /* finished this depth */
    pdepth++;
    }
}

