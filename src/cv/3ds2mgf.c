#ifndef lint
static const char	RCSid[] = "$Id: 3ds2mgf.c,v 2.1 2011/02/18 00:40:25 greg Exp $";
#endif
/*
      3DS2POV.C  by Steve Anger and Jeff Bowermaster
			MGF output added by Greg Ward

      Reads a 3D Studio .3DS file and writes a POV-Ray, Vivid,
      Polyray, MGF or raw scene file.

      Version 2.0 Written Feb/96

      Compiled with MSDOS GNU C++ 2.4.1 or generic ANSI-C compiler
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "rayopt3ds.h"

#ifdef __TURBOC__
#include <alloc.h>
extern unsigned _stklen = 16384;
#endif


#define FALSE 0
#define TRUE  1

/* Internal bounding modes */
#define OFF  0
#define ON   1
#define AUTO 2

#define MAX_LIB  10
#define ASPECT   1.333

/* Output formats */
#define POV10    0
#define POV20    1
#define VIVID    2
#define POLYRAY  3
#define MGF      4
#define RAW      99

#define DEG(x) ((double)(180.0/M_PI)*(x))
#define RAD(x) ((double)(M_PI/180.0)*(x))

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#ifndef MAXFLOAT
#define MAXFLOAT (1e37)
#endif
				/* RGB chromaticity definitions for MGF */
#define CIE_x_r		0.640
#define CIE_y_r		0.330
#define CIE_x_g		0.290
#define CIE_y_g		0.600
#define CIE_x_b		0.150
#define CIE_y_b		0.060
				/* computed luminances from above */
#define CIE_Y_r		0.265
#define CIE_Y_g		0.670
#define CIE_Y_b		0.065

/* A generic list type */
#define LIST_INSERT(root, node) list_insert ((List **)&root, (List *)node)
#define LIST_FIND(root, name)   list_find   ((List **)&root, name)
#define LIST_DELETE(root, node) list_delete ((List **)&root, (List *)node)
#define LIST_KILL(root)         list_kill   ((List **)&root)

#define LIST_FIELDS  \
    char name[80];   \
    void *next;


typedef unsigned char  byte;
typedef unsigned short word;
typedef unsigned long  dword;

typedef struct {
    LIST_FIELDS
} List;


typedef struct {
    int a, b, c;
} Face;


typedef struct {
    float red, green, blue;
} Colour;


/* Transformation command */
typedef struct {
    LIST_FIELDS

    Matrix matrix;
} Transform;


/* Morph command */
typedef struct {
    LIST_FIELDS

    int    count;          /* Number of objects in morph */
    char   names[4][80];   /* Name of n'th object in average */
    float  weight[4];      /* Weight applied to n'th object */

    Matrix matrix;
} Morph;


/* Omni light command */
typedef struct {
    LIST_FIELDS

    Vector pos;            /* Light position */
    Colour col;            /* Light colour */
} OmniLight;


/* Spotlight command */
typedef struct {
    LIST_FIELDS

    Vector pos;            /* Spotlight position */
    Vector target;         /* Spotlight target location */
    Colour col;            /* Spotlight colour */
    float  hotspot;        /* Hotspot angle (degrees) */
    float  falloff;        /* Falloff angle (degrees) */
    int    shadow_flag;    /* Shadow flag (not used) */
} Spotlight;


/* Camera command */
typedef struct {
    LIST_FIELDS

    Vector pos;            /* Camera location */
    Vector target;         /* Camera target */
    float  bank;           /* Banking angle (degrees) */
    float  lens;           /* Camera lens size (mm) */
} Camera;


/* Material list */
typedef struct {
    LIST_FIELDS

    int  external;         /* Externally defined material? */
} Material;


/* Object summary */
typedef struct {
    LIST_FIELDS

    Vector center;         /* Min value of object extents */
    Vector lengths;        /* Max value of object extents */
} Summary;


/* Material property */
typedef struct {
    LIST_FIELDS

    Colour ambient;
    Colour diffuse;
    Colour specular;
    float  shininess;
    float  transparency;
    float  reflection;
    int    self_illum;
    int    two_side;
    char   tex_map[40];
    float  tex_strength;
    char   bump_map[40];
    float  bump_strength;
} MatProp;


/* Default material property */
MatProp DefaultMaterial = { "Default", NULL, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0},
			     {1.0, 1.0, 1.0}, 70.0, 0.0, 0.0, FALSE, FALSE };

/* A mesh object */
typedef struct {
    LIST_FIELDS

    int  vertices;         /* Number of vertices */
    Vector *vertex;        /* List of object vertices */

    int  faces;            /* Number of faces */
    Face *face;            /* List of object faces */
    Material **mtl;        /* Materials for each face */

    Matrix matrix;         /* Local mesh matrix */
    Matrix invmatrix;
    Vector center;         /* Center of object */
    Vector lengths;        /* Dimensions of object */

    int hidden;            /* Hidden flag */
    int shadow;            /* Shadow flag */
} Mesh;


typedef struct {
    dword start;
    dword end;
    dword length;
    word  tag;
} Chunk;


typedef struct {
    byte red;
    byte green;
    byte blue;
} Colour_24;


Colour Black = {0.0, 0.0, 0.0};

OmniLight *omni_list  = NULL;
Spotlight *spot_list  = NULL;
Camera    *cam_list   = NULL;
Mesh      *mesh_list  = NULL;
Transform *trans_list = NULL;
Morph     *morph_list = NULL;
Material  *mtl_list   = NULL;
List      *excl_list  = NULL;
List      *box_list   = NULL;
MatProp   *mprop_list = NULL;
Summary   *summary    = NULL;


FILE   *in;
FILE   *out;
char   inname[80];
char   outname[80];
char   vuename[80];
char   obj_name[80] = "";
Colour fog_colour = {0.0, 0.0, 0.0};
Colour col        = {0.0, 0.0, 0.0};
Colour global_amb = {0.1, 0.1, 0.1};
Vector pos        = {0.0, 0.0, 0.0};
Vector target     = {0.0, 0.0, 0.0};
float  fog_distance = 0.0;
float  hotspot = -1;
float  falloff = -1;
Mesh   *mesh = NULL;
int    frame = -1;
char   libname[MAX_LIB][80];
float  smooth = 60.0;
int    bound = 0;
int    verbose = 0;
int    format = POV20;
int    internal_bounding = AUTO;
int    box_all = FALSE;
int    cameras = 0;
int    libs = 0;
float  vue_version = 1.0;
Matrix *ani_matrix = NULL;
int    no_opt = FALSE;
FILE   *meshf = NULL;


void process_args (int argc, char *argv[]);
void parse_option (char *option);
void list_insert (List **root, List *new_node);
void *list_find (List **root, char *name);
void list_delete (List **root, List *node);
void list_kill (List **root);
Material *update_materials (char *new_material, int ext);
MatProp *create_mprop (void);
void read_library (char *fname);
void write_intro (FILE *f);
void write_summary (FILE *f);
void write_bgsolid (FILE *f, Colour col);
void write_light (FILE *f, char *name, Vector pos, Colour col);
void write_spot (FILE *f, char *name, Vector pos, Vector target, Colour col,
	 float hotspot, float falloff);
void write_fog (FILE *f, Colour col, float dist);
void write_camera (FILE *f, char *name, Vector pos, Vector target, float lens,
	 float bank);
void write_material (FILE *f, char *mat);
void write_pov10_material (FILE *f, MatProp *m);
void write_pov20_material (FILE *f, MatProp *m);
void write_vivid_material (FILE *f, MatProp *m);
void write_polyray_material (FILE *f, MatProp *m);
void write_mgf_material (FILE *f, MatProp *m);
void write_mesh (FILE *f, Mesh *mesh);
Transform *parse_transform (char *string);
Morph *parse_morph (char *string);
OmniLight *parse_omnilight (char *string);
Spotlight *parse_spotlight (char *string);
Camera *parse_camera (char *string);
void read_frame (char *filename, int frame_no);
void find_frame (FILE *f, int frame_no);
void save_animation (void);
Mesh *create_mesh (char *name, int vertices, int faces);
Mesh *copy_mesh (Mesh *mesh);
void free_mesh_data (Mesh *mesh);
void update_limits (Mesh *mesh);
char *before (char *str, char *target);
char *after (char *str, char *target);
char *between (char *str, char *target1, char *target2);
char *parse_string (char *str);
char upcase (char c);
float colour_intens (Colour *colour);
void parse_file (void);
void parse_3ds (Chunk *mainchunk);
void parse_mdata (Chunk *mainchunk);
void parse_fog (Chunk *mainchunk);
void parse_fog_bgnd (void);
void parse_mat_entry (Chunk *mainchunk);
char *parse_mapname (Chunk *mainchunk);
void parse_named_object (Chunk *mainchunk);
void parse_n_tri_object (Chunk *mainchunk);
void parse_point_array (void);
void parse_face_array (Chunk *mainchunk);
void parse_msh_mat_group (void);
void parse_smooth_group (void);
void parse_mesh_matrix (void);
void parse_n_direct_light (Chunk *mainchunk);
void parse_dl_spotlight (void);
void parse_n_camera (void);
void parse_colour (Colour *colour);
void parse_colour_f (Colour *colour);
void parse_colour_24 (Colour_24 *colour);
float parse_percentage (void);
short parse_int_percentage (void);
float parse_float_percentage (void);
void start_chunk (Chunk *chunk);
void end_chunk (Chunk *chunk);
byte read_byte (void);
word read_word (void);
dword read_dword (void);
float read_float (void);
void read_point (Vector v);
char *read_string (void);
float findfov (float lens);
int read_mgfmatname (char *s, int n, FILE *f);

char *progname;


int main (int argc, char *argv[])
{
    char meshfname[128];
    Material *m;
    int i;

    process_args (argc, argv);

    if (!no_opt) {
	opt_set_format (format);
	opt_set_dec (4);
	opt_set_bound (bound);
	opt_set_smooth (smooth);
	opt_set_quiet (!verbose);
	opt_set_fname (outname, "");
    } else if (format == MGF) {
	strcpy(meshfname, outname);
	add_ext(meshfname, "inc", 1);
	if (!strcmp(meshfname, outname)) {
	    printf ("Output and mesh file names are identical!\n");
	    exit (1);
	}
	if ((meshf = fopen (meshfname, "w")) == NULL) {
	    printf ("Cannot open mesh output file %s!\n", meshfname);
	    exit (1);
	}
    }

    if ((in = fopen (inname, "rb")) == NULL) {
	printf ("Cannot open input file %s!\n", inname);
	exit (1);
    }

    if ((out = fopen (outname, "w")) == NULL) {
	printf ("Cannot open output file %s!\n", outname);
	exit (1);
    }

    /* Load the names of pre-defined materials */
    for (i = 0; i < MAX_LIB; i++) {
	if (strlen(libname[i]) > 0)
	    read_library (libname[i]);
    }

    /* Load the instructions for the current frame */
    if (strlen(vuename) > 0)
	read_frame (vuename, frame);

    printf("Output to: %s\n", outname);

    if (frame >= 0)
	printf ("Generating frame #%d\n", frame);

    printf("\nPlease wait; Processing...\n");

    write_intro(out);

    parse_file();

    fclose(in);

    for (m = mtl_list; m != NULL; m = m->next) {
	if (!m->external)
	    write_material (out, m->name);
    }

    if (frame >= 0)
	save_animation();

    if (!no_opt) {
        write_summary (out);
        fflush (out);

	opt_finish();
    } else if (meshf != NULL) {
	fclose(meshf);
	fprintf (out, "i %s\n", meshfname);
    }

    fclose (out);

    LIST_KILL (omni_list);
    LIST_KILL (spot_list);
    LIST_KILL (cam_list);
    LIST_KILL (mesh_list);
    LIST_KILL (trans_list);
    LIST_KILL (morph_list);
    LIST_KILL (mtl_list);
    LIST_KILL (excl_list);
    LIST_KILL (box_list);
    LIST_KILL (mprop_list);
    LIST_KILL (summary);

    return 0;
}


/* Handle the command line args */
void process_args (int argc, char *argv[])
{
    int i;
    char *env_opt, *option;

    printf("\n\nAutodesk 3D Studio to Raytracer file Translator. Feb/96\n");
    printf("Version 2.0 by Steve Anger and Jeff Bowermaster 1996\n");
    printf ("\n");

    if (argc < 2) {
	printf ("Usage: %s inputfile[.3ds] [outputfile] [options]\n\n", argv[0]);
	printf ("Options: -snnn        - Smooth triangles with angles < nnn\n");
	printf ("         -l<filename> - Specifies native material library\n");
	printf ("         -a<filename> - Use animation information in specified file\n");
	printf ("         -fnnn        - Generate frame nnn of animation\n");
	printf ("         -x<object>   - Exclude this object from conversion\n");
	printf ("         -b<object>   - Convert this object as a box\n");
	printf ("         +i, -i       - Turn internal bounding on or off\n");
	printf ("         +v, -v       - Turn verbose status messages on or off\n");
	printf ("         -op          - Output to POV-Ray 2.0 format\n");
	printf ("         -op1         - Output to POV-Ray 1.0 format\n");
	printf ("         -ov          - Output to Vivid format\n");
	printf ("         -ol          - Output to poLyray format\n");
	printf ("         -om          - Output to MGF\n");
	printf ("         -or          - Output to RAW triangle format\n\n");
	printf ("ex. %s birdshow +v -lmaterials.inc\n\n", argv[0]);
	exit(1);
    }
					/* figure default format from name */
    progname = strrchr(argv[0], '/');
    if (progname == NULL) progname = argv[0];
    else progname++;
    if (!strcmp(progname, "3ds2pov"))
	format = POV20;
    else if (!strcmp(progname, "3ds2viv"))
	format = VIVID;
    else if (!strcmp(progname, "3ds2pi"))
	format = POLYRAY;
    else if (!strcmp(progname, "3ds2mgf"))
	format = MGF;
    else if (!strcmp(progname, "3ds2raw"))
	format = RAW;
    else
	format = POV20;		/* default if program name strange */

    strcpy (inname, "");
    strcpy (outname, "");
    strcpy (vuename, "");

    for (i = 0; i < MAX_LIB; i++)
	strcpy (libname[i], "");

    frame = -1;
    smooth = 70.0;
    bound = 0;
    verbose = 0;
    internal_bounding = AUTO;
    box_all = FALSE;
    libs = 0;

    /* Parse the enviroment string options */
    env_opt = getenv ("3DS2POV");

    if (env_opt != NULL) {
	option = parse_string (env_opt);

	while (strlen(option) > 0) {
	    parse_option (option);
	    option = parse_string (NULL);
	}
    }

    /* Parse the command line options */
    for (i = 1; i < argc; i++)
	parse_option (argv[i]);

    if (strlen(inname) == 0)
	abortmsg ("No input file specified", 1);

    if (strlen(outname) == 0)
	strcpy (outname, inname);

    switch (format) {
	case POV10:
	case POV20:   add_ext (outname, "pov", 1); break;
	case VIVID:   add_ext (outname, "v",   1); break;
	case POLYRAY: add_ext (outname, "pi",  1); break;
	case MGF:     add_ext (outname, "mgf", 1); break;
	case RAW:     add_ext (outname, "raw", 1); break;
    }

    switch (internal_bounding) {
	case OFF:  bound = 2; break;
	case ON:   bound = 0; break;
	case AUTO: bound = (format == POV10) ? 0 : 2; break;
    }

    if ((strlen(vuename) > 0) != (frame >= 0))
	abortmsg ("The -a and -f parameters must be used together", 1);

    if (format == RAW || (format == MGF && smooth < 0.1))
	no_opt = TRUE;
}


void parse_option (char *option)
{
    List *excl, *box;
    char name[80];

    if (option[0] == '-' || option[0] == '+') {
	switch (upcase(option[1])) {
	    case 'A': strcpy (vuename, &option[2]);
		      break;

	    case 'B': strcpy (name, parse_string (&option[2]));
		      if (strlen(name) == 0)
			  box_all = TRUE;
		      else {
			  cleanup_name (name);

			  box = malloc (sizeof (*box));
			  strcpy (box->name, name);

			  LIST_INSERT (box_list, box);
		      }
		      break;

	    case 'F': if (option[2] != '\0')
			  frame = atoi (&option[2]);
		      break;

	    case 'I': if (option[0] == '-')
			  internal_bounding = OFF;
		      else
			  internal_bounding = ON;
		      break;

	    case 'L': if (libs == MAX_LIB)
			  abortmsg ("Too many libraries specified", 1);

		      strcpy (libname[libs++], &option[2]);
		      break;

	    case 'O': switch (upcase(option[2])) {
			  case 'P': if (option[3] == '1')
					format = POV10;
				    else
					format = POV20;
				    break;

			  case 'V': format = VIVID;
				    break;

			  case 'L': format = POLYRAY;
				    break;

			  case 'R': format = RAW;
				    break;

			  case 'M': format = MGF;
				    break;

			  default:  printf ("Invalid output format %s specified\n", option);
				    exit(1);
		      }
		      break;

	    case 'S': if (option[2] != '\0')
			  smooth = atof (&option[2]);
		      break;

	    case 'U': printf ("Warning: -u parameter no long has any effect\n");
		      printf ("         use +i or -i instead.\n");
		      break;

	    case 'V': if (option[0] == '-')
			  verbose = 0;
		      else
			  verbose = 1;
		      break;

	    case 'X': strcpy (name, parse_string (&option[2]));
		      cleanup_name (name);

		      excl = malloc (sizeof (*excl));
		      strcpy (excl->name, name);

		      LIST_INSERT (excl_list, excl);
		      break;

	    default : printf ("\nInvalid option %s specified\n", option);
		      exit (1);
	}
    }
    else if (strlen (inname) == 0) {
	strcpy (inname, option);
	add_ext (inname, "3ds", 0);
    }
    else if (strlen (outname) == 0)
	strcpy (outname, option);
    else
	abortmsg ("Too many file names specified.\n", 1);
}


/* Insert a new node into the list */
void list_insert (List **root, List *new_node)
{
    new_node->next = *root;

    *root = new_node;
}


/* Find the node with the specified name */
void *list_find (List **root, char *name)
{
    List *p;

    for (p = *root; p != NULL; p = p->next) {
	if (strcmp (p->name, name) == 0)
	    break;
    }

    return (void *)p;
}


/* Delete the indicated node from the list */
void list_delete (List **root, List *node)
{
    List *prev;

    prev = *root;
    while (prev != NULL && prev->next != node)
	prev = prev->next;

    if (prev == NULL)
	*root = node->next;
    else
	prev->next = node->next;

    free (node);
}


/* Delete the entire list */
void list_kill (List **root)
{
    List *temp;

    while (*root != NULL) {
	temp = *root;
	*root = (*root)->next;
	free (temp);
    }
}


/* Add a new material to the material list */
Material *update_materials (char *new_material, int ext)
{
    Material *p;

    p = LIST_FIND (mtl_list, new_material);

    if (p == NULL) {
	p = malloc (sizeof (*p));

	if (p == NULL)
	    abortmsg ("Out of memory adding material", 1);

	strcpy (p->name, new_material);
	p->external = ext;

	LIST_INSERT (mtl_list, p);
    }

    return p;
}


MatProp *create_mprop()
{
    MatProp *new_mprop;

    new_mprop = malloc (sizeof(*new_mprop));
    if (new_mprop == NULL)
	abortmsg ("Out of memory adding material", 1);

    strcpy (new_mprop->name, "");
    new_mprop->ambient = Black;
    new_mprop->diffuse = Black;
    new_mprop->specular = Black;
    new_mprop->shininess = 0.0;
    new_mprop->transparency = 0.0;
    new_mprop->reflection = 0.0;
    new_mprop->self_illum = FALSE;
    new_mprop->two_side = FALSE;

    strcpy (new_mprop->tex_map, "");
    new_mprop->tex_strength = 0.0;

    strcpy (new_mprop->bump_map, "");
    new_mprop->bump_strength = 0.0;

    return new_mprop;
}


/* Load in any predefined materials */
void read_library (char *fname)
{
    FILE *lib;
    char string[256], name[80];

    if ((lib = fopen (fname, "r")) == NULL) {
	printf ("Cannot open texture library file %s!\n", fname);
	exit(1);
    }

    switch (format) {
	case POV10:
	case POV20:
	    while (fgets (string, 256, lib) != NULL) {
		if (strstr (string, "#declare")) {
		    strcpy (name, between (string, "#declare", "="));
		    cleanup_name (name);
		    (void)update_materials (name, TRUE);
		}
	    }
	    break;

	case VIVID:
	    while (fgets (string, 256, lib) != NULL) {
		if (strstr (string, "#define")) {
		    (void)parse_string (string);
		    strcpy (name, parse_string (NULL));
		    cleanup_name (name);
		    (void)update_materials (name, TRUE);
		}
	    }
	    break;

	case POLYRAY:
	    while (fgets (string, 256, lib) != NULL) {
		if (strstr (string, "define")) {
		    (void)parse_string (string);
		    strcpy (name, parse_string (NULL));
		    cleanup_name (name);
		    (void)update_materials (name, TRUE);
		}
	    }
	    break;

	case MGF:
	    while (read_mgfmatname(name, 80, lib))
		(void)update_materials (name, TRUE);
	    break;
    }

    fclose (lib);
}


/* parse the next MGF material name from f, return FALSE if EOF */
int read_mgfmatname (char *s, int n, FILE *f)
{	
    char inpline[128];
    register char *cp, *cp2;
    register int i;
					/* find material */
    while (fgets(inpline, sizeof(inpline), f) != NULL) {
	for (cp = inpline; isspace(*cp); cp++)
	    ;
	if (*cp++ != 'm' || !isspace(*cp++))
	    continue;
	while (isspace(*cp))
	    cp++;
	if (!*cp)
	    continue;
	for (i=n, cp2=s; *cp && !isspace(*cp); cp++)	/* get name */
	    if (--i > 0)
		*cp2++ = *cp;
	*cp2 = '\0';
	while (isspace(*cp))
	    cp++;
	if (*cp++ != '=' || !isspace(*cp))		/* not defined? */
	    continue;
	return TRUE;
    }
    return FALSE;
}


void write_intro (FILE *f)
{
    int i;

    switch (format) {
	case POV10:
	case POV20:
	    fprintf (f, "#include \"colors.inc\"\n");
	    fprintf (f, "#include \"shapes.inc\"\n");
	    fprintf (f, "#include \"textures.inc\"\n");

	    for (i = 0; i < MAX_LIB; i++) {
		if (strlen(libname[i]) > 0)
		    fprintf (f, "#include \"%s\"\n", libname[i]);
	    }

	    fprintf (f, "\n");
	    break;

	case VIVID:
	    fprintf (f, "#include color.vc\n");

	    for (i = 0; i < MAX_LIB; i++) {
		if (strlen(libname[i]) > 0)
		    fprintf (f, "#include %s\n", libname[i]);
	    }

	    fprintf (f, "\n");
	    break;

	case POLYRAY:
	    fprintf (f, "include \"colors.inc\"\n");

	    for (i = 0; i < MAX_LIB; i++) {
		if (strlen(libname[i]) > 0)
		    fprintf (f, "include \"%s\"\n", libname[i]);
	    }

	    fprintf (f, "\n");
	    break;

	case MGF:
	    fprintf (f, "c R =\n\tcxy %.3f %.3f\n", CIE_x_r, CIE_y_r);
	    fprintf (f, "c G =\n\tcxy %.3f %.3f\n", CIE_x_g, CIE_y_g);
	    fprintf (f, "c B =\n\tcxy %.3f %.3f\n", CIE_x_b, CIE_y_b);

	    for (i = 0; i < MAX_LIB; i++) {
		if (strlen(libname[i]) > 0)
		    fprintf (f, "i %s\n", libname[i]);
	    }

	    fprintf (f, "\n");
	    break;
    }
}


/* Write the object summary */
void write_summary (FILE *f)
{
    char *comstr;
    Summary *s;

    if (summary == NULL)
        return;

    switch (format) {
	case POV10:
	case POV20:
	case VIVID:
	case POLYRAY:
	    comstr = "//";
	    break;
	case MGF:
	    comstr = "# ";
	    break;
	default:
		printf ("Illegal format in write_summary() '%c'\n", format);
		exit(1);
    }
    fprintf (f, "%s   Object    CenterX    CenterY    CenterZ    LengthX    LengthY    LengthZ\n", comstr);
    fprintf (f, "%s ---------- ---------- ---------- ---------- ---------- ---------- ----------\n", comstr);

    for (s = summary; s != NULL; s = s->next) {
        fprintf (f, "%s %-10s%11.2f%11.2f%11.2f%11.2f%11.2f%11.2f\n",
		 comstr, s->name, s->center[X], s->center[Y], s->center[Z],
		 s->lengths[X], s->lengths[Y], s->lengths[Z]);
    }

    fprintf (f, "\n");
}


/* Write background solid colour */
void write_bgsolid (FILE *f, Colour col)
{
    switch (format) {
	case POV10:
	    fprintf (f, "/* Background colour */\n");
	    fprintf (f, "object {\n");
	    fprintf (f, "   sphere { <0.0 0.0 0.0> 1e6 }\n");
	    fprintf (f, "   texture {\n");
	    fprintf (f, "      ambient 1.0\n");
	    fprintf (f, "      diffuse 0.0\n");
	    fprintf (f, "      color red %4.2f green %4.2f blue %4.2f\n",
				 col.red, col.green, col.blue);
	    fprintf (f, "   }\n");
	    fprintf (f, "}\n\n");
	    break;

	case POV20:
	    fprintf (f, "background { color red %4.2f green %4.2f blue %4.2f }\n\n",
			   col.red, col.green, col.blue);
	    break;

	case POLYRAY:
	    fprintf (f, "background <%4.2f, %4.2f, %4.2f>\n\n",
			   col.red, col.green, col.blue);
	    break;
    }
}


void write_light (FILE *f, char *name, Vector pos, Colour col)
{
    switch (format) {
	case POV10:
	    fprintf (f, "/* Light: %s */\n", name);
	    fprintf (f, "object {\n");
	    fprintf (f, "    light_source { <%.4f %.4f %.4f> color red %4.2f green %4.2f blue %4.2f }\n",
			pos[X], pos[Y], pos[Z], col.red, col.green, col.blue);
	    fprintf (f, "}\n\n");
	    break;

	case POV20:
	    fprintf (f, "/* Light: %s */\n", name);
	    fprintf (f, "light_source {\n");
	    fprintf (f, "    <%.4f, %.4f, %.4f> color rgb <%4.2f, %4.2f, %4.2f>\n",
		     pos[X], pos[Y], pos[Z], col.red, col.green, col.blue);
	    fprintf (f, "}\n\n");
	    break;

	case VIVID:
	    fprintf (f, "/* Light: %s */\n", name);
	    fprintf (f, "light {\n");
	    fprintf (f, "    type point\n");
	    fprintf (f, "    position %.4f %.4f %.4f\n",
			     pos[X], pos[Y], pos[Z]);
	    fprintf (f, "    color %4.2f %4.2f %4.2f\n",
			     col.red, col.green, col.blue);
	    fprintf (f, "}\n\n");
	    break;

	case POLYRAY:
	    fprintf (f, "// Light: %s\n", name);
	    fprintf (f, "light <%4.2f, %4.2f, %4.2f>, <%.4f, %.4f, %.4f>\n\n",
			 col.red, col.green, col.blue, pos[X], pos[Y], pos[Z]);
	    break;

	case MGF:
	    fprintf (f, "\n# Light\n");
	    if (name[0]) fprintf (f, "o %s\n", name);
	    fprintf (f, "m\n\tsides 1\n\tc\n\t\t\tcmix %.3f R %.3f G %.3f B\n\ted %e\n",
		    CIE_Y_r*col.red, CIE_Y_g*col.green, CIE_Y_b*col.blue,
	    100000.0*(CIE_Y_r*col.red + CIE_Y_g*col.green + CIE_Y_b*col.blue));
	    fprintf (f, "v c =\n\tp %.4f %.4f %.4f\nsph c .01\n",
		    pos[X], pos[Y], pos[Z]);
	    if (name[0]) fprintf (f, "o\n");
	    fprintf (f, "\n");
	    break;
    }
}


void write_spot (FILE *f, char *name, Vector pos, Vector target, Colour col,
			  float hotspot, float falloff)
{
    switch (format) {
	case POV10:
	    fprintf (f, "/* Spotlight: %s */\n", name);
	    fprintf (f, "object {\n");
	    fprintf (f, "    light_source {\n");
	    fprintf (f, "        <%.4f %.4f %.4f> color red %4.2f green %4.2f blue %4.2f\n",
				   pos[X], pos[Y], pos[Z],
				   col.red, col.green, col.blue);
	    fprintf (f, "        spotlight\n");
	    fprintf (f, "        point_at <%.4f %.4f %.4f>\n",
				   target[X], target[Y], target[Z]);
	    fprintf (f, "        tightness 0\n");
	    fprintf (f, "        radius %.2f\n", 0.5*hotspot);
	    fprintf (f, "        falloff %.2f\n", 0.5*falloff);
	    fprintf (f, "    }\n");
	    fprintf (f, "}\n\n");
	    break;

	case POV20:
	    fprintf (f, "/* Spotlight: %s */\n", name);
	    fprintf (f, "light_source {\n");
	    fprintf (f, "    <%.4f, %.4f, %.4f> color rgb <%4.2f, %4.2f, %4.2f>\n",
			     pos[X], pos[Y], pos[Z],
			     col.red, col.green, col.blue);
	    fprintf (f, "    spotlight\n");
	    fprintf (f, "    point_at <%.4f, %.4f, %.4f>\n",
			     target[X], target[Y], target[Z]);
	    fprintf (f, "    tightness 0\n");
	    fprintf (f, "    radius %.2f\n", 0.5*hotspot);
	    fprintf (f, "    falloff %.2f\n", 0.5*falloff);
	    fprintf (f, "}\n\n");
	    break;

	case VIVID:
	    fprintf (f, "/* Spotlight: %s */\n", name);
	    fprintf (f, "light {\n");
	    fprintf (f, "    type spot\n");
	    fprintf (f, "    position %.4f %.4f %.4f\n",
			     pos[X], pos[Y], pos[Z]);
	    fprintf (f, "    at %.4f %.4f %.4f\n",
			     target[X], target[Y], target[Z]);
	    fprintf (f, "    color %4.2f %4.2f %4.2f\n",
			     col.red, col.green, col.blue);
	    fprintf (f, "    min_angle %.2f\n", hotspot);
	    fprintf (f, "    max_angle %.2f\n", falloff);
	    fprintf (f, "}\n\n");
	    break;

	case POLYRAY:
	    fprintf (f, "// Spotlight: %s\n", name);
	    fprintf (f, "spot_light <%4.2f, %4.2f, %4.2f>, <%.4f, %.4f, %.4f>,\n",
			  col.red, col.green, col.blue, pos[X], pos[Y], pos[Z]);
	    fprintf (f, "           <%.4f, %.4f, %.4f>, 0.0, %.2f, %.2f\n\n",
			  target[X], target[Y], target[Z], hotspot/2.0, falloff/2.0);
	    break;

	case MGF:
	    fprintf (f, "\n# Spotlight\n");
	    if (name[0]) fprintf (f, "o %s\n", name);
	    fprintf (f, "# hotspot: %.2f\n# falloff: %.2f\n", hotspot, falloff);
	    fprintf (f, "m\n\tsides 1\n\tc\n\t\t\tcmix %.3f R %.3f G %.3f B\n\ted %e\n",
		    CIE_Y_r*col.red, CIE_Y_g*col.green, CIE_Y_b*col.blue,
	    100000.0*(CIE_Y_r*col.red + CIE_Y_g*col.green + CIE_Y_b*col.blue));
	    fprintf (f, "v c =\n\tp %.4f %.4f %.4f\n\tn %.4f %.4f %.4f\n",
		    pos[X], pos[Y], pos[Z],
		    target[X]-pos[X], target[Y]-pos[Y], target[Z]-pos[Z]);
	    fprintf (f, "ring c 0 .01\n");
	    if (name[0]) fprintf (f, "o\n");
	    fprintf (f, "\n");
	    break;
    }
}


void write_fog (FILE *f, Colour col, float dist)
{
    if (dist <= 0.0)
	return;

    switch (format) {
	case POV10:
	    fprintf (f, "fog {\n");
	    fprintf (f, "    color red %4.2f green %4.2f blue %4.2f %.4f\n",
			  col.red, col.green, col.blue, dist/2.0);
	    fprintf (f, "}\n\n");
	    break;

	case POV20:
	    fprintf (f, "fog {\n");
	    fprintf (f, "    color red %4.2f green %4.2f blue %4.2f distance %.4f\n",
			  col.red, col.green, col.blue, dist/2.0);
	    fprintf (f, "}\n\n");
	    break;
    }
}


void write_camera (FILE *f, char *name, Vector pos, Vector target,
			    float lens, float bank)
{
    float fov;

    cameras++;

    fov = findfov (lens);

    switch (format) {
	case POV10:
	    /* Comment out multiple cameras */
	    if (cameras > 1)
		fprintf (f, "/*\n");

	    fprintf (f, "/* Camera: %s */\n", name);
	    fprintf (f, "camera {\n");
	    fprintf (f, "   location <%.4f %.4f %.4f>\n",
			      pos[X], pos[Y], pos[Z]);
	    fprintf (f, "   direction <0 %.3f 0>\n", 0.60/tan(0.5*RAD(fov)) );
	    fprintf (f, "   up <0 0 1>\n");
	    fprintf (f, "   sky  <0 0 1>\n");
	    fprintf (f, "   right <%.3f 0 0>\n", ASPECT);
	    fprintf (f, "   look_at <%.4f %.4f %.4f>\n",
			      target[X], target[Y], target[Z]);
	    if (bank != 0.0)
		fprintf (f, "   /* Bank angle = %.2f */\n", bank);

	    fprintf (f, "}\n");

	    if (cameras > 1)
		fprintf (f, "*/\n");

	    fprintf (f, "\n");
	    break;

	case POV20:
	    /* Comment out multiple cameras */
	    if (cameras > 1)
		fprintf (f, "/*\n");

	    fprintf (f, "/* Camera: %s */\n", name);
	    fprintf (f, "camera {\n");
	    fprintf (f, "   location <%.4f, %.4f, %.4f>\n",
			      pos[X], pos[Y], pos[Z]);
	    fprintf (f, "   direction <0, %.3f, 0>\n", 0.60/tan(0.5*RAD(fov)) );
	    fprintf (f, "   up <0, 0, 1>\n");
	    fprintf (f, "   sky  <0, 0, 1>\n");
	    fprintf (f, "   right <%.3f, 0, 0>\n", ASPECT);
	    fprintf (f, "   look_at <%.4f, %.4f, %.4f>\n",
			      target[X], target[Y], target[Z]);
	    if (bank != 0.0)
		fprintf (f, "   /* Bank angle = %.2f */\n", bank);

	    fprintf (f, "}\n");

	    if (cameras > 1)
		fprintf (f, "*/\n");

	    fprintf (f, "\n");
	    break;

	case VIVID:
	    fprintf (f, "/* Camera: %s */\n", name);

	    if (cameras > 1)
		fprintf (f, "/*\n");

	    fprintf (f, "studio {\n");
	    fprintf (f, "    from %.4f %.4f %.4f\n",
			       pos[X], pos[Y], pos[Z]);
	    fprintf (f, "    at %.4f %.4f %.4f\n",
			       target[X], target[Y], target[Z]);
	    fprintf (f, "    up 0 0 1\n");
	    fprintf (f, "    angle %.2f\n", 1.1*fov);
	    fprintf (f, "    aspect %.3f\n", ASPECT);
	    fprintf (f, "    resolution 320 200\n");
	    fprintf (f, "    antialias none\n");
	    fprintf (f, "}\n");

	    if (cameras > 1)
		fprintf (f, "*/\n");

	    fprintf (f, "\n");
	    break;

	case POLYRAY:
	    if (cameras == 1) {
		fprintf (f, "// Camera: %s\n", name);
		fprintf (f, "viewpoint {\n");
		fprintf (f, "    from <%.4f, %.4f, %.4f>\n",
				   pos[X], pos[Y], pos[Z]);
		fprintf (f, "    at <%.4f, %.4f, %.4f>\n",
				   target[X], target[Y], target[Z]);
		fprintf (f, "    up <0, 0, 1>\n");
		fprintf (f, "    angle %.2f\n", 0.85*fov);
		fprintf (f, "    aspect %.3f\n", -(ASPECT));
		fprintf (f, "    resolution 320, 200\n");
		fprintf (f, "}\n");
	    }

	    fprintf (f, "\n");
	    break;

	case MGF:
	    fprintf (f, "# Camera %s\n", name);
	    fprintf (f, "# from: %.4f %.4f %.4f\n", pos[X], pos[Y], pos[Z]);
	    fprintf (f, "# to: %.4f %.4f %.4f\n", target[X], target[Y], target[Z]);
	    fprintf (f, "# lens length: %.2f\n", lens);
	    fprintf (f, "# bank: %.2f\n", bank);
	    break;
    }
}


void write_material (FILE *f, char *mat)
{
    MatProp *mprop = LIST_FIND (mprop_list, mat);

    if (mprop == NULL) {
       mprop = &DefaultMaterial;
       (void)strcpy(mprop->name, mat);
    }

    switch (format) {
	case POV10:
	    write_pov10_material (f, mprop);
	    break;

	case POV20:
	    write_pov20_material (f, mprop);
	    break;

	case VIVID:
	    write_vivid_material (f, mprop);
	    break;

	case POLYRAY:
	    write_polyray_material (f, mprop);
	    break;

	case MGF:
	    write_mgf_material (f, mprop);
	    break;
    }
}


void write_pov10_material (FILE *f, MatProp *m)
{
    float amb = 0.1, dif = 0.9, spec = 1.0;
    float dist_white, dist_diff, phong, phong_size;
    float red, green, blue;

    /* amb = get_ambient (m); */

    if (m->self_illum) {
	amb = 0.9;
	dif = 0.1;
    }

    dist_white = fabs(1.0 - m->specular.red) +
		 fabs(1.0 - m->specular.green) +
		 fabs(1.0 - m->specular.blue);

    dist_diff  = fabs(m->diffuse.red   - m->specular.red) +
		 fabs(m->diffuse.green - m->specular.green) +
		 fabs(m->diffuse.blue  - m->specular.blue);


    phong_size = 0.7*m->shininess;
    if (phong_size < 1.0) phong_size = 1.0;

    if (phong_size > 30.0)
	phong = 1.0;
    else
	phong = phong_size/30.0;

    fprintf (f, "#declare %s = texture {\n", m->name);
    fprintf (f, "    ambient %.2f\n", amb);
    fprintf (f, "    diffuse %.2f\n", dif);
    fprintf (f, "    phong %.2f\n", phong);
    fprintf (f, "    phong_size %.1f\n", phong_size);

    if (dist_diff < dist_white)
	fprintf (f, "    metallic\n");

    if (m->reflection > 0.0) {
	spec = (m->specular.red + m->specular.green + m->specular.blue)/3.0;
	fprintf (f, "    reflection %.3f\n", spec * m->reflection);
    }

    if (m->transparency > 0.0) {
	red   = m->diffuse.red;
	green = m->diffuse.green;
	blue  = m->diffuse.blue;

	/* Saturate the colour towards white as the transparency increases */
	red   = ((1.0 - m->transparency) * red)   + m->transparency;
	green = ((1.0 - m->transparency) * green) + m->transparency;
	blue  = ((1.0 - m->transparency) * blue)  + m->transparency;

	fprintf (f, "    color red %.3f green %.3f blue %.3f alpha %.3f\n",
			 red, green, blue, m->transparency);
	fprintf (f, "    ior 1.1\n");
	fprintf (f, "    refraction 1.0\n");
    }
    else
	fprintf (f, "    color red %.3f green %.3f blue %.3f\n",
			 m->diffuse.red, m->diffuse.green, m->diffuse.blue);

    if (strlen (m->tex_map) > 0) {
	fprintf (f, "    /* Image map: %s, Strength: %.2f */\n",
			 m->tex_map, m->tex_strength);
    }

    if (strlen (m->bump_map) > 0) {
	fprintf (f, "    /* Bump map: %s, Strength: %.2f */\n",
			 m->bump_map, m->bump_strength);
    }

    fprintf (f, "}\n\n");
}


void write_pov20_material (FILE *f, MatProp *m)
{
    float amb = 0.1, dif = 0.9, spec = 1.0;
    float dist_white, dist_diff, phong, phong_size;
    float red, green, blue;

    /* amb = get_ambient (m); */

    if (m->self_illum) {
	amb = 0.9;
	dif = 0.1;
    }

    dist_white = fabs(1.0 - m->specular.red) +
		 fabs(1.0 - m->specular.green) +
		 fabs(1.0 - m->specular.blue);

    dist_diff  = fabs(m->diffuse.red   - m->specular.red) +
		 fabs(m->diffuse.green - m->specular.green) +
		 fabs(m->diffuse.blue  - m->specular.blue);

    phong_size = 0.7*m->shininess;
    if (phong_size < 1.0) phong_size = 1.0;

    if (phong_size > 30.0)
	phong = 1.0;
    else
	phong = phong_size/30.0;

    fprintf (f, "#declare %s = texture {\n", m->name);
    fprintf (f, "    finish {\n");
    fprintf (f, "        ambient %.2f\n", amb);
    fprintf (f, "        diffuse %.2f\n", dif);
    fprintf (f, "        phong %.2f\n", phong);
    fprintf (f, "        phong_size %.1f\n", phong_size);

    if (dist_diff < dist_white)
	fprintf (f, "        metallic\n");

    if (m->reflection > 0.0) {
	spec = (m->specular.red + m->specular.green + m->specular.blue)/3.0;
	fprintf (f, "        reflection %.3f\n", spec * m->reflection);
    }

    if (m->transparency > 0.0) {
	fprintf (f, "        ior 1.1\n");
	fprintf (f, "        refraction 1.0\n");
    }

    fprintf (f, "    }\n");

    if (m->transparency > 0.0) {
	red   = m->diffuse.red;
	green = m->diffuse.green;
	blue  = m->diffuse.blue;

	/* Saturate the colour towards white as the transparency increases */
	red   = ((1.0 - m->transparency) * red)   + m->transparency;
	green = ((1.0 - m->transparency) * green) + m->transparency;
	blue  = ((1.0 - m->transparency) * blue)  + m->transparency;

	fprintf (f, "    pigment { rgbf <%.3f, %.3f, %.3f, %.3f> }\n",
		    red, green, blue, m->transparency);
    }
    else
	fprintf (f, "    pigment { rgb <%.3f, %.3f, %.3f> }\n",
		    m->diffuse.red, m->diffuse.green, m->diffuse.blue);

    if (strlen (m->tex_map) > 0) {
	fprintf (f, "    /* Image map: %s, Strength: %.2f */\n",
			 m->tex_map, m->tex_strength);
    }

    if (strlen (m->bump_map) > 0) {
	fprintf (f, "    /* Bump map: %s, Strength: %.2f */\n",
			 m->bump_map, m->bump_strength);
    }

    fprintf (f, "}\n\n");
}


void write_vivid_material (FILE *f, MatProp *m)
{
    float amb = 0.1, dif = 0.9;

    /* amb = get_ambient (m); */

    if (m->self_illum) {
	amb = 0.9;
	dif = 0.1;
    }

    if (m->transparency > 0.0) {
       dif = dif - m->transparency;
       if (dif < 0.0) dif = 0.0;
    }

    fprintf (f, "#define %s \\ \n", m->name);
    fprintf (f, "    surface {           \\ \n");
    fprintf (f, "        ambient %.3f %.3f %.3f \\ \n",
			 amb*m->ambient.red, amb*m->ambient.green, amb*m->ambient.blue);

    fprintf (f, "        diffuse %.3f %.3f %.3f \\ \n",
			 dif*m->diffuse.red, dif*m->diffuse.green, dif*m->diffuse.blue);

    fprintf (f, "        shine %.1f  %.3f %.3f %.3f \\ \n",
			 0.7*m->shininess, m->specular.red, m->specular.green, m->specular.blue);

    if (m->transparency > 0.0) {
	fprintf (f, "        transparent %.3f*white \\ \n", 1.0 - (1.0 - m->transparency)/14.0);
	fprintf (f, "        ior 1.1 \\ \n");
    }

    if (m->reflection > 0.0)
	fprintf (f, "        specular %.3f*white \\ \n", m->reflection);

    if (strlen (m->tex_map) > 0) {
	fprintf (f, "        /* Image map: %s, Strength: %.2f */ \\ \n",
			     m->tex_map, m->tex_strength);
    }

    if (strlen (m->bump_map) > 0) {
	fprintf (f, "        /* Bump map: %s, Strength: %.2f */ \\ \n",
			     m->bump_map, m->bump_strength);
    }

    fprintf (f, "    }\n\n");
}


void write_polyray_material (FILE *f, MatProp *m)
{
    float amb = 0.1, dif = 0.9, spec;

    /* amb = get_ambient (m); */

    if (m->self_illum) {
	amb = 0.9;
	dif = 0.1;
    }

    if (m->transparency > 0.0) {
       dif = dif - m->transparency;
       if (dif < 0.0) dif = 0.0;
    }

    if (m->shininess == 0.0)
	m->shininess = 0.1;

    if (m->shininess > 40.0)
	spec = 1.0;
    else
	spec = m->shininess/40.0;

    fprintf (f, "define %s\n", m->name);
    fprintf (f, "texture {\n");
    fprintf (f, "    surface {\n");
    fprintf (f, "        ambient <%.3f, %.3f, %.3f>, %.1f\n",
			 m->ambient.red, m->ambient.green, m->ambient.blue, amb);

    fprintf (f, "        diffuse <%.3f, %.3f, %.3f>, %.1f\n",
			 m->diffuse.red, m->diffuse.green, m->diffuse.blue, dif);

    fprintf (f, "        specular <%.3f, %.3f, %.3f>, %.2f\n",
			 m->specular.red, m->specular.green, m->specular.blue, spec);

    fprintf (f, "        microfacet Reitz %.1f\n", 400.0/m->shininess);

    if (m->transparency > 0.0)
	fprintf (f, "        transmission %.3f, 1.1\n", m->transparency);

    if (m->reflection > 0.0)
	fprintf (f, "        reflection %.3f\n", m->reflection);

    if (strlen (m->tex_map) > 0) {
	fprintf (f, "        // Image map: %s, Strength: %.2f\n",
			     m->tex_map, m->tex_strength);
    }

    if (strlen (m->bump_map) > 0) {
	fprintf (f, "        // Bump map: %s, Strength: %.2f\n",
			     m->bump_map, m->bump_strength);
    }

    fprintf (f, "    }\n");
    fprintf (f, "}\n\n");
}


void write_mgf_material (FILE *f, MatProp *m)
{
    float  dmag, smag, rdmag, rsmag, tdmag, tsmag, total;

    fprintf (f, "m %s =\n", m->name);
    fprintf (f, "\tsides %d\n", m->two_side ? 2 : 1);
    dmag = CIE_Y_r*m->diffuse.red + CIE_Y_g*m->diffuse.green
		+ CIE_Y_b*m->diffuse.blue;
    smag = CIE_Y_r*m->specular.red + CIE_Y_g*m->specular.green
		+ CIE_Y_b*m->specular.blue;
    rdmag = dmag;
    rsmag = smag * m->reflection;
    tdmag = 0.0;
    tsmag = m->transparency;
    total = rdmag + rsmag + tdmag + tsmag;
    if (total > 0.99) {
	total = 0.9/total;
	dmag *= total;
	smag *= total;
	rdmag *= total;
	rsmag *= total;
	tdmag *= total;
	tsmag *= total;
	total = 0.9;
    }
    if (dmag > 0.005) {
	fprintf (f, "\tc\n\t\tcmix %.3f R %.3f G %.3f B\n",
		CIE_Y_r*m->diffuse.red,
		CIE_Y_g*m->diffuse.green,
		CIE_Y_b*m->diffuse.blue);
	if (rdmag > 0.005)
	    fprintf (f, "\trd %.4f\n", rdmag);
	if (tdmag > 0.005)
	    fprintf (f, "\ttd %.4f\n", tdmag);
	if (m->self_illum)
	    fprintf (f, "\ted %.4f\n", dmag);
    }
    if (m->shininess > 1.1 && rsmag > 0.005) {
	fprintf (f, "\tc\n\t\tcmix %.3f R %.3f G %.3f B\n",
		CIE_Y_r*m->specular.red,
		CIE_Y_g*m->specular.green,
		CIE_Y_b*m->specular.blue);
	fprintf (f, "\trs %.4f %.4f\n", rsmag, 0.6/sqrt(m->shininess));
    }
    if (tsmag > 0.005)
	fprintf (f, "\tc\n\tts %.4f 0\n", tsmag);

    if (strlen (m->tex_map) > 0) {
	fprintf (f, "# image map: %s, strength: %.2f\n",
			 m->tex_map, m->tex_strength);
    }

    if (strlen (m->bump_map) > 0) {
	fprintf (f, "# bump map: %s, strength: %.2f\n",
			 m->bump_map, m->bump_strength);
    }

    fprintf (f, "\n");
}


/* Write a mesh file */
void write_mesh (FILE *f, Mesh *mesh)
{
    int i;
    char curmat[80];
    Vector va, vb, vc;
    Summary *new_summary;
    Matrix obj_matrix;

    if (mesh->hidden || LIST_FIND (excl_list, mesh->name))
	return;

    /* Add this object's stats to the summary */
    new_summary = malloc (sizeof(*new_summary));
    if (new_summary == NULL)
	abortmsg ("Out of memory adding summary", 1);

    strcpy (new_summary->name, mesh->name);
    vect_copy (new_summary->center,  mesh->center);
    vect_copy (new_summary->lengths, mesh->lengths);

    LIST_INSERT (summary, new_summary);

    /* Compute the object transformation matrix for animations */
    if (ani_matrix != NULL) {
	mat_copy (obj_matrix, *ani_matrix);
	if (vue_version > 2.0)
	    mat_mult (obj_matrix, mesh->invmatrix, obj_matrix);
    }

    switch (format) {
	case MGF:
	    if (no_opt) {
		if (mesh->name[0]) fprintf (meshf, "o %s\n", mesh->name);
		for (i = 0; i < mesh->vertices; i++) {
		    vect_copy(va, mesh->vertex[i]);
		    if (ani_matrix != NULL)
			vect_transform (va, va, obj_matrix);
		    fprintf (meshf, "v v%d =\n\tp %.5f %.5f %.5f\n",
				i, va[X], va[Y], va[Z]);
		}
		curmat[0] = '\0';
		for (i = 0; i < mesh->faces; i++) {
		    if (strcmp(mesh->mtl[i]->name, curmat)) {
			strcpy(curmat, mesh->mtl[i]->name);
			fprintf (meshf, "m %s\n", curmat);
		    }
		    fprintf (meshf, "f v%d v%d v%d\n", mesh->face[i].a,
				mesh->face[i].b, mesh->face[i].c);
		}
		if (mesh->name[0]) fprintf (meshf, "o\n");
		break;
	    }
	    /*FALL THROUGH*/
	case POV10:
	case POV20:
	case VIVID:
	case POLYRAY:
	    opt_set_vert (mesh->vertices);

	    for (i = 0; i < mesh->faces; i++) {
		vect_copy (va, mesh->vertex[mesh->face[i].a]);
		vect_copy (vb, mesh->vertex[mesh->face[i].b]);
		vect_copy (vc, mesh->vertex[mesh->face[i].c]);

		opt_set_texture (mesh->mtl[i]->name);

		opt_add_tri (va[X], va[Y], va[Z], vc[X], vc[Y], vc[Z],
			     vb[X], vb[Y], vb[Z]);
	    }

	    fflush (f);

	    if (ani_matrix != NULL)
		opt_set_transform (obj_matrix);

	    if (box_all || LIST_FIND (box_list, mesh->name))
		opt_write_box (mesh->name);
	    else
		opt_write_file (mesh->name);

	    break;

	case RAW:
	    fprintf (f, "%s\n", mesh->name);

	    for (i = 0; i < mesh->faces; i++) {
		vect_copy (va, mesh->vertex[mesh->face[i].a]);
		vect_copy (vb, mesh->vertex[mesh->face[i].b]);
		vect_copy (vc, mesh->vertex[mesh->face[i].c]);

		if (ani_matrix != NULL) {
		    vect_transform (va, va, obj_matrix);
		    vect_transform (vb, vb, obj_matrix);
		    vect_transform (vc, vc, obj_matrix);
		}

		fprintf (f, "%f %f %f   %f %f %f   %f %f %f\n",
			    va[X], va[Y], va[Z], vb[X], vb[Y], vb[Z],
			    vc[X], vc[Y], vc[Z]);
	    }

	    break;
    }
}


/* Parses an object transformation and returns a pointer to the
   newly allocated transformation */
Transform *parse_transform (char *string)
{
    Transform *t;
    char      *token;
    int       token_no;

    t = (Transform *)malloc (sizeof(*t));
    if (t == NULL)
	abortmsg ("Out of memory allocating transform", 1);

    mat_identity (t->matrix);

    token = parse_string (string);
    token_no = 0;

    while (strlen(token) > 0) {
	 switch (token_no) {
	     case  0: break;
	     case  1: strcpy (t->name, token); break;
	     case  2: t->matrix[0][0] = atof(token); break;
	     case  3: t->matrix[0][1] = atof(token); break;
	     case  4: t->matrix[0][2] = atof(token); break;
	     case  5: t->matrix[1][0] = atof(token); break;
	     case  6: t->matrix[1][1] = atof(token); break;
	     case  7: t->matrix[1][2] = atof(token); break;
	     case  8: t->matrix[2][0] = atof(token); break;
	     case  9: t->matrix[2][1] = atof(token); break;
	     case 10: t->matrix[2][2] = atof(token); break;
	     case 11: t->matrix[3][0] = atof(token); break;
	     case 12: t->matrix[3][1] = atof(token); break;
	     case 13: t->matrix[3][2] = atof(token); break;

	     default: abortmsg ("Error parsing transform", 1);
	 }

	 token = parse_string (NULL);
	 token_no++;
    }

    t->matrix[0][3] = 0.0;
    t->matrix[1][3] = 0.0;
    t->matrix[2][3] = 0.0;
    t->matrix[3][3] = 1.0;

    cleanup_name (t->name);

    return t;
}


/* Parses a morph command and returns a pointer to the
   newly allocated morph */
Morph *parse_morph (char *string)
{
    Morph  *m;
    char   *token;
    int    i, token_no;

    m = (Morph *)malloc (sizeof(*m));
    if (m == NULL)
	abortmsg ("Out of memory allocating morph", 1);

    mat_identity (m->matrix);

    token = parse_string (string);

    token = parse_string (NULL);
    strcpy (m->name, token);

    token = parse_string (NULL);
    m->count = atoi (token);

    if (strlen (m->name) == 0 || m->count < 1 || m->count > 4)
	abortmsg ("Error parsing morph command", 1);

    cleanup_name (m->name);

    for (i = 0; i < m->count; i++) {
	token = parse_string (NULL);
	strcpy (m->names[i], token);

	token = parse_string (NULL);
	m->weight[i] = atof (token);

	if (strlen (m->names[i]) == 0)
	    abortmsg ("Error parsing morph command", 1);

	cleanup_name (m->names[i]);
    }

    token = parse_string (NULL);
    token_no = 0;

    while (strlen(token) > 0) {
	 switch (token_no) {
	     case  0: m->matrix[0][0] = atof(token); break;
	     case  1: m->matrix[0][1] = atof(token); break;
	     case  2: m->matrix[0][2] = atof(token); break;
	     case  3: m->matrix[1][0] = atof(token); break;
	     case  4: m->matrix[1][1] = atof(token); break;
	     case  5: m->matrix[1][2] = atof(token); break;
	     case  6: m->matrix[2][0] = atof(token); break;
	     case  7: m->matrix[2][1] = atof(token); break;
	     case  8: m->matrix[2][2] = atof(token); break;
	     case  9: m->matrix[3][0] = atof(token); break;
	     case 10: m->matrix[3][1] = atof(token); break;
	     case 11: m->matrix[3][2] = atof(token); break;

	     default: abortmsg ("Error parsing morph command", 1);
	 }

	 token = parse_string (NULL);
	 token_no++;
    }

    m->matrix[0][3] = 0.0;
    m->matrix[1][3] = 0.0;
    m->matrix[2][3] = 0.0;
    m->matrix[3][3] = 1.0;

    return m;
}


/* Parses an omni light and returns a pointer to the
   newly allocated light */
OmniLight *parse_omnilight (char *string)
{
    OmniLight *o;
    char      *token;
    int       token_no;

    o = (OmniLight *)malloc (sizeof(*o));
    if (o == NULL)
	abortmsg ("Out of memory allocating omnilight", 1);

    token = parse_string (string);
    token_no = 0;

    while (strlen(token) > 0) {
	 switch (token_no) {
	     case 0: break;
	     case 1: strcpy (o->name, token); break;
	     case 2: o->pos[X] = atof (token); break;
	     case 3: o->pos[Y] = atof (token); break;
	     case 4: o->pos[Z] = atof (token); break;
	     case 5: o->col.red   = atof (token); break;
	     case 6: o->col.green = atof (token); break;
	     case 7: o->col.blue  = atof (token); break;

	     default: abortmsg ("Error parsing omnilight", 1);
	 }

	 token = parse_string (NULL);
	 token_no++;
    }

    cleanup_name (o->name);

    return o;
}


/* Parses a spotlight and returns a pointer to the
   newly allocated spotlight */
Spotlight *parse_spotlight (char *string)
{
    Spotlight *s;
    char      *token;
    int       token_no;

    s = (Spotlight *)malloc (sizeof(*s));
    if (s == NULL)
	abortmsg ("Out of memory allocating spotlight", 1);

    token = parse_string (string);
    token_no = 0;

    while (strlen(token) > 0) {
	 switch (token_no) {
	     case  0: break;
	     case  1: strcpy (s->name, token); break;
	     case  2: s->pos[X] = atof (token); break;
	     case  3: s->pos[Y] = atof (token); break;
	     case  4: s->pos[Z] = atof (token); break;
	     case  5: s->target[X] = atof (token); break;
	     case  6: s->target[Y] = atof (token); break;
	     case  7: s->target[Z] = atof (token); break;
	     case  8: s->col.red   = atof (token); break;
	     case  9: s->col.green = atof (token); break;
	     case 10: s->col.blue  = atof (token); break;
	     case 11: s->hotspot   = atof (token); break;
	     case 12: s->falloff   = atof (token); break;
	     case 13: break;

	     default: abortmsg ("Error parsing spotlight", 1);
	 }

	 token = parse_string (NULL);
	 token_no++;
    }

    cleanup_name (s->name);

    return s;
}


/* Parses a camera command and returns a pointer to the
   newly allocated camera */
Camera *parse_camera (char *string)
{
    Camera *c;
    char   *token;
    int    token_no;

    c = (Camera *)malloc (sizeof(*c));
    if (c == NULL)
	abortmsg ("Out of memory allocating camera", 1);

    token = parse_string (string);
    token_no = 0;

    while (strlen(token) > 0) {
	 switch (token_no) {
	     case 0: break;
	     case 1: c->pos[X] = atof (token); break;
	     case 2: c->pos[Y] = atof (token); break;
	     case 3: c->pos[Z] = atof (token); break;
	     case 4: c->target[X] = atof (token); break;
	     case 5: c->target[Y] = atof (token); break;
	     case 6: c->target[Z] = atof (token); break;
	     case 7: c->bank = atof (token); break;
	     case 8: c->lens = atof (token); break;

	     default: abortmsg ("Error parsing camera", 1);
	 }

	 token = parse_string (NULL);
	 token_no++;
    }

    return c;
}


/* Load the transforms, camera movements, etc for the specified frame */
void read_frame (char *filename, int frame_no)
{
    FILE  *f;
    char  fname[80];
    char  string[256];
    char  *token;

    /* Open the .vue file */
    strcpy (fname, filename);   /* Make a copy we can mess with */
    add_ext (fname, "vue", 0);

    f = fopen (fname, "r");
    if (f == NULL) {
	printf ("Error opening file '%s'\n", fname);
	exit(1);
    }

    /* Load the specified frame */
    find_frame (f, frame_no);

    while (fgets (string, 256, f) != NULL) {
	token = parse_string (string);

	if (strcmp (token, "frame") == 0)
	    break;
	else if (strcmp (token, "transform") == 0) {
	    LIST_INSERT (trans_list, parse_transform (string));
	}
	else if (strcmp (token, "morph") == 0) {
	    LIST_INSERT (morph_list, parse_morph (string));
	}
	else if (strcmp (token, "light") == 0) {
	    LIST_INSERT (omni_list, parse_omnilight (string));
	}
	else if (strcmp (token, "spotlight") == 0) {
	    LIST_INSERT (spot_list, parse_spotlight (string));
	}
	else if (strcmp (token, "camera") == 0) {
	    if (cam_list != NULL)
		abortmsg ("ERROR - Multiple cameras in .vue file", 1);

	    LIST_INSERT (cam_list, parse_camera (string));
	}
	else if (strcmp (token, "top") == 0)
	    abortmsg ("ERROR - Orthogonal viewports are not supported", 1);
	else if (strcmp (token, "bottom") == 0)
	    abortmsg ("ERROR - Orthogonal viewports are not supported", 1);
	else if (strcmp (token, "left") == 0)
	    abortmsg ("ERROR - Orthogonal viewports are not supported", 1);
	else if (strcmp (token, "right") == 0)
	    abortmsg ("ERROR - Orthogonal viewports are not supported", 1);
	else if (strcmp (token, "front") == 0)
	    abortmsg ("ERROR - Orthogonal viewports are not supported", 1);
	else if (strcmp (token, "back") == 0)
	    abortmsg ("ERROR - Orthogonal viewports are not supported", 1);
	else if (strcmp (token, "user") == 0)
	    abortmsg ("ERROR - User viewports are not supported", 1);
    }

    fclose(f);
}


void find_frame (FILE *f, int frame_no)
{
    char  string[256];
    char  *token;
    int   frame = 0;

    /* Search the .vue file for the required frame */
    while (1) {
	/* Read the next line in the file */
	if (fgets (string, 256, f) == NULL) {
	    printf ("Unable to locate frame #%d in .vue file\n", frame_no);
	    exit(1);
	}

	token = parse_string (string);

	if (strcmp (token, "frame") == 0) {
	    token = parse_string (NULL);

	    if (strlen(token) == 0) {
		printf ("Unable to locate frame #%d in .vue file\n", frame_no);
		exit(1);
	    }

	    frame = atoi (token);

	    if (frame == frame_no)
		break;
	}
	else if (strcmp (token, "VERSION") == 0) {
	    token = parse_string (NULL);

	    vue_version = atoi(token) / 100.0;
	}
    }
}


void save_animation()
{
    Mesh      *mesh, *master;
    Transform *t;
    Morph     *m;
    Vector    temp;
    int       i, j;

    printf ("\n");

    for (t = trans_list; t != NULL; t = t->next) {
	 printf ("Transforming object: %s\n", t->name);

	 ani_matrix = &(t->matrix);

	 mesh = LIST_FIND (mesh_list, t->name);

	 if (mesh == NULL) {
	     printf ("Unable to locate mesh object %s\n", t->name);
	     exit(1);
	 }

	 write_mesh (out, mesh);
    }

    for (m = morph_list; m != NULL; m = m->next) {
	printf ("Morphing object: %s\n", m->name);

	ani_matrix = &(m->matrix);

	mesh = LIST_FIND (mesh_list, m->name);
	if (mesh == NULL) {
	    printf ("Unable to locate mesh object %s\n", m->name);
	    exit(1);
	}

	/* Make a copy to mess with */
	master = copy_mesh (mesh);
	master->hidden = FALSE;

	strcpy (master->name, m->name);

	for (i = 0; i < master->vertices; i++)
	    vect_init (master->vertex[i], 0.0, 0.0, 0.0);

	for (i = 0; i < m->count; i++) {
	    mesh = LIST_FIND (mesh_list, m->names[i]);
	    if (mesh == NULL) {
		printf ("Unable to locate mesh object %s\n", m->names[0]);
		exit(1);
	    }

	    if (mesh->vertices != master->vertices)
		abortmsg ("Morphed objects do not contain the same number of vertices", 1);

	    if (mesh->faces != master->faces)
		abortmsg ("Morphed objects do not contain the same number of faces", 1);

	    for (j = 0; j < master->vertices; j++) {
		vect_transform (temp, mesh->vertex[j], mesh->invmatrix);
		vect_scale (temp, temp, m->weight[i]);
		vect_add (master->vertex[j], master->vertex[j], temp);
	    }
	}

	for (i = 0; i < master->vertices; i++)
	    vect_transform (master->vertex[i], master->vertex[i], master->matrix);

	write_mesh (out, master);

	free_mesh_data (master);
	free (master);
    }

    for (mesh = mesh_list; mesh != NULL; mesh = mesh->next)
	free_mesh_data (mesh);
}


/* Create a new mesh */
Mesh *create_mesh (char *name, int vertices, int faces)
{
    Mesh *new_mesh;

    new_mesh = malloc (sizeof(*new_mesh));
    if (new_mesh == NULL)
	abortmsg ("Out of memory allocating mesh", 1);

    strcpy (new_mesh->name, name);

    new_mesh->vertices = vertices;

    if (vertices <= 0)
	new_mesh->vertex = NULL;
    else {
	new_mesh->vertex = malloc (vertices * sizeof(*new_mesh->vertex));
	if (new_mesh->vertex == NULL)
	    abortmsg ("Out of memory allocating mesh", 1);
    }

    new_mesh->faces = faces;

    if (faces <= 0) {
	new_mesh->face = NULL;
	new_mesh->mtl = NULL;
    }
    else {
	new_mesh->face = malloc (faces * sizeof(*new_mesh->face));
	if (new_mesh->face == NULL)
	    abortmsg ("Out of memory allocating mesh", 1);

	new_mesh->mtl = malloc (faces * sizeof(*new_mesh->mtl));
	if (new_mesh->mtl == NULL)
	    abortmsg ("Out of memory allocating mesh", 1);
    }

    vect_init (new_mesh->center,  0.0, 0.0, 0.0);
    vect_init (new_mesh->lengths, 0.0, 0.0, 0.0);

    mat_identity (new_mesh->matrix);
    mat_identity (new_mesh->invmatrix);

    new_mesh->hidden = FALSE;
    new_mesh->shadow = TRUE;

    return new_mesh;
}


/* Creates a duplicate copy of a mesh */
Mesh *copy_mesh (Mesh *mesh)
{
    Mesh *new_mesh;
    int  i;

    new_mesh = create_mesh (mesh->name, mesh->vertices, mesh->faces);

    if (new_mesh == NULL)
	abortmsg ("Out of memory allocating mesh", 1);

    for (i = 0; i < mesh->vertices; i++)
	vect_copy (new_mesh->vertex[i], mesh->vertex[i]);

    for (i = 0; i < mesh->faces; i++) {
	new_mesh->face[i] = mesh->face[i];
	new_mesh->mtl[i]  = mesh->mtl[i];
    }

    mat_copy (new_mesh->matrix, mesh->matrix);
    mat_copy (new_mesh->invmatrix, mesh->invmatrix);

    vect_copy (new_mesh->center,  mesh->center);
    vect_copy (new_mesh->lengths, mesh->lengths);

    new_mesh->hidden = mesh->hidden;
    new_mesh->shadow = mesh->shadow;

    return new_mesh;
}


/* Free all data associated with mesh object */
void free_mesh_data (Mesh *mesh)
{
    if (mesh->vertex != NULL)
	free (mesh->vertex);

    if (mesh->face != NULL)
	free (mesh->face);

    if (mesh->mtl != NULL)
	free (mesh->mtl);
}


/* Updates the center (pivot) point of the mesh */
void update_limits (Mesh *mesh)
{
    Vector vmin = {+MAXFLOAT, +MAXFLOAT, +MAXFLOAT};
    Vector vmax = {-MAXFLOAT, -MAXFLOAT, -MAXFLOAT};
    int    i;

    for (i = 0; i < mesh->vertices; i++) {
	vect_min (vmin, vmin, mesh->vertex[i]);
	vect_max (vmax, vmax, mesh->vertex[i]);
    }

    vect_add  (mesh->center, vmin, vmax);
    vect_scale (mesh->center, mesh->center, 0.5);

    vect_sub (mesh->lengths, vmax, vmin);
}


/* Return the sub-string of 'str' that is before 'target' */
char *before (char *str, char *target)
{
    static char result[256];
    char   *search;

    strncpy (result, str, 256);
    result[255] = '\0';

    search = strstr (result, target);

    if (search != NULL)
	*search = '\0';

    return result;
}


/* Return the sub-string of 'str' that is after 'target' */
char *after (char *str, char *target)
{
    static char result[256];
    char   *search;

    search = strstr (str, target);

    if (search == NULL)
	strncpy (result, "", 256);
    else
	strncpy (result, search + strlen(target), 256);

    result[255] = '\0';

    return result;
}


/* Return the sub-string of 'str' that is between 'target1' and 'target2' */
char *between (char *str, char *target1, char *target2)
{
    static char result[256];

    strcpy (result, after (str, target1));
    strcpy (result, before (result, target2));

    return result;
}


/* Works like the C strtok() function except that it can handle */
/* tokens enclosed in double quotes */
char *parse_string (char *str)
{
    static char result[256];
    static char *p;
    char QUOTE = '\"';
    int  index;

    strcpy (result, "");
    index = 0;

    if (str != NULL)
	p = str;

    /* Find the start of the next token */
    while (isspace (*p))
	p++;

    if (*p == QUOTE) {
	p++;

	while (*p != '\0' && *p != QUOTE)
	    result[index++] = *p++;

	if (*p == QUOTE)
	    p++;
    }
    else {
	while (*p != '\0' && !isspace(*p))
	    result[index++] = *p++;
    }

    result[index] = '\0';

    return result;
}


/* Convert character 'c' to upper case */
char upcase (char c)
{
    if (c >= 'a' && c <= 'z')
	c = c - 'a' + 'A';

    return c;
}


float colour_intens (Colour *colour)
{
    return sqrt (colour->red   * colour->red +
		 colour->green * colour->green +
		 colour->blue  * colour->blue);
}


void parse_file()
{
    Chunk chunk;

    start_chunk(&chunk);

    if (chunk.tag == 0x4D4D)
	parse_3ds (&chunk);
    else
	abortmsg ("Error: Input file is not .3DS format", 1);

    end_chunk (&chunk);
}


void parse_3ds (Chunk *mainchunk)
{
    Chunk chunk;

    do  {
	start_chunk (&chunk);
	if (feof(in)) {
		fprintf(stderr, "%s: unexpected EOF\n", progname);
		break;
	}
	if (chunk.end <= mainchunk->end) {
	    switch (chunk.tag) {
		case 0x3D3D: parse_mdata (&chunk);
			     break;
	    }
	}

	end_chunk (&chunk);
    } while (chunk.end <= mainchunk->end);
}


void parse_mdata (Chunk *mainchunk)
{
    Chunk chunk;
    Colour bgnd_colour;

    do  {
	start_chunk (&chunk);

	if (chunk.end <= mainchunk->end) {
	    switch (chunk.tag) {
		case 0x2100: parse_colour (&global_amb);
			     break;
		case 0x1200: parse_colour (&bgnd_colour);
			     break;
		case 0x1201: write_bgsolid (out, bgnd_colour);
			     break;
		case 0x2200: parse_fog (&chunk);
			     break;
		case 0x2210: parse_fog_bgnd();
			     break;
		case 0x2201: write_fog (out, fog_colour, fog_distance);
			     break;
		case 0xAFFF: parse_mat_entry (&chunk);
			     break;
		case 0x4000: parse_named_object (&chunk);
			     break;
	    }
	}

	end_chunk (&chunk);
    } while (chunk.end <= mainchunk->end);
}


void parse_fog (Chunk *mainchunk)
{
    Chunk chunk;

    (void)read_float();
    (void)read_float();
    fog_distance = read_float();
    (void)read_float();

    parse_colour (&fog_colour);

    do  {
	start_chunk (&chunk);

	if (chunk.end <= mainchunk->end) {
	    switch (chunk.tag) {
		case 0x2210: parse_fog_bgnd();
			     break;
	    }
	}

	end_chunk (&chunk);
    } while (chunk.end <= mainchunk->end);
}


void parse_fog_bgnd()
{

}


void parse_mat_entry (Chunk *mainchunk)
{
    Chunk chunk;
    MatProp *mprop;

    mprop = create_mprop();

    do  {
	start_chunk (&chunk);

	if (chunk.end <= mainchunk->end) {
	    switch (chunk.tag) {
		case 0xA000: strcpy (mprop->name, read_string());
			     cleanup_name (mprop->name);
			     break;

		case 0xA010: parse_colour (&mprop->ambient);
			     break;

		case 0xA020: parse_colour (&mprop->diffuse);
			     break;

		case 0xA030: parse_colour (&mprop->specular);
			     break;

		case 0xA040: mprop->shininess = 100.0*parse_percentage();
			     break;

		case 0xA050: mprop->transparency = parse_percentage();
			     break;

		case 0xA080: mprop->self_illum = TRUE;
			     break;

		case 0xA081: mprop->two_side = TRUE;
			     break;

		case 0xA220: mprop->reflection = parse_percentage();
			     (void)parse_mapname (&chunk);
			     break;

		case 0xA310: if (mprop->reflection == 0.0)
				 mprop->reflection = 1.0;
			     break;

		case 0xA200: mprop->tex_strength = parse_percentage();
			     strcpy (mprop->tex_map, parse_mapname (&chunk));
			     break;

		case 0xA230: mprop->bump_strength = parse_percentage();
			     strcpy (mprop->bump_map, parse_mapname (&chunk));
			     break;
	    }
	}

	end_chunk (&chunk);
    } while (chunk.end <= mainchunk->end);

    LIST_INSERT (mprop_list, mprop);
}


char *parse_mapname (Chunk *mainchunk)
{
    static char name[80] = "";
    Chunk chunk;

    do  {
	start_chunk (&chunk);

	if (chunk.end <= mainchunk->end) {
	    switch (chunk.tag) {
		case 0xA300: strcpy (name, read_string());
			     break;
	    }
	}

	end_chunk (&chunk);
    } while (chunk.end <= mainchunk->end);

    return name;
}


void parse_named_object (Chunk *mainchunk)
{
    Chunk chunk;

    strcpy (obj_name, read_string());
    cleanup_name (obj_name);

    printf ("Working on: %s\n", obj_name);

    mesh = NULL;

    do  {
	start_chunk (&chunk);

	if (chunk.end <= mainchunk->end) {
	    switch (chunk.tag) {
		case 0x4100: parse_n_tri_object (&chunk);
			     break;
		case 0x4600: parse_n_direct_light (&chunk);
			     break;
		case 0x4700: parse_n_camera();
			     break;
		case 0x4010: if (mesh != NULL) mesh->hidden = TRUE;
			     break;
		case 0x4012: if (mesh != NULL) mesh->shadow = FALSE;
			     break;
	    }
	}

	end_chunk (&chunk);
    } while (chunk.end <= mainchunk->end);

    if (mesh != NULL) {
	update_limits (mesh);

	if (frame >= 0)
	    LIST_INSERT (mesh_list, mesh);
	else {
	    write_mesh (out, mesh);

	    free_mesh_data (mesh);
	    free (mesh);
	}
    }
}


void parse_n_tri_object (Chunk *mainchunk)
{
    Chunk chunk;

    mesh = create_mesh (obj_name, 0, 0);

    do  {
	start_chunk (&chunk);

	if (chunk.end <= mainchunk->end) {
	    switch (chunk.tag) {
		case 0x4110: parse_point_array();
			     break;
		case 0x4120: parse_face_array (&chunk);
			     break;
		case 0x4160: parse_mesh_matrix();
			     break;
	    }
	}

	end_chunk (&chunk);
    } while (chunk.end <= mainchunk->end);
}


void parse_point_array()
{
    int i;

    mesh->vertices = read_word();
    mesh->vertex = malloc (mesh->vertices * sizeof(*(mesh->vertex)));
    if (mesh->vertex == NULL)
	abortmsg ("Out of memory allocating mesh", 1);

    for (i = 0; i < mesh->vertices; i++)
	read_point (mesh->vertex[i]);
}


void parse_face_array (Chunk *mainchunk)
{
    Chunk chunk;
    int i;

    mesh->faces = read_word();
    mesh->face = malloc (mesh->faces * sizeof(*(mesh->face)));
    if (mesh->face == NULL)
	abortmsg ("Out of memory allocating mesh", 1);

    mesh->mtl = malloc (mesh->faces * sizeof(*(mesh->mtl)));
    if (mesh->mtl == NULL)
	abortmsg ("Out of memory allocating mesh", 1);

    for (i = 0; i < mesh->faces; i++) {
	mesh->face[i].a = read_word();
	mesh->face[i].b = read_word();
	mesh->face[i].c = read_word();
	(void)read_word();

	mesh->mtl[i] = NULL;
    }

    do  {
	start_chunk (&chunk);

	if (chunk.end <= mainchunk->end) {
	    switch (chunk.tag) {
		case 0x4130: parse_msh_mat_group();
			     break;
		case 0x4150: parse_smooth_group();
			     break;
	    }
	}

	end_chunk (&chunk);
    } while (chunk.end <= mainchunk->end);

    for (i = 0; i < mesh->faces; i++) {
	if (mesh->mtl[i] == NULL)
	    mesh->mtl[i] = update_materials ("Default", 0);
    }
}


void parse_msh_mat_group()
{
    Material *new_mtl;
    char mtlname[80];
    int  mtlcnt;
    int  i, face;

    strcpy (mtlname, read_string());
    cleanup_name (mtlname);

    new_mtl = update_materials (mtlname, 0);

    mtlcnt = read_word();

    for (i = 0; i < mtlcnt; i++) {
	face = read_word();
	mesh->mtl[face] = new_mtl;
    }
}


void parse_smooth_group()
{

}


void parse_mesh_matrix()
{
    int i, j;

    if (mesh != NULL) {
	for (i = 0; i < 4; i++) {
	    for (j = 0; j < 3; j++)
		mesh->matrix[i][j] = read_float();
	}

	mat_inv (mesh->invmatrix, mesh->matrix);
    }
}


void parse_n_direct_light (Chunk *mainchunk)
{
    Chunk chunk;
    Spotlight *s;
    OmniLight *o;
    int light_off = FALSE;
    int spot_flag = FALSE;

    read_point (pos);
    parse_colour (&col);

    do  {
	start_chunk (&chunk);

	if (chunk.end <= mainchunk->end) {
	    switch (chunk.tag) {
		case 0x4620: light_off = TRUE;
			     break;
		case 0x4610: parse_dl_spotlight();
			     spot_flag = TRUE;
			     break;
	    }
	}

	end_chunk (&chunk);
    } while (chunk.end <= mainchunk->end);

    if (light_off)
	return;

    if (!spot_flag) {
	if (frame >= 0) {
	    o = LIST_FIND (omni_list, obj_name);

	    if (o != NULL) {
		pos[X] = o->pos[X];
		pos[Y] = o->pos[Y];
		pos[Z] = o->pos[Z];
		col    = o->col;
	    }
	}

	write_light (out, obj_name, pos, col);
    }
    else {
	if (frame >= 0) {
	    s = LIST_FIND (spot_list, obj_name);

	    if (s != NULL) {
		pos[X]    = s->pos[X];
		pos[Y]    = s->pos[Y];
		pos[Z]    = s->pos[Z];
		target[X] = s->target[X];
		target[Y] = s->target[Y];
		target[Z] = s->target[Z];
		col       = s->col;
		hotspot   = s->hotspot;
		falloff   = s->falloff;
	    }
	}

	if (falloff <= 0.0)
	    falloff = 180.0;

	if (hotspot <= 0.0)
	    hotspot = 0.7*falloff;

	write_spot (out, obj_name, pos, target, col, hotspot, falloff);
    }
}


void parse_dl_spotlight()
{
    read_point (target);

    hotspot = read_float();
    falloff = read_float();
}


void parse_n_camera()
{
    float  bank;
    float  lens;

    read_point (pos);
    read_point (target);
    bank = read_float();
    lens = read_float();

    if (frame >= 0 && cam_list != NULL) {
	pos[X]    = cam_list->pos[X];
	pos[Y]    = cam_list->pos[Y];
	pos[Z]    = cam_list->pos[Z];
	target[X] = cam_list->target[X];
	target[Y] = cam_list->target[Y];
	target[Z] = cam_list->target[Z];
	lens      = cam_list->lens;
	bank      = cam_list->bank;
    }

    write_camera (out, obj_name, pos, target, lens, bank);
}


void parse_colour (Colour *colour)
{
    Chunk chunk;
    Colour_24 colour_24;

    start_chunk (&chunk);

    switch (chunk.tag) {
	case 0x0010: parse_colour_f (colour);
		     break;

	case 0x0011: parse_colour_24 (&colour_24);
		     colour->red   = colour_24.red/255.0;
		     colour->green = colour_24.green/255.0;
		     colour->blue  = colour_24.blue/255.0;
		     break;

	default:     abortmsg ("Error parsing colour", 1);
    }

    end_chunk (&chunk);
}


void parse_colour_f (Colour *colour)
{
    colour->red   = read_float();
    colour->green = read_float();
    colour->blue  = read_float();
}


void parse_colour_24 (Colour_24 *colour)
{
    colour->red   = read_byte();
    colour->green = read_byte();
    colour->blue  = read_byte();
}


float parse_percentage()
{
    Chunk chunk;
    float percent = 0.0;

    start_chunk (&chunk);

    switch (chunk.tag) {
	case 0x0030: percent = parse_int_percentage()/100.0;
		     break;

	case 0x0031: percent = parse_float_percentage();
		     break;

	default:     printf ("WARNING: Error parsing percentage");
    }

    end_chunk (&chunk);

    return percent;
}


short parse_int_percentage()
{
    word percent = read_word();

    return percent;
}


float parse_float_percentage()
{
    float percent = read_float();

    return percent;
}


void start_chunk (Chunk *chunk)
{
    chunk->start  = ftell(in);
    chunk->tag    = read_word();
    chunk->length = read_dword();
    if (chunk->length < sizeof(word)+sizeof(dword))
	chunk->length = sizeof(word) + sizeof(dword);
    chunk->end    = chunk->start + chunk->length;
}


void end_chunk (Chunk *chunk)
{
    fseek (in, chunk->end, 0);
}


byte read_byte()
{
    byte data;

    data = fgetc (in);

    return data;
}


word read_word()
{
    word data;

    data = fgetc (in);
    data |= fgetc (in) << 8;

    return data;
}


dword read_dword()
{
    dword data;

    data = read_word();
    data |= read_word() << 16;

    return data;
}


float read_float()
{
    union { dword i; char c[8]; } u;
    dword data;

    data = read_dword();

    if (sizeof(dword) == sizeof(float))
	return *(float *)&data;

    u.i = 1;
    if (u.c[0] == 0)
	return *(float *)&data;	/* assume big-endian */

    if (sizeof(dword) != 2*sizeof(float)) {
    	fputs("Unsupported word length\n", stderr);
	exit(1);
    }
    u.i = data;
    return *(float *)&u.c[4];
}


void read_point (Vector v)
{
    v[X] = read_float();
    v[Y] = read_float();
    v[Z] = read_float();
}


char *read_string()
{
    static char string[80];
    int i;

    for (i = 0; i < 80; i++) {
	string[i] = read_byte();

	if (string[i] == '\0')
	    break;
    }

    return string;
}


float findfov (float lens)
{
    static float lens_table[13] =
		 { 15.0, 17.0, 24.0, 35.0, 50.0, 85.0, 100.0, 135.0, 200.0,
		   500.0, 625.0, 800.0, 1000.0 };
    static float fov_table[13] =
		 { 115.0, 102.0, 84.0, 63.0, 46.0, 28.0, 24.0, 18.0,
		   12.0, 5.0, 4.0, 3.125, 2.5 };

    float fov, f1, f2, l1, l2;
    int   i;

    if (lens < 15.0)
	lens = 15.0;
    else if (lens > 1000.0)
	lens = 1000.0;

    for (i = 0; i < 13; i++)
	if (lens < lens_table[i])
	    break;

    if (i == 13)
	i = 12;
    else if (i == 0)
	i = 1;

    f1 = fov_table[i-1];
    f2 = fov_table[i];
    l1 = lens_table[i-1];
    l2 = lens_table[i];

    fov = f1 + (lens - l1) * (f2 - f1) / (l2 - l1);

    return fov;
}


