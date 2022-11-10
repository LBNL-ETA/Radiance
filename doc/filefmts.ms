.EQ
delim @@
.EN
.SH 1
\f[BI]Radiance\f[B] File Formats
.pdfhref O 1 "Radiance File Formats"
.pdfhref M "radiance-file-formats"
.LP
This chapter discusses the standard file formats specific to
\f[I]Radiance\f[R], and gives their internal structure, with pointers to
routines for reading and writing them.
The following file formats (listed with their conventional suffixes) are
covered:
.IP "\f[B]Scene Description (.rad suffix)\f[R]" 3
This is the main input file type, describing materials and geometry for
the rendering programs, and must be compiled into an octree by oconv
prior to ray-tracing.
It is an ASCII text format, and is often translated from a CAD
description, but may be created or edited by a text editor as well.
.RS 3
.RE
.IP "\f[B]Function File (.cal suffix)\f[R]" 3
Also a text format, these files describe mathematical patterns,
textures, and surface shapes.
In the case of patterns and textures, the functions serve as input
directly to the rendering programs.
In the case of surfaces, the functions serve as input to one of the
generator programs, \f[B]gensurf\f[R], \f[B]genrev\f[R] or
\f[B]genworm\f[R].
Additionally, \f[B]pcomb\f[R] may be used to perform math on
\f[I]Radiance\f[R] pictures and the \f[B]rcalc\f[R] utility may be used
in creative ways to manipulate data for scene generation and data
analysis.
.RS 3
.RE
.IP "\f[B]Data File (.dat suffix)\f[R]" 3
Another ASCII format, data files are used directly by the rendering
programs to interpolate values for luminaire photometry, among other
things.
.RS 3
.RE
.IP "\f[B]Font File (.fnt suffix)\f[R]" 3
A simple, polygonal font representation for rendering text patterns.
This ASCII format describes each character \[lq]glyph\[rq] as a sequence
of vertices in rectangular, integer coordinates ranging from 0 to 255.
.RS 3
.RE
.IP "\f[B]Octree (.oct suffix)\f[R]" 3
A binary data structure computed from one or more scene description
files by the oconv program.
It may contain frozen scene data in binary form, or merely references to
the original scene files.
.RS 3
.RE
.IP "\f[B]Picture (.hdr suffix)\f[R]" 3
A binary image file containing calibrated, real radiance values at each
pixel.
\f[I]Radiance\f[R] pictures may be displayed, analyzed, and converted to
other image formats.
.RS 3
.RE
.IP "\f[B]Z-buffer (.zbf suffix)\f[R]" 3
A binary file with the distances to each pixel in a corresponding
picture.
.RS 3
.RE
.IP "\f[B]Ambient File (.amb suffix)\f[R]" 3
A binary file used to store diffuse interreflection values, which are
shared between cooperating rendering processes running sequentially or
in parallel.
Since these values are view-independent, sharing this information across
multiple runs is highly economical.
.RS 3
.RE
.LP
We will discuss each of these formats in turn, giving examples and
pointers to routines in the source code for reading and writing them,
and the programs that use them.
In general, the ASCII text formats have no standard routines for writing
them, since they generally originate outside of \f[I]Radiance\f[R] or
are created with simple \f[I]printf(3)\f[R] statements.
Most binary file formats are machine and system independent, meaning
they can be moved safely from one place to another and
\f[I]Radiance\f[R] will still understand them (provided no unintentional
character translation takes place along the way)\**.
.FS
The single exception to this rule is the Z-buffer file, whose contents
are dictated by the floating point representation and byte order of the
originating machine.
This choice was made for economic reasons, and is rarely a problem.
.FE
Most binary files also include a standard ASCII header at their
beginning that may be read by the \f[B]getinfo\f[R] program.
This offers a convenient method for identifying the file contents when
the file name is ambiguous.
.SH 2
Scene Description Format (.rad suffix)
.pdfhref O 2 "Scene Description Format \(.rad suffix\)"
.pdfhref M "scene-description-format-.rad-suffix"
.LP
The semantics of the \f[I]Radiance\f[R] scene description format are
covered in the Reference Manual.
We will therefore focus on the file syntax and structure, which are
simple and straightforward.
In fact, some would say that the \f[I]Radiance\f[R] scene description
format is brain-dead, in the sense that it offers few language amenities
and requires the awkward counting of string and real arguments (not to
mention those non-existent integer arguments).
We have little to offer in its defense.
.PP
The truth is, the scene format was designed to grow with
\f[I]Radiance\f[R], and we wanted to keep it as simple as possible so as
to encourage others to write translators to and from it.
Specifically, we wanted to be able to read files using the
\f[I]scanf(3)\f[R] library function and write files using
\f[I]printf(3)\f[R].
Furthermore, we wanted everyone\[cq]s parsers to be stable over time,
which meant no primitive-specific syntax.
We also decided that a flat file structure was most practical, since
hierarchies are typically lost on the first translation, and sufficient
structure could be provided by the file system itself.
Since we did not intend text editing to be the primary input method, we
felt the effects of these programming decisions on the human readability
and writability of the format were less important.
.PP
Even so, the format is relatively easy to read and write once you get
used to it, and with the \f[I]Radiance\f[R] generator programs and
in-line command expansion, the text editor becomes a powerful modeling
tool in the hands of an experienced user.
Together with the need for editing material descriptions, our assumption
that users would rarely edit these files turned out to be mistaken.
Consequently, it is a good idea for all users to familiarize themselves
with the scene description format, awkward or not.
.SH 3
Basic File Structure
.pdfhref O 3 "Basic File Structure"
.pdfhref M "basic-file-structure"
.LP
There are four statement types in a scene description file: comments,
commands, primitives and aliases.
These may be interspersed in the file, and the only structural
requirement is that modifiers precede the primitives they modify.
.SH 4
Comments
.pdfhref O 4 "Comments"
.pdfhref M "comments"
.LP
The simplest statement type is a comment statement begins with a pound
sign (`#') and continues to the end of line:
.IP
.nf
\f[C]
# This is a comment.
\f[]
.fi
.SH 4
Commands
.pdfhref O 4 "Commands"
.pdfhref M "commands"
.LP
An in-line command, which begins with an exclamation mark (`!') and
continues to the end of line:
.IP
.nf
\f[C]
!xform -n chair1 -t 10 5 8 chair.rad
\f[]
.fi
.LP
The command is executed during file parsing, and its output is read as
more input.
Long commands may be continued on multiple lines by escaping the newline
character with a backslash (`\[rs]'):
.IP
.nf
\f[C]
!gensurf marble sink \[aq]15.5+x(theta(s),phi(t))\[aq] \[rs]
     \[aq]10.5+y(theta(s),phi(t))\[aq] \[rs]
     \[aq]30.75+z(theta(s),phi(t))\[aq] \[rs]
     8 29 -f basin.cal -s
\f[]
.fi
.LP
Since the command is executed by the shell, pipes and other facilities
are available as well.
The following command creates a counter with a precisely cut hole for
the sink basin just given:
.IP
.nf
\f[C]
!( echo marble polygon sink_top 0 0 108 31 \[rs]
     10.5 30.75 31 22 30.75 0 22 30.75 0 0 \[rs]
     30.75 31 0 30.75 31 10.5 30.75 ; \[rs]
     cnt 30 | rcalc \[rs]
          -e \[aq]$1=15.5+x(theta(0),phi(1-$1/29))\[aq] \[rs]
          -e \[aq]$2=10.5+y(theta(0),phi(1-$1/29))\[aq] \[rs]
          -e \[aq]$3=30.75\[aq] -f basin.cal )
\f[]
.fi
.LP
Note in the above example that two commands are executed in sequence.
The first creates the counter perimeter, and the second cuts the hole.
The two commands are enclosed in parentheses, so if a final
transformation is added by \f[B]xform\f[R] with the \f[B]-c\f[R] option,
it will be applied to both commands, not just the last one.
.SH 4
Primitives
.pdfhref O 4 "Primitives"
.pdfhref M "primitives"
.LP
A primitive can be thought of as an indivisible unit of scene
information.
It can be a surface, material, pattern, texture or mixture.
The basic structure of a primitive is as follows:
.IP
.nf
\f[C]
modifier type identifier
n  S1  S2  S3  ..Sn
0
m  R1  R2  R3  ..Rm
\f[]
.fi
.LP
The \f[C]modifier\f[R] is the \f[C]identifier\f[R] of a previously
defined primitive, or \f[C]void\f[R] if no modifier is appropriate.
The type is one of the supported \f[I]Radiance\f[R] primitive keywords,
such as \f[I]polygon\f[R] or \f[I]plastic\f[R].
Following the \f[C]modifier\f[R], type and identifier are the string
arguments, preceded by the number of string arguments and separated by
white space.
If there are no string arguments, then 0 should be given for
\f[C]n\f[R].
The string arguments are followed by the integer arguments in the same
fashion.
(Since there are no \f[I]Radiance\f[R] primitives currently using
integer arguments, the count is always 0.)
Finally, the number of real arguments is given, followed by the real
arguments.
.PP
The location of the primitive in the scene description has no
importance, except that its modifier refers to the most recently defined
primitive with that identifier.
If no such modifier was defined, an error results.
In fact, \[lq]undefined modifier\[rq] is the most frequently reported
error when parsing an invalid scene description, since any random bit of
junk encountered where a statement is expected will be interpreted as a
modifier.
One final note about modifiers \[em] since surfaces never modify
anything, their identifiers are neither stored nor referenced in the
parser\[cq]s modifier list, and serve only for debugging purposes during
editing and rendering.
.PP
Within a primitive, white space serves only to separate words, and
multiple spaces, tabs, form feeds, returns, and newlines are all
considered as one separator.
Consequently, it is not possible for a string argument to contain any
white space, which is OK because no \f[I]Radiance\f[R] primitive needs
this.
.SH 4
Aliases
.pdfhref O 4 "Aliases"
.pdfhref M "aliases"
.LP
An alias simply associates a new modifier and identifier with a
previously defined primitive.
The syntax is as follows:
.IP
.nf
\f[C]
modifier alias new_identifier old_identifier
\f[]
.fi
.LP
The \f[C]old_identifier\f[R] should be associated with some modifier
primitive (i.e., non-surface) given earlier.
The \f[C]modifier\f[R], if different from the original, will be used
instead in later applications of \f[C]new_identifier\f[R].
.PP
Aliases are most often used to give new names to previously defined
materials.
They may also be used to associate different patterns or textures with
the same material.
.SH 3
Scene Hierarchy
.pdfhref O 3 "Scene Hierarchy"
.pdfhref M "scene-hierarchy"
.LP
Hierarchical scene descriptions are achieved through expansion of
in-line \f[B]xform\f[R] commands.
The \f[B]xform\f[R] command is used to read and place other
\f[I]Radiance\f[R] scene description files in the calling file, and
these other descriptions may in turn read others, and so on down the
tree.
No check is made to assure that none of the calling files is called
again, even by itself.
If this happens, commands open commands until the system runs out of
processes, which is a very nasty business and to be avoided.
.SH 3
\f[BI]Radiance\f[B] Programs
.pdfhref O 3 "Radiance Programs"
.pdfhref M "radiance-programs"
.LP
The following table shows programs in the main \f[I]Radiance\f[R]
distribution that read and write scene description files.
Additionally, there are other translators that write scene files, which
are available separately as free contributions or as part of other (CAD)
programs.
.PP
.na
.TS
delim(@@) tab(	);
l l l l.
T{
Program
T}	T{
Read
T}	T{
Write
T}	T{
Function
T}
_
T{
\f[B]arch2rad\f[R]
T}	T{
T}	T{
X
T}	T{
Convert Architrion text file to \f[I]Radiance\f[R]
T}
T{
\f[B]genblinds\f[R]
T}	T{
T}	T{
X
T}	T{
Generate curved venetian blinds
T}
T{
\f[B]genbox\f[R]
T}	T{
T}	T{
X
T}	T{
Generate parallelepiped
T}
T{
\f[B]genclock\f[R]
T}	T{
T}	T{
X
T}	T{
Generate analog clock
T}
T{
\f[B]genprism\f[R]
T}	T{
T}	T{
X
T}	T{
Generate extruded polygon
T}
T{
\f[B]genrev\f[R]
T}	T{
T}	T{
X
T}	T{
Generate surface of revolution
T}
T{
\f[B]gensky\f[R]
T}	T{
T}	T{
X
T}	T{
Generate CIE sky distribution
T}
T{
\f[B]gensurf\f[R]
T}	T{
T}	T{
X
T}	T{
Generate arbitrary surface patch
T}
T{
\f[B]genworm\f[R]
T}	T{
T}	T{
X
T}	T{
Generate varying diameter curved path
T}
T{
\f[B]ies2rad\f[R]
T}	T{
T}	T{
X
T}	T{
Convert IES luminaire file to \f[I]Radiance\f[R]
T}
T{
\f[B]mgf2rad\f[R]
T}	T{
T}	T{
X
T}	T{
Convert MGF file to \f[I]Radiance\f[R]
T}
T{
\f[B]mkillum\f[R]
T}	T{
X
T}	T{
X
T}	T{
Compute \f[I]illum\f[R] secondary sources
T}
T{
\f[B]nff2rad\f[R]
T}	T{
T}	T{
X
T}	T{
Convert NFF file to \f[I]Radiance\f[R]
T}
T{
\f[B]objline\f[R]
T}	T{
T}	T{
X
T}	T{
Generate line drawing of \f[I]Radiance\f[R] file
T}
T{
\f[B]objview\f[R]
T}	T{
T}	T{
X
T}	T{
Quick view of \f[I]Radiance\f[R] object
T}
T{
\f[B]oconv\f[R]
T}	T{
T}	T{
X
T}	T{
Compile \f[I]Radiance\f[R] scene description
T}
T{
\f[B]obj2rad\f[R]
T}	T{
T}	T{
X
T}	T{
Convert Wavefront .OBJ file to \f[I]Radiance\f[R]
T}
T{
\f[B]rad\f[R]
T}	T{
X
T}	T{
T}	T{
Render \f[I]Radiance\f[R] scene
T}
T{
\f[B]rad2mgf\f[R]
T}	T{
X
T}	T{
T}	T{
Convert \f[I]Radiance\f[R] file to MGF
T}
T{
\f[B]raddepend\f[R]
T}	T{
X
T}	T{
T}	T{
Determine scene file dependencies
T}
T{
\f[B]replmarks\f[R]
T}	T{
X
T}	T{
X
T}	T{
Replace triangular markers with objects
T}
T{
\f[B]rpict\f[R]
T}	T{
X
T}	T{
T}	T{
Batch rendering program
T}
T{
\f[B]rtrace\f[R]
T}	T{
X
T}	T{
T}	T{
Customizable ray-tracer
T}
T{
\f[B]rview\f[R]
T}	T{
X
T}	T{
T}	T{
Interactive renderer
T}
T{
\f[B]thf2rad\f[R]
T}	T{
T}	T{
X
T}	T{
Convert GDS things file to \f[I]Radiance\f[R]
T}
T{
\f[B]tmesh2rad\f[R]
T}	T{
T}	T{
X
T}	T{
Convert triangle mesh file to \f[I]Radiance\f[R]
T}
T{
\f[B]xform\f[R]
T}	T{
X
T}	T{
X
T}	T{
Transform Radiance objects
T}
.TE
.ad
.LP
\f[B]Table 1.\f[R] Radiance programs that read and write scene
descriptions.
.SH 3
\f[BI]Radiance\f[B] C Library
.pdfhref O 3 "Radiance C Library"
.pdfhref M "radiance-c-library"
.LP
The principal library function for reading scene description files is
\f[C]readobj(inpspec)\f[R], defined in \f[C]src/common/readobj.c\f[R].
This routine takes the name of a file, or command beginning with `!', or
\f[C]NULL\f[R] if standard input is to be read, and loads the
\f[I]Radiance\f[R] data structures defined in
\f[C]src/common/object.h\f[R].
If loading \f[I]Radiance\f[R] data structures is not the action desired,
then a more custom approach is necessary, such as that used in
\f[C]src/gen/xform.c\f[R].
If using \f[I]Radiance\f[R] data structures is acceptable, but the data
need not remain resident in memory, then follow the lead in
\f[C]src/ot/getbbox.c\f[R] and use \f[C]src/ot/readobj2.c\f[R] instead.
In any case, the list of defined primitive types in
\f[C]src/common/otypes.h\f[R] is crucial.
.SH 2
Function File Format (.cal suffix)
.pdfhref O 2 "Function File Format \(.cal suffix\)"
.pdfhref M "function-file-format-.cal-suffix"
.LP
Function files are used throughout \f[I]Radiance\f[R] to specify
mathematical formulas and relations for procedural textures, patterns
and surfaces.
They are also used by filter programs such as rcalc to manipulate data,
and pcomb to manipulate pictures.
.PP
Function file syntax is simple and should be familiar to most
programmers, as it is based on fairly standard algebraic expressions.
Here is an example, which corresponds to the in-line commands given in
the previous section:
.IP
.nf
\f[C]
{
basin.cal - calculate coordinates for basin sink.
}

theta(s) = PI*(0.5+0.5*s);
phi(t) = 2*PI*t;

R(th,p) = 5 + ( 3.25*cos(p)\[ha]2 +
                1.75*sin(p)\[ha]2 ) * sin(th)\[ha]2;

x(th,p) = R(th,p)*sin(th)*cos(p);
y(th,p) = R(th,p)*sin(th)*sin(p);
z(th,p) = R(th,p)*cos(th);
\f[]
.fi
.LP
In contrast to the usual semantics in programs where each statement
corresponds to an evaluation, statements in function files correspond to
\f[I]definitions\f[R].
Once a function or variable has been defined, it may be used in later
definitions, along with predefined functions such as \f[C]sin(x)\f[R]
and \f[C]cos(x)\f[R] and constants such as PI \**.
.FS
TBD - There once was a footnote here
.FE
(All math functions use standard C conventions, hence trigonometry is
done in radians rather than degrees.)
.PP
Evaluation order (operator precedence) follows standard rules of
algebra.
Exponentiation is performed first \f[C](x\[ha]y)\f[R], followed by
multiplication and division \f[C](x*y, x/y)\f[R], then addition and
subtraction \f[C](x+y, x-y)\f[R].
Unary minus is most tightly associated \f[C](-x)\f[R], and parentheses
override operator precedence in the usual way.
Semicolons separate statements, and white space is generally ignored.
Comments are enclosed by curly braces, which may be nested.
.PP
The above file does not actually \f[I]do\f[R] anything, it merely
defines functions that are useful by a program that does.
Taking our \f[B]gensurf\f[R] example from the previous section:
.IP
.nf
\f[C]
!gensurf marble sink \[aq]15.5+x(theta(s),phi(t))\[aq] \[rs]
     \[aq]10.5+y(theta(s),phi(t))\[aq] \[rs]
     \[aq]30.75+z(theta(s),phi(t))\[aq] \[rs]
     8 29 -f basin.cal -s
\f[]
.fi
.LP
The \f[B]-f\f[R] option loads in our file, which is then used to
evaluate expressions such as \f[C]15.5+x(theta(s),phi(t))\f[R] for
specific values of \f[C]s\f[R] and \f[C]t\f[R].
These variables range from 0 to 1 over the surface patch in increments
of @1 / 8@ and @1 / 29@, respectively.
(See the \f[B]gensurf\f[R] manual page for details.)
The entire expression for each evaluation could have been written in the
command line, but it is much more convenient to create a function file.
.SH 3
Language Features
.pdfhref O 3 "Language Features"
.pdfhref M "language-features"
.LP
Subtle features of the functional language provide much greater power
than first meets the eye.
One of these is the ability to define recursive functions.
The following example defines the factorial function (\f[I]n!\f[R]):
.IP
.nf
\f[C]
fact(n) : if(n-1.5, n*fact(n-1), 1);
\f[]
.fi
.LP
This uses the library function \f[C]if(cond,e1,e0)\f[R], which returns
\f[C]e1\f[R] if cond is greater than zero, and \f[C]e0\f[R] otherwise.
Since only one of these expressions is evaluated, \f[C]fact(n)\f[R] will
call itself until \f[C]n\f[R] is less than 2, when the \f[C]if\f[R]
expression returns 1\**.
.FS
Note that we compare \f[C]n\f[R] to 1.5, so as to avoid any round-off
problems caused by floating point math.
Caution is advised because all expressions are evaluated as
double-precision real, and comparisons to zero are unreliable.
.FE
The colon (`:') is used in place of the usual equals assignment (`=')
because we want this function to have the constant attribute, which
means any later appearance in an expression of \f[C]fact(ce)\f[R] where
ce is also a constant expression will be replaced by its value.
This can be an important savings in cases where an expression or
subexpression is expensive to evaluate, and only needs to be computed
once.
All of the standard library functions have the constant attribute.
(See the following section for a complete list.)
.PP
Another handy language feature is the ability to pass function names as
arguments.
A simple example of this is the following function, which computes the
numerical derivative of a function given as its first argument:
.IP
.nf
\f[C]
FTINY : 1e-7;
d1(f,x) = (f(x+FTINY)-f(x-FTINY))/FTINY/2;
\f[]
.fi
.LP
Evaluating \f[C]d1(sin,1.1)\f[R] using this formula yields 0.4536, which
is fairly close to the true derivative, which is \f[C]cos(1.1)\f[R].
.PP
A third language feature, which is normally transparent to the user, is
the notion of \f[I]contexts\f[R].
Identifiers may be composed of parts, starting with a name and
continuing with one or more context names.
Each name is delimited by a back-quote (`\[ga]').
Names themselves begin with a letter and continue with any sequence of
letters, digits, underscores and decimal points.
The following are examples of valid identifiers:
.IP
.nf
\f[C]
v1, V.W, x_rand\[ga]local, \[ga]A_, Qa_5\[ga]
\f[]
.fi
.LP
If a context mark appears at the beginning of the identifier, then its
reference must be local.
If it appears at the end, then its reference must be global.
A local reference must be resolved in the local context, i.e., no
contexts above this one will be searched.
A global reference must correspond to the original context, ignoring any
local redefinitions.
.PP
The reason that contexts are normally transparent is that they are
controlled by the calling program \[em] there are no explicit language
features for establishing contexts.
A new context is established automatically for each function file loaded
by the rendering programs.
That way, it is safe to reuse variable names that have been used in
other files, and even in the main initialization file,
\f[C]rayinit.cal\f[R].
.PP
Although not strictly necessary, there are two good reasons to define
variables and functions before referencing them in a function file.
One is related to contexts.
If a previous definition of a variable name is given in an enclosing
context (e.g., \f[C]rayinit.cal\f[R]), then that reference will be used
rather than a later one in the current context, unless the reference is
made explicitly local by starting the identifier with a context mark.
The second reason for defining before referencing is constant
expressions.
If a variable or function has the constant attribute (i.e., defined with
`:' instead of `='), then a later subexpression referencing it can be
replaced by its evaluated result during compilation.
If the constant is not defined until after it is referenced, it remains
as a reference, which takes time to evaluate each time.
.PP
Other features of the language are truly transparent, but knowledge of
them can help one to write more efficient function files:
.IP \[bu] 3
Once a variable has been evaluated, the result is cached and it is not
reevaluated unless the client program changes an internal counter
(\f[C]eclock\f[R]), which indicates that something has changed.
This means that using variables to hold frequently used values will not
only simplify the function file, it will save time during evaluation.
.IP \[bu] 3
An argument passed in a function call is not evaluated until the
function calls for it specifically, and the result is also cached to
avoid redundant calculation.
The conditional evaluation feature is actually a requirement for
recursive functions to work, but caching is not.
Argument value caching means it is more efficient to pass an
expensive-to-compute expression than to have the function compute it
internally if it appears more than once in the function definition.
This is especially true for recursive functions with deep call trees.
.SH 3
Standard Definitions (Library)
.pdfhref O 3 "Standard Definitions \(Library\)"
.pdfhref M "standard-definitions-library"
.LP
The following are always defined:
.IP "\f[B]\f[CB]if(a, b, c)\f[B]\f[R]" 3
Conditional expression.
If a is positive, return b, else return c.
.RS 3
.RE
.IP "\f[B]\f[CB]select(N, a1, a2, ..)\f[B]\f[R]" 3
Return Nth argument.
If N is 0, then return the count of arguments excluding the first.
This provides basic array functionality.
.RS 3
.RE
.IP "\f[B]\f[CB]sqrt(x)\f[B]\f[R]" 3
Return square root of \f[C]x\f[R], where \f[C]x >= 0\f[R].
.RS 3
.RE
.IP "\f[B]\f[CB]sin(x), cos(x), tan(x), asin(x), acos(x), atan(x), atan2(y,x)\f[B]\f[R]" 3
Standard trigonometry functions.
.RS 3
.RE
.IP "\f[B]\f[CB]floor(x), ceil(x)\f[B]\f[R]" 3
Greatest lower bound and least upper bound (integer).
.RS 3
.RE
.IP "\f[B]\f[CB]exp(x), log(x), log10(x)\f[B]\f[R]" 3
Exponent and logarithm functions.
.RS 3
.RE
.IP "\f[B]\f[CB]rand(x)\f[B]\f[R]" 3
Return pseudo-random number in the range [0,1) for any argument x.
The same return value is guaranteed for the same argument.
.RS 3
.RE
.LP
The following are sometimes defined, depending on the program:
.IP "\f[B]PI\f[R]" 3
The ratio of a circle\[cq]s circumference to its diameter.
.RS 3
.RE
.IP "\f[B]\f[CB]erf(z), erfc(z)\f[B]\f[R]" 3
Error function and complimentary error function.
.RS 3
.RE
.IP "\f[B]\f[CB]j0(x), j1(x), jn(n,x), y0(x), y1(x), yn(n,x)\f[B]\f[R]" 3
Bessel functions.
.RS 3
.RE
.IP "\f[B]\f[CB]hermite(p0,p1,r0,r1,t)\f[B]\f[R]" 3
One-dimensional Hermite polynomial.
.RS 3
.RE
.LP
The rendering programs also define the following noise functions:
.IP "\f[B]\f[CB]noise3(x,y,z), noise3x(x,y,z), noise3y(x,y,z), noise3z(x,y,z)\f[B]\f[R]" 3
Perlin noise function and its gradient [Perlin85][Arvo91,p.396].
.RS 3
.RE
.IP "\f[B]\f[CB]fnoise3(x,y,z)\f[B]\f[R]" 3
Fractal noise function, ranging from -1 to 1.
.RS 3
.RE
.LP
Interaction with the renderer is achieved via special purpose variables
and functions whose values correspond to the current ray intersection
and the calling primitive.
Unlike the above functions, none of these have the constant attribute
since their values change from one ray to the next:
.IP "\f[B]\f[CB]Dx, Dy, Dz\f[B]\f[R]" 3
ray direction
.RS 3
.RE
.IP "\f[B]\f[CB]Nx, Ny, Nz\f[B]\f[R]" 3
surface normal
.RS 3
.RE
.LP
\f[C]Px, Py, Pz\f[R] :intersection point
.IP "\f[B]\f[CB]T\f[B]\f[R]" 3
distance from start
.RS 3
.RE
.IP "\f[B]\f[CB]Ts\f[B]\f[R]" 3
single ray (shadow) distance
.RS 3
.RE
.IP "\f[B]\f[CB]Rdot\f[B]\f[R]" 3
ray dot product
.RS 3
.RE
.IP "\f[B]\f[CB]S\f[B]\f[R]" 3
world scale
.RS 3
.RE
.IP "\f[B]\f[CB]Tx, Ty, Tz\f[B]\f[R]" 3
world origin
.RS 3
.RE
.IP "\f[B]\f[CB]Ix, Iy, Iz\f[B]\f[R]" 3
world i unit vector
.RS 3
.RE
.IP "\f[B]\f[CB]Jx, Jy, Jz\f[B]\f[R]" 3
world j unit vector
.RS 3
.RE
.IP "\f[B]\f[CB]Kx, Ky, Kz\f[B]\f[R]" 3
world k unit vector
.RS 3
.RE
.IP "\f[B]\f[CB]arg(n)\f[B]\f[R]" 3
real arguments, arg(0) is count
.RS 3
.RE
.LP
For BRDF primitives, the following variables are also available:
.IP "\f[B]\f[CB]NxP, NyP, NzP\f[B]\f[R]" 3
perturbed surface normal
.RS 3
.RE
.IP "\f[B]\f[CB]RdotP\f[B]\f[R]" 3
perturbed ray dot product
.RS 3
.RE
.IP "\f[B]\f[CB]CrP, CgP, CbP\f[B]\f[R]" 3
perturbed material color
.RS 3
.RE
.LP
For prism1 and prism2 primitives, the following are available:
.IP "\f[B]\f[CB]DxA, DyA, DzA\f[B]\f[R]" 3
direction to target light source
.RS 3
.RE
.LP
Other functions, variables and constants are defined as well in the file
\f[C]src/rt/rayinit.cal\f[R], which gets installed in the standard
\f[I]Radiance\f[R] library directory and can be modified or appended as
desired\**.
.FS
It is usually a good idea to store any such customized files in a
personal library location and set the \f[C]RAYPATH\f[R] environment
variable to search there first.
This way, it does not affect other users or get overwritten during the
next system installation.
.FE
.SH 3
\f[BI]Radiance\f[B] Programs
.pdfhref O 3 "Radiance Programs"
.pdfhref M "radiance-programs-1"
.LP
Table 2 shows \f[I]Radiance\f[R] programs that read and write function
files.
.PP
.na
.TS
delim(@@) tab(	);
l l l l.
T{
Program
T}	T{
Read
T}	T{
Write
T}	T{
Function
T}
_
T{
\f[B]calc\f[R]
T}	T{
X
T}	T{
X
T}	T{
Interactive calculator
T}
T{
\f[B]genrev\f[R]
T}	T{
X
T}	T{
T}	T{
Generate surface of revolution
T}
T{
\f[B]gensurf\f[R]
T}	T{
X
T}	T{
T}	T{
Generate arbitrary surface patch
T}
T{
\f[B]genworm\f[R]
T}	T{
X
T}	T{
T}	T{
Generate varying diameter curved path
T}
T{
\f[B]macbethcal\f[R]
T}	T{
T}	T{
X
T}	T{
Compute image color & contrast correction
T}
T{
\f[B]pcomb\f[R]
T}	T{
X
T}	T{
T}	T{
Perform arbitrary math on picture(s)
T}
T{
\f[B]rcalc\f[R]
T}	T{
X
T}	T{
T}	T{
Record stream calculator
T}
T{
\f[B]rpict\f[R]
T}	T{
X
T}	T{
T}	T{
Batch rendering program
T}
T{
\f[B]rtrace\f[R]
T}	T{
X
T}	T{
T}	T{
Customizable ray-tracer
T}
T{
\f[B]rview\f[R]
T}	T{
X
T}	T{
T}	T{
Interactive renderer
T}
T{
\f[B]tabfunc\f[R]
T}	T{
T}	T{
X
T}	T{
Create function file from tabulated data
T}
.TE
.ad
.LP
\f[B]Table 2.\f[R] Programs in the \f[I]Radiance\f[R] distribution that
read and write function files.
.PP
In addition, the program \f[B]ev\f[R] evaluates expressions given as
command line arguments, though it does not handle function files or
definitions.
There are also a number of 2-d plotting routines that use a slightly
modified statement syntax, called \f[B]bgraph\f[R], \f[B]dgraph\f[R],
\f[B]gcomp\f[R], and \f[B]igraph\f[R].
Additional utility programs are useful in combination with rcalc for
data analysis and scene generation.
The program \f[B]cnt\f[R] generates simple records to drive
\f[B]rcalc\f[R], and the \f[B]total\f[R] program is handy for adding up
results.
The \f[B]histo\f[R] program computes histograms needed for certain types
of statistical analysis.
The \f[B]lam\f[R] program concatenates columns from multiple input
files, and \f[B]neat\f[R] neatens up columns for better display.
.SH 3
\f[BI]Radiance\f[B] C Library
.pdfhref O 3 "Radiance C Library"
.pdfhref M "radiance-c-library-1"
.LP
The standard routines for loading and evaluating function files are
divided into three modules, \f[C]src/common/calexpr.c\f[R] for
expression parsing and evaluation, \f[C]src/common/caldefn.c\f[R] for
variable and function storage and lookup, and
\f[C]src/common/calfunc.c\f[R] for library function storage and function
evaluation.
There is a fourth module for writing out expressions called
\f[C]src/common/calprnt.c\f[R], which we will not discuss.
They all use the header file \f[C]src/common/calcomp.h\f[R], which
defines common data structures and evaluation macros.
Of these, the three most often used declarations for external routines
are:
.IP "\f[B]\f[CB]typedef struct epnode EPNODE;\f[B]\f[R]" 3
Expression parse tree node.
Some routines return pointers to this structure type, and the main
evaluation macro, \f[C]evalue(ep)\f[R], takes an argument of this type.
.RS 3
.RE
.IP "\f[B]\f[CB](double) evalue(ep);\f[B]\f[R]" 3
Evaluate an expression parse tree.
Uses node type table to access appropriate function depending on root
node type.
(E.g., an addition node calls \f[C]eadd(ep)\f[R].)
.RS 3
.RE
.IP "\f[B]\f[CB]extern unsigned long eclock;\f[B]\f[R]" 3
This is a counter used to determine when variable values need updating.
The initial value is 0, which tells the routines always to reevaluate
variables.
Once incremented to 1, variable evaluations are cached and not
recomputed until \f[C]eclock\f[R] is incremented again.
Usually, client programs increment \f[C]eclock\f[R] each time
definitions or internal states affecting returned values change.
This assures the quickest evaluation of correct results.
.RS 3
.RE
.LP
The usual approach to handling definitions is to compile them into the
central lookup table; variable and function references are later
evaluated by traversing the stored parse trees.
Syntax errors and undefined symbol errors during evaluation result in
calls to the user-definable routine \f[C]eputs(s)\f[R] to report the
error and \f[C]quit(status)\f[R] to exit the program.
Domain and range errors during evaluation set \f[C]errno\f[R], then call
the user-definable routine \f[C]wputs(s)\f[R] to report the error and
return zero as the expression result.
Following are standard routines provided for parsing from a file and
parsing from a string:
.IP "\f[B]\f[CB]EPNODE *eparse(char *expr);\f[B]\f[R]" 3
Compile the string expr into a parse tree for later evaluation with
evalue(ep).
.RS 3
.RE
.IP "\f[B]\f[CB]epfree(EPNODE *ep);\f[B]\f[R]" 3
Free memory associated with ep, including any variables referenced if
they are no longer defined.
.RS 3
.RE
.IP "\f[B]\f[CB]double eval(char *expr);\f[B]\f[R]" 3
Immediately parse, evaluate and free the given string expression.
.RS 3
.RE
.IP "\f[B]\f[CB]fcompile(char *fname);\f[B]\f[R]" 3
Compile definitions from the given file, or standard input if fname is
NULL.
.RS 3
.RE
.IP "\f[B]\f[CB]scompile(char *str, char *fn, int ln);\f[B]\f[R]" 3
Compile definitions from the string str, taken from file fn at line
number ln.
If no file is associated, fn can be NULL, and ln can be 0.
.RS 3
.RE
.LP
The following routines allow one to control the current context for
definition lookup and storage:
.IP "\f[B]\f[CB]char *setcontext(char *ctx);\f[B]\f[R]" 3
Set the current context to ctx.
If ctx is NULL, then simply return the current context.
An empty string sets the global context.
.RS 3
.RE
.IP "\f[B]\f[CB]char *pushcontext(char *name);\f[B]\f[R]" 3
Push another context onto the stack.
Return the new (full) context.
.RS 3
.RE
.IP "\f[B]\f[CB]char *popcontext();\f[B]\f[R]" 3
Pop the top context name off the stack.
Return the new (full) context.
.RS 3
.RE
.LP
The following functions permit the explicit setting of variable and
function values:
.IP "\f[B]\f[CB]varset(char *vname, int assign, double val);\f[B]\f[R]" 3
Set the specified variable to the given value, using a constant
assignment if assign is `:' or a regular one if it is `='.
This is always faster than compiling a string to do the same thing.
.RS 3
.RE
.IP "\f[B]\f[CB]funset(char *fname, int nargs, int assign, double (*fptr)(char *fn));\f[B]\f[R]" 3
Assign the specified library function, which takes a minimum of nargs
arguments.
The function will have the constant attribute if assign is `:', or not
if it is `='.
The only argument to the referenced function pointer is the function
name, which will equal fname.
(This string must therefore be declared static.)
This offers a convenient method to identify calls to an identical
function assigned multiple tasks.
Argument values are obtained with calls back to the argument(n) library
function.
.RS 3
.RE
.LP
The following functions are provided for evaluating a function or
variable in the current definition set:
.IP "\f[B]\f[CB]double varvalue(char *vname);\f[B]\f[R]" 3
Evaluate the given variable and return the result.
Since a hash lookup is necessary to resolve the reference, this is
slightly less efficient than evaluating a compiled expression via
evalue(ep), which uses soft links generated and maintained during
compilation.
.RS 3
.RE
.IP "\f[B]\f[CB]double funvalue(char *fn, int n, double a);\f[B]\f[R]" 3
Evaluate the function fn, passing n real arguments in the array a.
There is currently no mechanism for passing functions or function name
arguments from client programs.
.RS 3
.RE
.LP
These routines can be used to check the existence of a specific function
or variable:
.IP "\f[B]\f[CB]int vardefined(char *vname);\f[B]\f[R]" 3
Return non-zero if the named variable is defined.
(If the symbol is defined as a function, zero is returned.)
.RS 3
.RE
.IP "\f[B]\f[CB]int fundefined(char *fname);\f[B]\f[R]" 3
Return the number of required arguments for the named function if it is
defined, or zero if it is not defined.
(If the symbol is defined as a variable, zero is returned.)
.RS 3
.RE
.LP
These routines allow definitions to be cleared:
.IP "\f[B]\f[CB]dclear(char *dname);\f[B]\f[R]" 3
Clear the given variable or function, unless it is a constant
expression.
.RS 3
.RE
.IP "\f[B]\f[CB]dremove(char *dname);\f[B]\f[R]" 3
Clear the given variable or function, even if it is a constant
expression.
Library definitions cannot be removed, except by calling funset with a
\f[C]NULL\f[R] pointer for the function argument.
.RS 3
.RE
.IP "\f[B]\f[CB]dcleanup(int level);\f[B]\f[R]" 3
Clear definitions.
If level is 0, then just clear variable definitions.
If level is 2, then clear constants as well.
If the current context is local, then only local definitions will be
affected.
If global, all definitions in all contexts will be affected.
.RS 3
.RE
.LP
These routines may be used during library function evaluation:
.IP "\f[B]\f[CB]int nargum();\f[B]\f[R]" 3
Determine the number of arguments available in the current function
evaluation context.
.RS 3
.RE
.IP "\f[B]\f[CB]double argument(int n);\f[B]\f[R]" 3
Evaluate and return the nth argument.
.RS 3
.RE
.IP "\f[B]\f[CB]char *argfun(n);\f[B]\f[R]" 3
Get the name of the function passed as argument n.\ (Generates an error
if the nth argument is not a function.)
.RS 3
.RE
.LP
Other, even more specialized routines are provided for controlling the
parsing process, printing out expressions and sifting through stored
definitions, but these are not accessed by most client programs.
Worth noting are the various compile flags that affect which features of
the expression language are included.
The standard library sets the flags \f[C]-DVARIABLE\f[R],
\f[C]-DFUNCTION\f[R], \f[C]-DRCONST\f[R], and \f[C]-DBIGLIB\f[R].
Here is a list of compile flags and their meanings:
.IP "\f[B]\f[CB]-DVARIABLE\f[B]\f[R]" 3
Allow user-defined variables and (if \f[C]-DFUNCTION\f[R]) user- defined
functions.
.RS 3
.RE
.IP "\f[B]\f[CB]-DFUNCTION\f[B]\f[R]" 3
Compile in library functions and (if \f[C]-DVARIABLE\f[R]) allow
user-supplied function definitions.
.RS 3
.RE
.IP "\f[B]\f[CB]-DBIGLIB\f[B]\f[R]" 3
Include larger library of standard functions, i.e., standard C math
library.
Otherwise, only minimal library is compiled in, and other functions may
be added using \f[C]funset\f[R].
.RS 3
.RE
.IP "\f[B]\f[CB]-DRCONST\f[B]\f[R]" 3
Reduce constant subexpressions during compilation.
This can result in substantial savings during later evaluation, but the
original user-supplied expressions are lost.
.RS 3
.RE
.IP "\f[B]\f[CB]-DREDEFW\f[B]\f[R]" 3
Issue a warning via \f[C]wputs(s)\f[R] if a new definition hides a
constant definition or library function, or replaces an existing,
distinct definition for the same symbol.
(The \f[C]varset\f[R] routine never generates warnings, however.)
.RS 3
.RE
.IP "\f[B]\f[CB]-DINCHAN\f[B]\f[R]" 3
Provide for \f[C]\[rs]$N\f[R] syntax for input channels, which result in
callbacks to client-supplied chanvalue(n) routine on each evaluation.
.RS 3
.RE
.IP "\f[B]\f[CB]-DOUTCHAN\f[B]\f[R]" 3
Provide for \f[C]\[rs]$N\f[R] lvalue syntax for output channels, which
are evaluated via the chanout(cs) library function, which calls
\f[C](*cs)(n, value)\f[R]for each assigned channel definition.
.RS 3
.RE
.SH 2
Data File Format (.dat suffix)
.pdfhref O 2 "Data File Format \(.dat suffix\)"
.pdfhref M "data-file-format-.dat-suffix"
.LP
Although it is possible to store tabular data in a function file using
the select library function, it is more convenient and efficient to
devote a special file format to this purpose.
\f[I]Radiance\f[R] data files store scalar values on an N-dimensional
rectangular grid.
Grid (independent) axes may be regularly or irregularly divided, as
shown in Figure 1.
This data is interpolated during rendering (using N-dimensional linear
interpolation) to compute the desired values.
.PSPIC  "figure1.eps"
.ce 1000
Division of axes in .dat file
.ce 0
.PP
\f[B]Figure 1.\f[R] A 2-dimensional grid with one regularly divided axis
and one irregularly divided axis.
Each intersection corresponds to a data value that appears in the file.
.PP
Data files are broken into two sections, the header and the body.
The header specifies the grid, and the body contains the data values in
a standard order.
The first value in the file is a positive integer indicating the number
of dimensions.
Next comes that number of axis specifications, in one of two formats.
For a regularly divided axis, the starting and ending value is given,
followed by the number of divisions.
For an irregularly divided axis, two zeros are followed by the number of
divisions then that number of division values.
The two zeros are merely there to indicate an irregular spacing is being
specified.
Once all the axes have been given, the header is done and the body of
the file begins, which consists of one data value after another.
The ordering of the data is such that the last axis given is the one
being traversed most rapidly, similar to a static array assignment in C.
.PP
A file corresponding to the topology shown in Figure 1 is:
.IP
.nf
\f[C]
######### Header ########
2                   # Two-dimensional data array
0.5 0.1 5           # The regularly spaced axis
0 0 5 3 5 10 16 20  # The irregularly spaced axis
########## Body #########
# The data values, starting with the
# upper left, moving right then down:
 19.089   7.001  14.647   6.3671  8.0003
  3.8388 11.873  19.294  16.605   2.7435
 16.699   6.387   2.8123 16.195  17.615
 14.36   14.413  16.184  15.635   4.5403
  3.6740 14.550  10.332  15.932   1.2678
\f[]
.fi
.LP
Comments begin with a pound sign (`#') and continue to the end of the
line.
White space is ignored except as a data separator, thus the position of
header and data values on each line is irrelevant except to improve
readability.
.SH 3
\f[BI]Radiance\f[B] Programs
.pdfhref O 3 "Radiance Programs"
.pdfhref M "radiance-programs-2"
.LP
Table 3 shows \f[I]Radiance\f[R] programs that read and write data
files.
.PP
.na
.TS
delim(@@) tab(	);
l l l l.
T{
Program
T}	T{
Read
T}	T{
Write
T}	T{
Function
T}
_
T{
\f[B]ies2rad\f[R]
T}	T{
T}	T{
X
T}	T{
Convert IES luminaire file to \f[I]Radiance\f[R]
T}
T{
\f[B]mgf2rad\f[R]
T}	T{
T}	T{
X
T}	T{
Convert MGF file to \f[I]Radiance\f[R]
T}
T{
\f[B]rpict\f[R]
T}	T{
X
T}	T{
T}	T{
Batch rendering program
T}
T{
\f[B]rtrace\f[R]
T}	T{
X
T}	T{
T}	T{
Customizable ray-tracer
T}
T{
\f[B]rview\f[R]
T}	T{
X
T}	T{
T}	T{
Interactive renderer
T}
.TE
.ad
.LP
\f[B]Table 3.\f[R] Programs in the \f[I]Radiance\f[R] distribution that
read and write data files.
.SH 3
\f[BI]Radiance\f[B] C Library
.pdfhref O 3 "Radiance C Library"
.pdfhref M "radiance-c-library-2"
.LP
The header file \f[C]src/rt/data.h\f[R] gives the standard data
structures used by the routines in \f[C]src/rt/data.c\f[R] for reading
and interpolating data files.
The main data type is \f[C]DATARRAY\f[R], which is a structure
containing the grid specification and a pointer to the data array, which
is of the type \f[C]DATATYPE\f[R] (normally \f[B]float\f[R] to save
space).
.PP
The main routine for reading data files is \f[C]getdata(dname)\f[R],
which searches the standard \f[I]Radiance\f[R] library locations set by
the \f[C]RAYPATH\f[R] environment variable.
The return value is a pointer to the loaded \f[C]DATARRAY\f[R], which
may have been loaded by a previous call.
(The routine keeps a hash table of loaded files to minimize time and
memory requirements.)
The \f[C]freedata(dname)\f[R] routine frees memory associated with the
named data file, or all data arrays if \f[C]dname\f[R] is
\f[C]NULL\f[R].
.PP
The routine that interpolates data values is \f[C]datavalue(dp,pt)\f[R],
which takes a \f[C]DATARRAY\f[R] pointer and an array of
\f[B]double\f[R]s of the appropriate length (the number of dimensions in
\f[C]dp\f[R]).
The \f[B]double\f[R] returned is the interpolated value at that point in
the scalar field.
If the requested point lies outside the data\[cq]s grid, it is
extrapolated from the perimeter values up to a distance of one division
in each dimension, and falls off harmonically to zero outside of that.
This was selected as the most robust compromise, but to be safe it is
generally best to avoid giving out-of-domain points to
\f[C]datavalue\f[R].
.SH 2
Font File Format (.fnt suffix)
.pdfhref O 2 "Font File Format \(.fnt suffix\)"
.pdfhref M "font-file-format-.fnt-suffix"
.LP
Font files are used for text patterns and mixtures, and by the
\f[B]psign\f[R] program to generate simple text labels.
Each character glyph is set up on a simple rectangular coordinate system
running from [0,255] in x and y, and the glyph itself is a polygon.
Figure 2 shows an example of the letter \[lq]A\[rq].
.PSPIC  "figure2.eps"
.ce 1000
Drawing of an A, based on coordinates
.ce 0
.PP
\f[B]Figure 2.\f[R] A glyph for an \[lq]A\[rq] character in standard
font coordinates.
Note that the hole is made via a seam, just as with \f[I]Radiance\f[R]
scene polygons.
The actual aspect and spacing of the character will be determined by the
client program.
.PP
Each glyph begins with the decimal value of that character\[cq]s index,
which is 65 for \[lq]A\[rq] according to the ASCII standard.
This is followed by the number of vertices, then the vertices themselves
in @left ( x sub 1 , y sub 1 right ) , left ( x sub 2 , y sub 2 right )@
order.
White space again serves as a separator, and comments may begin with a
pound sign (`#') and continue to the end of line.
Here is the glyph entry for the letter \[lq]A\[rq] corresponding to
Figure 2:
.IP
.nf
\f[C]
65   15        # Helvetica \[dq]A\[dq]
     155  222  242  48   185  48   168  86
     83   86   65   48   12   48   101  222
     155  222  128  179  126  179  97   116
     155  116  128  179  155  222
\f[]
.fi
.LP
If the number of vertices given is zero, then the character is a space.
This is not the same as no entry, which means there is no valid glyph
for that character.
Glyphs may appear in any order, with indices ranging from 0 to 255.
The maximum number of vertices for a single glyph is 32767.
.PP
Two standard font files are provided, \f[C]helvet.fnt\f[R] and
\f[C]hexbit4x1.fnt\f[R].
The former is a Helvetica font from the public domain Hershey font set.
The second is a simple bit pattern font for hexadecimal encodings of
bitmaps.
.SH 3
\f[BI]Radiance\f[B] Programs
.pdfhref O 3 "Radiance Programs"
.pdfhref M "radiance-programs-3"
.LP
Table 4 shows \f[I]Radiance\f[R] programs that read and write font
files.
.PP
.na
.TS
delim(@@) tab(	);
l l l l.
T{
Program
T}	T{
Read
T}	T{
Write
T}	T{
Function
T}
_
T{
\f[B]pcompos\f[R]
T}	T{
X
T}	T{
T}	T{
Compose \f[I]Radiance\f[R] pictures
T}
T{
\f[B]psign\f[R]
T}	T{
X
T}	T{
T}	T{
Generate \f[I]Radiance\f[R] picture label
T}
T{
\f[B]rpict\f[R]
T}	T{
X
T}	T{
T}	T{
Batch rendering program
T}
T{
\f[B]rtrace\f[R]
T}	T{
X
T}	T{
T}	T{
Customizable ray-tracer
T}
T{
\f[B]rview\f[R]
T}	T{
X
T}	T{
T}	T{
Interactive renderer
T}
.TE
.ad
.LP
\f[B]Table 4.\f[R] Programs in the \f[I]Radiance\f[R] distribution that
read and write font files.
.SH 3
\f[BI]Radiance\f[B] C Library
.pdfhref O 3 "Radiance C Library"
.pdfhref M "radiance-c-library-3"
.LP
Similar to data files, font files are usually read and stored in a
lookup table.
The data structures for fonts are in src/common/font.h, and the routines
for reading and spacing them are in \f[C]src/common/font.c\f[R].
The main structure type is \f[C]FONT\f[R].
The routine \f[C]getfont(fname)\f[R] loads a font file from the
\f[I]Radiance\f[R] library (set by the \f[C]RAYPATH\f[R] environment
variable), and returns a pointer to the resulting \f[C]FONT\f[R]
structure.
If the file has been previously loaded, a pointer to the stored
structure is returned.
The \f[C]freefont(fname)\f[R] routine frees memory associated with the
named font file and deletes it from the table, or frees all font data if
\f[C]fname\f[R] is \f[C]NULL\f[R].
.PP
Three different routines are available for text spacing.
The \f[C]uniftext(sp,tp,f\f[R]) function takes the nul-terminated string
\f[C]tp\f[R] and computes uniform per-character spacing for the font
\f[C]f\f[R], returned in the short integer array \f[C]sp\f[R].
(This is a fairly simple process, and all spacing values will be 255
unless a character has no corresponding glyph.)
The \f[C]squeeztext(sp,tp,f,cis)\f[R] performs a similar function, but
puts only \f[C]ci\f[R]s units between adjacent characters, based on the
actual width of each font glyph.
The most sophisticated spacing function is
\f[C]proptext(sp,tp,f,cis,nsi)\f[R], which produces a total line length
equal to what it would be with uniform spacing, while maintaining equal
inter-character spacing throughout (i.e., proportional spacing).
The \f[C]nsi\f[R] argument is the number of spaces (zero-vertex glyphs)
considered as an indent.
That is, if this many or more adjacent spaces occur in \f[C]tp\f[R], the
indented text following will appear at the same point as it would have
had the spacing been uniform.
This maintains columns in tabulated text despite the proportional
spacing.
Tabs are not understood or interpreted by any of these routines, and
must be expanded to the appropriate number of spaces via
\f[B]expand\f[R].
.SH 3
Octree Format (.oct suffix)
.pdfhref O 3 "Octree Format \(.oct suffix\)"
.pdfhref M "octree-format-.oct-suffix"
.LP
In \f[I]Radiance\f[R], octrees are used to accelerate ray intersection
calculations as described by Glassner [Glassner84].
This data structure is computed by the \f[B]oconv\f[R] program, which
produces a binary file as its output.
An octree file contains a list of \f[I]Radiance\f[R] scene description
files (which may be empty), some information to guarantee portability
between systems and different versions of the code, followed by the
octree data itself.
If the octree file is \[lq]frozen,\[rq] then it will also contain the
scene data, compiled into a binary format for quick loading.
This is most convenient for octrees that are used in \f[I]instance\f[R]
primitives, which may be moved to a different (library) location from
the originating scene files.
.PP
An octree recursively subdivides 3-dimensional space into 8 subtrees,
hence its name.
Each \[lq]leaf\[rq] node contains between zero and \f[C]MAXSET\f[R]
surface primitives, indicating that section of space contains part or
all of those surfaces.
(Surface primitives may appear more than once in the octree.)
The goal of \f[B]oconv\f[R] is to build an octree that contains no more
than N surfaces in each leaf node, where N is set by the \f[B]-n\f[R]
option (5 by default).
It may allow more surfaces in places where the octree has reached its
maximum resolution (depth), set by the \f[B]-r\f[R] option (1024 \[em]
depth 10 by default).
Figure 3 shows a quadtree dividing 2-dimensional space, analogous to our
3-dimensional octree.
.PSPIC  "figure3.eps"
.ce 1000
A quadtree dividing two-dimensional space
.ce 0
.PP
\f[B]Figure 3.\f[R] An example quadtree divided so that no leaf node
contains more than 2 objects.
A three-dimensional octree works the same way.
Each leaf node is either empty, or contains a list of intersecting
surfaces.
.SH 3
Basic File Structure
.pdfhref O 3 "Basic File Structure"
.pdfhref M "basic-file-structure-1"
.LP
An octree file is divided into five sections: the information header,
the scene boundaries, the scene file names, the octree data structure,
and the compiled scene data.
If the octree is frozen, then the compiled scene data is included and
the scene file names are not.
Otherwise, the scene data is left off.
.SH 4
Information Header
.pdfhref O 4 "Information Header"
.pdfhref M "information-header"
.LP
As with other binary \f[I]Radiance\f[R] formats, the beginning of an
octree file is the information header.
The first line is \[lq]#?RADIANCE\[rq] to aid in identification by the
UNIX \f[B]file\f[R] program.
Following this is the \f[B]oconv\f[R] command (or commands) used to
produce the octree, then a line indicating the format,
\f[C]FORMAT=Radiance_octree\f[R].
The end of the information header is always an empty line.
Here is an example of an octree information header, as reported by
\f[B]getinfo\f[R]:
.IP
.nf
\f[C]
#?RADIANCE
oconv model.b90 desk misc
oconv -f -i modelb.oct window blinds lights lamp
FORMAT=Radiance_octree
\f[]
.fi
.LP
The actual content of this header is ignored when an octree is read
except for the \f[C]FORMAT\f[R] line, which if it appears must match the
one shown above.
.SH 4
Scene Boundaries
.pdfhref O 4 "Scene Boundaries"
.pdfhref M "scene-boundaries"
.LP
After the information header, there is a magic number indicating the
format version and the size of object indices (in bytes per index).
This is a two-byte quantity, which must be one of the following in the
current release:
.PP
.na
.nr LLold \n[LL]
.TS
delim(@@) tab(	);
lw(35.0n) lw(35.0n).
T{
.nr LL 35.0n
285
T}	T{
.nr LL 35.0n
Two-byte object indices.
T}
T{
.nr LL 35.0n
287
T}	T{
.nr LL 35.0n
Four-byte object indices.
T}
T{
.nr LL 35.0n
291
T}	T{
.nr LL 35.0n
Eight-byte object indices.
(Only supported on architectures with 64-bit \f[B]longs\f[R].)
T}
.TE
.nr LL \n[LLold]
.ad
.LP
Technically, the code will also support odd-sized integers, but they are
not found on any existing machine architectures so we can forget about
them.
.PP
Following the octree magic number, we have the enclosing cube for the
scene, which defines the dimensions of the octree\[cq]s root node.
The cube is aligned along the world coordinate axes, so may be defined
by one corner (the 3-dimensional minimum) and the side length.
For historical reasons, these four values are stored as ASCII-encoded
real values in nul-terminated strings.
(The octree boundaries may also be read using \f[B]getinfo\f[R] with the
\f[B]-d\f[R] option.)
.SH 4
Scene File Names
.pdfhref O 4 "Scene File Names"
.pdfhref M "scene-file-names"
.LP
Following the octree dimensions, the names of the scene description
files are given, each stored a nul-terminated string.
The end of this file list is indicated by an empty string.
If the octree is \[lq]frozen,\[rq] meaning it contains the compiled
scene information as well, then no file names will be present (i.e., the
first string will be empty).
.SH 4
Octree Data Structure
.pdfhref O 4 "Octree Data Structure"
.pdfhref M "octree-data-structure"
.LP
After the scene description files, an N-byte integer indicates the total
number of primitives given to \f[B]oconv\f[R], where N is the size
derived from the magic number as we described.
This object count will be used to verify that the files have not changed
significantly since the octree was written\**.
.FS
Small changes that do not affect geometry will not cause problems, but
if the primitive count changes, so does the indexing of surfaces, and
with that the octree data structure becomes invalid.
A second check is made to insure that no non-surface primitives appear
in any leaf nodes, and this at least guarantees that the renderer will
not dump core from an outdated octree, even if the results are wrong.
.FE
.PP
After the primitive count, the actual octree is stored, using the
following recursive procedure:
.IP
.nf
\f[C]
puttree(ot) begin
if ot is a tree then
               write the character \[aq]\[rs]002\[aq]
call puttree on each child node (0-7) else if ot is empty then
               write the character \[aq]\[rs]000\[aq]
          else
               write the character \[aq]\[rs]001\[aq]
               write out the number of surfaces
               write out each surface\[aq]s index
end
end puttree
\f[]
.fi
.LP
The number of surfaces and the surface indices are each N-byte integers,
and the tree node types are single bytes.
Reading the octree is accomplished with a complementary procedure.
.SH 4
Compiled Scene Data
.pdfhref O 4 "Compiled Scene Data"
.pdfhref M "compiled-scene-data"
.LP
If the octree is frozen, then this data structure is followed by a
compiled version of the scene.
This avoids the problems of changes to scene files, and allows an octree
to be moved conveniently from one location and one system to another
without worrying about the associated scene files.
.PP
The scene data begins with a listing of the defined primitive types.
This list consists of the name of each type as a nul-terminated string,
followed by an empty string once the list has been exhausted.
This permits the indexing of primitive types with a single byte later
on, without concern about changes to \f[I]Radiance\f[R] involving
\f[C]src/common/otypes.h\f[R].
.PP
The scene primitives are written one at a time.
First is a single byte with the primitive type index, as determined from
the list given above.
Second is the N-byte modifier index, followed by the primitive\[cq]s
identifier as a nul-terminated string.
String arguments start with a 2-byte integer indicating the argument
count, followed by the strings themselves, which are nul-terminated.
Real arguments next have a 2-byte count followed by the real values,
each stored as a 4-byte mantissa followed by a 1-byte (signed) exponent.
(The mantissa is the numerator of a fraction of @2 sup 31 \[u2212] 1@.)
The end of data is indicated with a -1 value for the object type
(byte=255).
.SH 3
\f[BI]Radiance\f[B] Programs
.pdfhref O 3 "Radiance Programs"
.pdfhref M "radiance-programs-4"
.LP
Table 5 shows \f[I]Radiance\f[R] programs that read and write octree
files.
.PP
.na
.TS
delim(@@) tab(	);
l l l l.
T{
Program
T}	T{
Read
T}	T{
Write
T}	T{
Function
T}
_
T{
\f[B]getinfo\f[R]
T}	T{
X
T}	T{
T}	T{
Print information header from binary file
T}
T{
\f[B]oconv\f[R]
T}	T{
X
T}	T{
X
T}	T{
Compile \f[I]Radiance\f[R] scene description
T}
T{
\f[B]rad\f[R]
T}	T{
X
T}	T{
X
T}	T{
Render \f[I]Radiance\f[R] scene
T}
T{
\f[B]rpict\f[R]
T}	T{
X
T}	T{
T}	T{
Batch rendering program
T}
T{
\f[B]rpiece\f[R]
T}	T{
X
T}	T{
T}	T{
Parallel batch rendering program
T}
T{
\f[B]rtrace\f[R]
T}	T{
X
T}	T{
T}	T{
Customizable ray-tracer
T}
T{
\f[B]rview\f[R]
T}	T{
X
T}	T{
T}	T{
Interactive renderer
T}
.TE
.ad
.LP
\f[B]Table 5.\f[R] Programs in the \f[I]Radiance\f[R] distribution that
read and write octree files.
.SH 3
\f[BI]Radiance\f[B] C Library
.pdfhref O 3 "Radiance C Library"
.pdfhref M "radiance-c-library-4"
.LP
Since reading an octree file also may involve reading a
\f[I]Radiance\f[R] scene description, some of the same library routines
are called indirectly.
The header file \f[C]src/common/octree.h\f[R] is needed in addition to
the \f[C]src/common/object.\f[R]h file.
The module \f[C]src/ot/writeoct.c\f[R] contains the main routines for
writing an octree to stdout, while \f[C]src/common/readoct.c\f[R]
contains the corresponding routines for reading an octree from a file.
Both modules access routines in \f[C]src/common/portio.c\f[R] for
reading and writing portable binary data.
.PP
Here is the main call for writing out an octree:
.IP "\f[B]\f[CB]writeoct(int store, CUBE *scene, char *ofn[]);\f[B]\f[R]" 3
Write the octree stored in scene to stdout, assuming the header has
already been written.
The flags in store determine what will be included.
Normally, this variable is one of \f[C]IO_ALL\f[R] or
\f[C](IO_ALL & \[ti]IO_FILES)\f[R] correspondingto writing a normal or a
frozen octree, respectively.
.RS 3
.RE
.LP
Here is the main call for reading in an octree:
.IP "\f[B]\f[CB]readoct(char *fname, int load, CUBE *scene, char *ofn[]);\f[B]\f[R]" 3
Read the octree file fname into scene, saving scene file names in the
ofn array.
What is loaded depends on the flags in load,which may be one or more of
\f[C]IO_CHECK\f[R], \f[C]IO_INFO\f[R], \f[C]IO_SCENE\f[R],
\f[C]IO_TREE\f[R], \f[C]IO_FILES\f[R] and \f[C]IO_BOUNDS\f[R].
These correspond to checking file type and consistency, transferring the
information header to stdout, loading the scene data, loading the octree
structure, assigning the scene file names to ofn, and assigning the
octree cube boundaries.
The macro \f[C]IO_ALL\f[R] includes all of these flags, for convenience.
.RS 3
.RE
.SH 2
Picture File Format (.hdr suffix)
.pdfhref O 2 "Picture File Format \(.hdr suffix\)"
.pdfhref M "picture-file-format-.hdr-suffix"
.LP
\f[I]Radiance\f[R] pictures\**
.FS
The picture filename extension used to be .pic, but that conflicted with
too many other programs.
It was replaced with .hdr, an abbreviation of \[lq]high dynamic
range.\[rq]
.FE
differ from standard computer graphics images inasmuch as they contain
real physical data, namely radiance values at each pixel.
To do this, it is necessary to maintain floating point information, and
we use a 4-byte/pixel encoding described in Chapter II.5 of Graphics
Gems II [Arvo91,p.80].
The basic idea is to store a 1-byte mantissa for each of three
primaries, and a common 1-byte exponent.
The accuracy of these values will be on the order of 1% (±1 in 200) over
a dynamic range from 10\*{-38\*} to 10\*{38\*}.
.PP
Although \f[I]Radiance\f[R] pictures \f[I]may\f[R] contain physical
data, they do not \f[I]necessarily\f[R] contain physical data.
If the rendering did not use properly modeled light sources, or the
picture was converted from some other format, or custom filtering was
applied, then the physical data will be invalid.
Table 6 lists programs that read and write \f[I]Radiance\f[R] pictures,
with pluses next to the X-marks indicating where physical data is
preserved (or at least understood).
Specifically, if the picture file read or written by a program has an
\[lq]X+\[rq], then it has maintained the physical validity of the pixels
by keeping track of any exposure or color corrections in the appropriate
header variables, described below.
.SH 3
Basic File Structure
.pdfhref O 3 "Basic File Structure"
.pdfhref M "basic-file-structure-2"
.LP
\f[I]Radiance\f[R] picture files are divided into three sections: the
information header, the resolution string, and the scanline records.
All of these must be present or the file is incomplete.
.SH 4
Information Header
.pdfhref O 4 "Information Header"
.pdfhref M "information-header-1"
.LP
The information header begins with the usual \f[C]#?RADIANCE\f[R]
identifier, followed by one or more lines containing the programs used
to produce the picture.
These commands may be interspersed with variables describing relevant
information such as the view, exposure, color correction, and so on.
Variable assignments begin on a new line, and the variable name (usually
all upper case) is followed by an equals sign (`='), which is followed
by the assigned value up until the end of line.
Some variable assignments override previous assignments in the same
header, where other assignments are cumulative.
Here are the most important variables for \f[I]Radiance\f[R] pictures:
.IP "\f[B]\f[CB]FORMAT\f[B]\f[R]" 3
A line indicating the file\[cq]s format.
At most one \f[C]FORMAT\f[R] line is allowed, and it must be assigned a
value of either \f[C]32- bit_rle_rgbe\f[R] or \f[C]32-bit_rle_xyze\f[R]
to be a valid \f[I]Radiance\f[R] picture.
.RS 3
.RE
.IP "\f[B]\f[CB]EXPOSURE\f[B]\f[R]" 3
A single floating point number indicating a multiplier that has been
applied to all the pixels in the file.
\f[C]EXPOSURE\f[R] values are cumulative, so the original pixel values
(i.e., radiances in w/sr/m\*{2\*}) must be derived by taking the values
in the file and dividing by all the \f[C]EXPOSURE\f[R] settings
multiplied together.
No \f[C]EXPOSURE\f[R] setting implies that no exposure changes have
taken place.
.RS 3
.RE
.IP "\f[B]\f[CB]COLORCORR\f[B]\f[R]" 3
A color correction multiplier that has been applied to this picture.
Similar to the \f[C]EXPOSURE\f[R] variable except given in three parts
for the three primaries.
In general, the value should have a brightness of unity, so that it does
not affect the actual brightness of pixels, which should be tracked by
\f[C]EXPOSURE\f[R] changes instead.
(This variable is also cumulative.)
.RS 3
.RE
.IP "\f[B]\f[CB]SOFTWARE\f[B]\f[R]" 3
The software version used to create the picture, usually something like
\f[C]RADIANCE 3.04 official release July 16, 1996\f[R].
.RS 3
.RE
.IP "\f[B]\f[CB]PIXASPECT\f[B]\f[R]" 3
The pixel aspect ratio, expressed as a decimal fraction of the height of
each pixel to its width.
This is not to be confused with the image aspect ratio, which is the
total height over width.
(The image aspect ratio is actually equal to the height in pixels over
the width in pixels, \f[I]multiplied\f[R] by the pixel aspect ratio.)
These assignments are cumulative, so the actual pixel aspect ratio is
all ratios multiplied together.
If no \f[C]PIXASPECT\f[R] assignment appears, the ratio is assumed to be
1.
.RS 3
.RE
.IP "\f[B]\f[CB]VIEW\f[B]\f[R]" 3
The \f[I]Radiance\f[R] view parameters used to create this picture.
Multiple assignments are cumulative inasmuch as new view options add to
or override old ones.
.RS 3
.RE
.LP
\f[C]PRIMARIES\f[R] The CIE (x,y) chromaticity coordinates of the three
(RGB) primaries and the white point used to standardize the
picture\[cq]s color system.
This is used mainly by the \f[B]ra_xyze\f[R] program to convert between
color systems.
If no \f[C]PRIMARIES\f[R] line appears, we assume the standard primaries
defined in \f[C]src/common/color.h\f[R], namely
\f[C]0.640 0.330 0.290 0.600 0.150 0.060 0.333 0.333\f[R] for red,
green, blue and white, respectively.
.PP
As always, the end of the header is indicated by an empty line.
.SH 4
Resolution String
.pdfhref O 4 "Resolution String"
.pdfhref M "resolution-string"
.LP
All \f[I]Radiance\f[R] pictures have a standard coordinate system, which
is shown in Figure 4.
The origin is always at the lower left corner, with the X coordinate
increasing to the right, and the Y coordinate increasing in the upward
direction.
The actual ordering of pixels in the picture file, however, is addressed
by the resolution string.
.PSPIC  "figure4.eps"
.ce 1000
Radiance picture coordinate system
.ce 0
.PP
\f[B]Figure 4.\f[R] The standard coordinate system for an MxN picture.
.PP
The resolution string is given as one of the following:
.IP "\f[B]\f[CB]-Y N +X M\f[B]\f[R]" 3
The standard orientation produced by the renderers, indicating that Y is
decreasing in the file, and X is increasing.
X positions are increasing in each scanline, starting with the upper
left position in the picture and moving to the upper right initially,
then on down the picture.
Some programs will only handle pictures with this ordering.
.RS 3
.RE
.IP "\f[B]\f[CB]-Y N -X M\f[B]\f[R]" 3
The X ordering has been reversed, effectively flipping the image left to
right from the standard ordering.
.RS 3
.RE
.IP "\f[B]\f[CB]+Y N -X M\f[B]\f[R]" 3
The image has been flipped left to right and also top to bottom, which
is the same as rotating it by 180 degrees.
.RS 3
.RE
.IP "\f[B]\f[CB]+Y N +X M\f[B]\f[R]" 3
The image has been flipped top to bottom from the standard ordering.
.RS 3
.RE
.IP "\f[B]\f[CB]+X M +Y N\f[B]\f[R]" 3
The image has been rotated 90 degrees clockwise.
.RS 3
.RE
.IP "\f[B]\f[CB]-X M +Y N\f[B]\f[R]" 3
The image has been rotated 90 degrees clockwise, then flipped top to
bottom.
.RS 3
.RE
.IP "\f[B]\f[CB]-X M -Y N\f[B]\f[R]" 3
The image has been rotated 90 degrees counter-clockwise.
.RS 3
.RE
.IP "\f[B]\f[CB]+X M -Y N\f[B]\f[R]" 3
The image has been rotate 90 degrees counter-clockwise, then flipped top
to bottom.
.RS 3
.RE
.LP
The reason for tracking all these changes in picture orientation is so
programs that compute ray origin and direction from the \f[C]VIEW\f[R]
variable in the information header will work despite such changes.
Also, it can reduce memory requirements on converting from other image
formats that have a different scanline ordering, such as Targa.
.SH 4
Scanline Records
.pdfhref O 4 "Scanline Records"
.pdfhref M "scanline-records"
.LP
\f[I]Radiance\f[R] scanlines come in one of three flavors, described
below:
.IP "\f[B]Uncompressed\f[R]" 3
Each scanline is represented by M pixels with 4 bytes per pixel, for a
total length of 4xM bytes.
This is the simplest format to read and write, since it has a one-to-one
correspondence to an array of \f[C]COLR\f[R] values.
.RS 3
.RE
.IP "\f[B]Old run-length encoded\f[R]" 3
Repeated pixel values are indicated by an illegal (i.e., unnormalized)
pixel that has 1\[cq]s for all three mantissas, and an exponent that
corresponds to the number of times the previous pixel is repeated.
Consecutive repeat indicators contain higher-order bytes of the count.
.RS 3
.RE
.IP "\f[B]New run-length encoded\f[R]" 3
In this format, the four scanline components (three primaries and
exponent) are separated for better compression using adaptive run-length
encoding (described by Glassner in Chapter II.8 of Graphics Gems II
[Arvo91,p.89]).
The record begins with an unnormalized pixel having two bytes equal to
2, followed by the upper byte and the lower byte of the scanline length
(which must be less than 32768).
A run is indicated by a byte with its high-order bit set, corresponding
to a count with excess 128.
A non-run is indicated with a byte less than 128.
The maximum compression ratio using this scheme is better than 100:1,
but typical performance for \f[I]Radiance\f[R] pictures is more like
2:1.
.RS 3
.RE
.LP
The physical values these scanlines correspond to depend on the format
and other information contained in the information header.
If the \f[C]FORMAT\f[R] string indicates RGB data, then the units for
each primary are spectral radiances over the corresponding waveband,
such that a pixel value of @left ( 1 , 1 , 1 right )@ corresponds to a
total energy of 1 w/sr/m\*{2\*} over the visible spectrum.
The actual luminance value (in lm/sr/m\*{2\*}) can be computed from the
following formula for the standard \f[I]Radiance\f[R] RGB primaries:
.PP
.EQ
L sub v = 179 left ( 0.265 r + 0.670 g + 0.065 b right )
.EN
.PP
The value of 179 lm/w is the standard \f[I]luminous efficacy\f[R] of
equal-energy white light that is defined and used by \f[I]Radiance\f[R]
specifically for this conversion.
This and the other values above are defined in
\f[C]src/common/color.h\f[R], and the above formula is given as a macro,
\f[C]luminance(col)\f[R].
.PP
If the \f[C]FORMAT\f[R] string indicates XYZ data, then the units of the
Y primary are already lm/st/m\*{2\*}, so the above conversion is
unnecessary.
.SH 3
\f[BI]Radiance\f[B] programs
.pdfhref O 3 "Radiance programs"
.pdfhref M "radiance-programs-5"
.LP
Table 6 shows the many programs that read and write \f[I]Radiance\f[R]
pictures.
.PP
.na
.TS
delim(@@) tab(	);
l l l l.
T{
Program
T}	T{
Read
T}	T{
Write
T}	T{
Function
T}
_
T{
\f[B]falsecolor\f[R]
T}	T{
X+
T}	T{
X
T}	T{
Create false color image
T}
T{
\f[B]findglare\f[R]
T}	T{
X+
T}	T{
T}	T{
Find sources of discomfort glare
T}
T{
\f[B]getinfo\f[R]
T}	T{
X
T}	T{
T}	T{
Print information header from binary file
T}
T{
\f[B]macbethcal\f[R]
T}	T{
X
T}	T{
X
T}	T{
Compute image color & contrast correction
T}
T{
\f[B]normpat\f[R]
T}	T{
X
T}	T{
X
T}	T{
Normalize picture for use as pattern tile
T}
T{
\f[B]objpict\f[R]
T}	T{
T}	T{
X
T}	T{
Generate composite picture of object
T}
T{
\f[B]pcomb\f[R]
T}	T{
X+
T}	T{
X
T}	T{
Perform arbitrary math on picture(s)
T}
T{
\f[B]pcond\f[R]
T}	T{
X+
T}	T{
X
T}	T{
Condition a picture for display
T}
T{
\f[B]pcompos\f[R]
T}	T{
X
T}	T{
X
T}	T{
Composite pictures
T}
T{
\f[B]pextrem\f[R]
T}	T{
X+
T}	T{
T}	T{
Find minimum and maximum pixels
T}
T{
\f[B]pfilt\f[R]
T}	T{
X+
T}	T{
X+
T}	T{
Filter and anti-alias picture
T}
T{
\f[B]pflip\f[R]
T}	T{
X+
T}	T{
X+
T}	T{
Flip picture left-right and/or top-bottom
T}
T{
\f[B]pinterp\f[R]
T}	T{
X+
T}	T{
X+
T}	T{
Interpolate/extrapolate picture views
T}
T{
\f[B]protate\f[R]
T}	T{
X+
T}	T{
X+
T}	T{
Rotate picture 90 degrees clockwise
T}
T{
\f[B]psign\f[R]
T}	T{
T}	T{
X
T}	T{
Create text picture
T}
T{
\f[B]pvalue\f[R]
T}	T{
X+
T}	T{
X+
T}	T{
Convert picture to/from simpler formats
T}
T{
\f[B]ra_avs\f[R]
T}	T{
X
T}	T{
X
T}	T{
Convert to/from AVS image format
T}
T{
\f[B]ra_pict\f[R]
T}	T{
X
T}	T{
X
T}	T{
Convert to/from Macintosh PICT2 format
T}
T{
\f[B]ra_ppm\f[R]
T}	T{
X
T}	T{
X
T}	T{
Convert to/from Poskanzer Port.
Pixmap
T}
T{
\f[B]ra_pr\f[R]
T}	T{
X
T}	T{
X
T}	T{
Convert to/from Sun 8-bit rasterfile
T}
T{
\f[B]ra_pr24\f[R]
T}	T{
X
T}	T{
X
T}	T{
Convert to/from Sun 24-bit rasterfile
T}
T{
\f[B]ra_ps\f[R]
T}	T{
X
T}	T{
T}	T{
Convert to B&W or color PostScript
T}
T{
\f[B]ra_rgbe\f[R]
T}	T{
X
T}	T{
X
T}	T{
Convert to/from uncompressed picture
T}
T{
\f[B]ra_t8\f[R]
T}	T{
X
T}	T{
X
T}	T{
Convert to/from Targa 8-bit format
T}
T{
\f[B]ra_t16\f[R]
T}	T{
X
T}	T{
X
T}	T{
Convert to/from Targa 16-bit and 24-bit
T}
T{
\f[B]ra_tiff\f[R]
T}	T{
X
T}	T{
X
T}	T{
Convert to/from TIFF 8-bit and 24-bit
T}
T{
\f[B]ra_xyze\f[R]
T}	T{
X
T}	T{
X
T}	T{
Convert to/from CIE primary picture
T}
T{
\f[B]rad\f[R]
T}	T{
T}	T{
X+
T}	T{
Render Radiance scene
T}
T{
\f[B]ranimate\f[R]
T}	T{
T}	T{
X+
T}	T{
Animate Radiance scene
T}
T{
\f[B]rpict\f[R]
T}	T{
X
T}	T{
X+
T}	T{
Batch rendering program
T}
T{
\f[B]rpiece\f[R]
T}	T{
X
T}	T{
X+
T}	T{
Parallel batch rendering program
T}
T{
\f[B]rtrace\f[R]
T}	T{
X
T}	T{
X+
T}	T{
Customizable ray-tracer
T}
T{
\f[B]rview\f[R]
T}	T{
X
T}	T{
X+
T}	T{
Interactive renderer
T}
T{
\f[B]vwright\f[R]
T}	T{
X
T}	T{
T}	T{
Get view parameters and shift them
T}
T{
\f[B]xglaresrc\f[R]
T}	T{
X
T}	T{
T}	T{
Display glare sources from \f[B]findglare\f[R]
T}
T{
\f[B]ximage\f[R]
T}	T{
X+
T}	T{
T}	T{
Display \f[I]Radiance\f[R] picture under \f[I]X11\f[R]
T}
T{
\f[B]xshowtrace\f[R]
T}	T{
X
T}	T{
T}	T{
Show ray traces on \f[I]X11\f[R] display
T}
.TE
.ad
.LP
\f[B]Table 6.\f[R] \f[I]Radiance\f[R] programs that read and write
picture files.
Pluses indicate when a program makes use of or preserves physical pixel
values.
.SH 3
\f[BI]Radiance\f[B] C Library
.pdfhref O 3 "Radiance C Library"
.pdfhref M "radiance-c-library-5"
.LP
There are a fair number of routines for reading, writing and
manipulating \f[I]Radiance\f[R] pictures.
If you want to write a converter to or from a 24-bit image format, you
can follow the skeletal example in \f[C]src/px/ra_skel.c\f[R].
This has all of the basic functionality of the other \f[B]ra_*\f[R]
image conversion programs, with the actual code for the destination type
removed (or simplified).
The method in \f[C]ra_skel\f[R] uses the routines in
\f[C]src/common/colrops.c\f[R] to avoid conversion to machine floating
point, which can slow things down and is not necessary in this case.
.PP
Below we describe routines for reading and writing pictures, which rely
heavily on definitions in \f[C]src/common/color.h\f[R].
We start with the calls for manipulating information headers, followed
by the calls for resolution strings, then the calls for scanline
records.
.PP
Information headers are manipulated with the routines in
\f[C]src/common/header.c\f[R] and the macros in \f[C]color.h\f[R].
Features for handing views are defined in \f[C]src/common/view.h\f[R]
with routines in \f[C]src/common/image.c\f[R].
Here are the relevant calls for reading and copying information headers:
.IP "\f[B]\f[CB]int checkheader(FILE *fin, char *fmt, FILE *fout);\f[B]\f[R]" 3
Read the header information from \f[C]fin\f[R], copying to
\f[C]fou\f[R]t (unless fout is \f[C]NULL\f[R]), checking any
\f[C]FORMAT\f[R] line against the string \f[C]fmt\f[R].
The \f[C]FORMAT\f[R] line (if it exists) will not be copied to
\f[C]fout\f[R].
The function returns 1 if the header was OK and the format matched, 0 if
the header was OK but there was no format line, and -1 if the format
line did not match or there was some problem reading the header.
Wildcard characters (`*' and `?') may appear in \f[C]fmt\f[R], in which
case a globbing match is applied, and the matching format value will be
copied to fmt upon success.
The normal \f[C]fmt\f[R] values for pictures are \f[C]COLRFMT\f[R] for
\f[I]Radiance\f[R] RGB, \f[C]CIEFMT\f[R] for 4-byte XYZ pixels, or a
copy of \f[C]PICFMT\f[R] for glob matches to either.
(Do not pass \f[C]PICFMT\f[R] directly, as this will cause an illegal
memory access on systems that store static strings in read-only memory.)
.RS 3
.RE
.IP "\f[B]\f[CB]int getheader(FILE *fp, int (*f)(), char *p);\f[B]\f[R]" 3
For those who need more control, \f[C]getheader\f[R] reads the header
from \f[C]fp\f[R], calling the function \f[C]f\f[R] (if
not\f[C]NULL\f[R]) with each input line and the client data pointer
\f[C]p\f[R].
A simple call to skip the header is \f[C]getheader(fp,NULL,NULL)\f[R].
To copy the header unconditionally to \f[C]stdout\f[R], call
\f[C]getheader(fp,fputs,stdout)\f[R].
More often, \f[C]getheader\f[R] is called with a client function, which
checks each line for specific variable settings.
.RS 3
.RE
.IP "\f[B]\f[CB]int isformat(char *s);\f[B]\f[R]" 3
Returns non-zero if the line \f[C]s\f[R] is a \f[C]FORMAT\f[R]
assignment.
.RS 3
.RE
.IP "\f[B]\f[CB]int formatval(char *r, char *s);\f[B]\f[R]" 3
Returns the \f[C]FORMAT\f[R] value from line \f[C]s\f[R] in the string
\f[C]r\f[R].
Returns non-zero if \f[C]s\f[R] is a valid format line.
.RS 3
.RE
.IP "\f[B]\f[CB]fputformat(char *s, FILE *fp);\f[B]\f[R]" 3
Put format assignment \f[C]s\f[R] to the stream \f[C]fp\f[R].
.RS 3
.RE
.IP "\f[B]\f[CB]isexpos(s)\f[B]\f[R]" 3
Macro returns non-zero if the line \f[C]s\f[R] is an \f[C]EXPOSURE\f[R]
setting.
.RS 3
.RE
.IP "\f[B]exposval(s)\f[R]" 3
Macro returns \f[B]double\f[R] exposure value from line \f[C]s\f[R].
.RS 3
.RE
.IP "\f[B]\f[CB]fputexpos(ex,fp)\f[B]\f[R]" 3
Macro** | puts real exposure value \f[C]ex\f[R] to stream \f[C]fp\f[R].
.RS 3
.RE
.LP
\f[C]iscolcor(s)\f[R] :Macro returns non-zero if the line \f[C]s\f[R] is
a \f[C]COLORCORR\f[R] setting.
.IP "\f[B]\f[CB]colcorval(cc,s)\f[B]\f[R]" 3
Macro assign color correction value from line \f[C]s\f[R] in the
\f[C]COLOR\f[R] variable \f[C]cc\f[R].
.RS 3
.RE
.IP "\f[B]\f[CB]fputcolcor(cc,fp)\f[B]\f[R]" 3
Macro puts correction \f[C]COLOR\f[R] \f[C]cc\f[R] to stream
\f[C]fp\f[R].
.RS 3
.RE
.IP "\f[B]\f[CB]isaspect(s)\f[B]\f[R]" 3
Macro returns non-zero if the line \f[C]s\f[R] is a \f[C]PIXASPECT\f[R]
setting.
.RS 3
.RE
.IP "\f[B]\f[CB]aspectval(s)\f[B]\f[R]" 3
Macro returns double pixel aspect value from line \f[C]s\f[R].
.RS 3
.RE
.IP "\f[B]\f[CB]fputaspect(pa,fp)\f[B]\f[R]" 3
Macro puts real pixel aspect value \f[C]pa\f[R] to stream \f[C]fp\f[R].
.RS 3
.RE
.IP "\f[B]\f[CB]int isview(char *s);\f[B]\f[R]" 3
Returns non-zero if header line \f[C]s\f[R] contains view parameters.
Note that \[ga]\f[C]s\f[R] could be either a \f[C]VIEW\f[R] assignment
or a rendering command.
.RS 3
.RE
.IP "\f[B]\f[CB]int sscanview(VIEW *vp, char *s);\f[B]\f[R]" 3
Scan view options from the string \f[C]s\f[R] into the \f[C]VIEW\f[R]
structure pointed to by \f[C]vp\f[R].
.RS 3
.RE
.IP "\f[B]\f[CB]fprintview(VIEW *vp, FILE *fp);\f[B]\f[R]" 3
Print view options in \f[C]vp\f[R] to the stream \f[C]fp\f[R].
Note that this does not print out \[lq]VIEW=\[rq] first, or end the
line.
Therefore, one usually calls \f[C]fputs(VIEWSTR,fp)\f[R] followed by
\f[C]fprintview(vp,fp)\f[R], then \f[C]putc(\[aq]\[rs]n\[aq],fp)\f[R].
.RS 3
.RE
.IP "\f[B]\f[CB]isprims(s)\f[B]\f[R]" 3
Macro returns non-zero if the line \f[C]s\f[R] is a \f[C]PRIMARIES\f[R]
setting.
.RS 3
.RE
.IP "\f[B]\f[CB]primsval(p,s)\f[B]\f[R]" 3
Macro assign color primitives from line \f[C]s\f[R] in the
\f[C]RGBPRIMS\f[R] variable \f[C]p\f[R].
.RS 3
.RE
.IP "\f[B]\f[CB]fputprims(p,fp)\f[B]\f[R]" 3
Macro puts color primitives \f[C]p\f[R] to stream \f[C]fp\f[R].
.RS 3
.RE
.LP
The header file \f[C]src/common/resolu.h\f[R] has macros for resolution
strings, which are handled by routines in \f[C]src/common/resolu.c\f[R].
Here are the relevant calls:
.IP "\f[B]\f[CB]fgetsresolu(rs,fp)\f[B]\f[R]" 3
Macro to get a resolution string from the stream \f[C]fp\f[R] and put it
in the \f[C]RESOLU\f[R] structure pointed to by \f[C]rs\f[R].
The return value is non-zero if the resolution was successfully loaded.
.RS 3
.RE
.IP "\f[B]\f[CB]fputsresolu(rs,fp)\f[B]\f[R]" 3
Macro to write the \f[C]RESOLU\f[R] structure pointed to by rs to the
stream fp.
.RS 3
.RE
.IP "\f[B]\f[CB]scanlen(rs)\f[B]\f[R]" 3
Macro to get the scanline length from the \f[C]RESOL\f[R]U structure
pointed to by \f[C]rs\f[R].
.RS 3
.RE
.IP "\f[B]\f[CB]numscans(rs)\f[B]\f[R]" 3
Macro to get the number of scanlines from the \f[C]RESOLU\f[R] structure
pointed to by \f[C]rs\f[R].
.RS 3
.RE
.IP "\f[B]\f[CB]fscnresolu(slp,nsp,fp)\f[B]\f[R]" 3
Macro to read in a resolution string from \f[C]f\f[R]p and assign the
scanline length and number of scanlines to the integers pointed to by
\f[C]slp\f[R] and \f[C]nsp\f[R], respectively.
This call expects standard English-text ordered scanlines, and returns
non-zero only if the resolution string agrees.
.RS 3
.RE
.IP "\f[B]\f[CB]fprtresolu(sl,ns,fp)\f[B]\f[R]" 3
Macro to print out a resolution string for \f[C]ns\f[R] scanlines of
length \f[C]sl\f[R] in standard English-text ordering to \f[C]fp\f[R].
.RS 3
.RE
.LP
The file \f[C]src/common/color.c\f[R] contains the essential routines
for reading and writing scanline records.
Here are the relevant calls and macros:
.IP "\f[B]\f[CB]int freadcolrs(COLR *scn, int sl, FILE *fp);\f[B]\f[R]" 3
Read a scanline record of length \f[C]sl\f[R] from stream \f[C]fp\f[R]
into the \f[C]COLR\f[R] array \f[C]scn\f[R].
Interprets uncompressed, old, and new run-length encoded records.
Returns 0 on success, -1 on failure.
.RS 3
.RE
.IP "\f[B]\f[CB]int fwritecolrs(COLR *scn, int sl, FILE *fp);\f[B]\f[R]" 3
Write the scanline record stored in the \f[C]COLR\f[R] array
\f[C]scn\f[R], length \f[C]sl\f[R], to the stream \f[C]fp\f[R].
Uses the new run-length encoding unless \f[C]sl\f[R] is less than 8 or
greater than 32767, when an uncompressed record is written.
Returns 0 on success, -1 if there was an error.
.RS 3
.RE
.IP "\f[B]\f[CB]int freadscan(COLOR *fscn, int sl, FILE *fp);\f[B]\f[R]" 3
Reads a scanline of length \f[C]sl\f[R] from the stream \f[C]fp\f[R] and
converts to machine floating-point format in the array \f[C]fscn\f[R].
Recognizes all scanline record encodings.
Returns 0 on success, or -1 on failure.
.RS 3
.RE
.IP "\f[B]\f[CB]int fwritescan(COLOR *fscn, int sl, FILE *fp);\f[B]\f[R]" 3
Write the floating-point scanline of length \f[C]sl\f[R] stored in the
array fscn to the stream \f[C]fp\f[R].
Converts to 4-byte/pixel scanline format and calls \f[C]fwritecolrs\f[R]
to do the actual write.
Returns 0 on success, or -1 if there was an error.
.RS 3
.RE
.IP "\f[B]\f[CB]colval(col,p)\f[B]\f[R]" 3
Macro to get primary p from the floating-point \f[C]COLOR col\f[R].
The primary index may be one of \f[C]RED\f[R], \f[C]GRN\f[R] or
\f[C]BLU\f[R] for RGB colors, or CIEX, CIEY or CIEZ for XYZ colors.
This macro is a valid lvalue, so can be used on the left of assignment
statements as well.
.RS 3
.RE
.IP "\f[B]\f[CB]colrval(clr,p)\f[B]\f[R]" 3
Macro to get primary \f[C]p\f[R] from the 4-byte \f[C]COLR\f[R] pixel
\f[C]clr\f[R].
The primary index may be one of RED, GRN or BLU for RGB colors, or CIEX,
CIEY or CIEZ for XYZ colors.
Unless just one primary is needed, it is more efficient to call
\f[C]colr_color\f[R] and use the \f[C]colval\f[R] macro on the result.
.RS 3
.RE
.IP "\f[B]\f[CB]colr_color(COLOR col, COLR clr);\f[B]\f[R]" 3
Convert the 4-byte pixel \f[C]clr\f[R] to a machine floating-point
representation and store the result in \f[C]col\f[R].
.RS 3
.RE
.IP "\f[B]\f[CB]setcolr(COLR clr, double p1, p2, p3);\f[B]\f[R]" 3
Assign the 4-byte pixel \f[C]clr\f[R] according to the three primary
values \f[C]p1\f[R], \f[C]p2\f[R] and \f[C]p3\f[R].
These can be either \f[I]Radiance\f[R] RGB values or CIE XYZ values.
.RS 3
.RE
.SH 2
Z-buffer Format (.zbf suffix)
.pdfhref O 2 "Z-buffer Format \(.zbf suffix\)"
.pdfhref M "z-buffer-format-.zbf-suffix"
.LP
The Z-buffer format used in \f[I]Radiance\f[R] hardly qualifies as a
format at all.
It is in fact a straight copy on the 4-byte machine floating point
values of each pixel in standard scanline order.
There is no information header or resolution string that would make the
file independently useful.
This is usually OK, because Z-buffer files are almost always created and
used in conjunction with a picture file, which has this identifying
information.
.PP
The major shortcoming of this format is that the machine representation
and byte ordering is not always the same from one system to another,
which limits the portability of Z-buffer files.
Since they are primarily used for interpolating animation frames, and
this usually occurs on networks with similar architectures, there is
usually no problem.
Also, since the adoption of IEEE standard floating-point calculations,
different machine representations are becoming quite rare.
[TBD: is this necessary at this point?]
.SH 3
\f[BI]Radiance\f[B] programs
.pdfhref O 3 "Radiance programs"
.pdfhref M "radiance-programs-6"
.LP
Table 7 shows the programs that read and write \f[I]Radiance\f[R]
Z-buffer files.
The pvalue program may be used to convert Z-buffer files to
\f[I]Radiance\f[R] pictures for the purpose of visualizing the values
using falsecolor.
For example, the following command converts the Z-buffer file
\f[C]frame110.zbf\f[R] associated with the picture
\f[C]frame110.hdr\f[R] to a viewable image:
.IP
.nf
\f[C]
% pvalue -h \[ga]getinfo -d < frame110.hdr\[ga] -r -b -df
frame110.zbf | falsecolor -m 1 -s 40 -l Meters >
frame110z.hdr
\f[]
.fi
.LP
The \f[B]getinfo\f[R] program appearing in back-quotes was used to get
the dimensions associated with the Z-buffer from its corresponding
picture file.
.PP
.na
.TS
delim(@@) tab(	);
l l l l.
T{
Program
T}	T{
Read
T}	T{
Write
T}	T{
Function
T}
_
T{
\f[B]pinterp\f[R]
T}	T{
X
T}	T{
X
T}	T{
Interpolate/extrapolate picture views
T}
T{
\f[B]pvalue\f[R]
T}	T{
X
T}	T{
X
T}	T{
Convert picture to/from simpler formats
T}
T{
\f[B]rad\f[R]
T}	T{
T}	T{
X
T}	T{
Render Radiance scene
T}
T{
\f[B]ranimate\f[R]
T}	T{
X
T}	T{
X
T}	T{
Animate Radiance scene
T}
T{
\f[B]rpict\f[R]
T}	T{
T}	T{
X
T}	T{
Batch rendering program
T}
T{
\f[B]rtrace\f[R]
T}	T{
T}	T{
X
T}	T{
Customizable ray-tracer
T}
.TE
.ad
.LP
\f[B]Table 7.\f[R] \f[I]Radiance\f[R] programs that read and write
Z-buffer files.
.SH 3
\f[BI]Radiance\f[B] C Library
.pdfhref O 3 "Radiance C Library"
.pdfhref M "radiance-c-library-6"
.LP
There are no library functions devoted to reading and writing Z-buffer
files in particular.
The normal method is to read and write Z-buffer scanlines with the
standard \f[C]fread\f[R] and \f[C]fwrite\f[R] library functions using an
appropriate float array.
.SH 2
Ambient File Format (.amb suffix)
.pdfhref O 2 "Ambient File Format \(.amb suffix\)"
.pdfhref M "ambient-file-format-.amb-suffix"
.LP
\f[I]Radiance\f[R] can store its diffuse interreflection cache in an
\f[I]ambient file\f[R] for reuse by other processes.
This file is in a system-independent binary format, similar to an octree
or picture file, with an information header that can be read using
\f[B]getinfo\f[R].
Following the header, there is a magic number specific to this file
type, then the ambient value records themselves in encoded form.
.SH 4
Information Header
.pdfhref O 4 "Information Header"
.pdfhref M "information-header-2"
.LP
The information header begins with the usual \[lq]#?RADIANCE\[rq]
identifier, followed by the originating program and the ambient
calculation parameters (and octree name).
After this is the line:
.IP
.nf
\f[C]
 FORMAT=Radiance_ambval
\f[]
.fi
.LP
This identifies the general file type, followed by an empty line ending
the header.
As with most information headers, this exact sequence need not be
followed, so long as there is no inconsistent \f[C]FORMAT\f[R] setting.
.SH 4
Magic Number
.pdfhref O 4 "Magic Number"
.pdfhref M "magic-number"
.LP
Following the information header is the two-byte magic number, which for
the current ambient file format is 557.
This number may change later should the file format be altered in
incompatible ways.
.SH 4
Ambient Value Records
.pdfhref O 4 "Ambient Value Records"
.pdfhref M "ambient-value-records"
.LP
Ambient values are written to the file in no particular order.
Each diffuse interreflection value in \f[I]Radiance\f[R] has the
following members:
.IP "\f[B]Level\f[R]" 3
The number of reflections between the primary (eye) ray and this
surface.
A value with fewer reflections may be used in place of one with more,
but not the other way around.
.RS 3
.RE
.IP "\f[B]Weight\f[R]" 3
The weighting value associated with this ray or ambient value.
Similar to the level to avoid using inappropriate values from the cache.
.RS 3
.RE
.IP "\f[B]Position\f[R]" 3
The origin point of this interreflection calculation.
.RS 3
.RE
.IP "\f[B]Direction\f[R]" 3
The surface normal indicating the zenith of the sample hemisphere for
this value.
.RS 3
.RE
.IP "\f[B]Value\f[R]" 3
The calculated indirect irradiance at this point, in w/m\*{2\*} (RGB
color).
.RS 3
.RE
.IP "\f[B]Radius\f[R]" 3
The cosine-weighted, harmonic mean distance to other surfaces visible
from this point, used to decide point spacing.
.RS 3
.RE
.IP "\f[B]Posgradient\f[R]" 3
The position-based gradient vector, indicating how brightness changes as
a function of position in the sample plane.
.RS 3
.RE
.IP "\f[B]Dirgradient\f[R]" 3
The direction-based gradient vector, indicating how brightness changes
as a function of surface orientation.
.RS 3
.RE
.LP
The members are stored one after the other in the above order using
system-independent binary representations.
The Level member takes 1 byte, Weight takes 5, Position takes 15,
Direction another 15, Value is 4 bytes (using the same color format as
\f[I]Radiance\f[R] pictures), Radius takes 5 bytes, and Posgradient and
Dirgradient each take 15 bytes, for a total size of 75 bytes per record.
.SH 3
\f[BI]Radiance\f[B] Programs
.pdfhref O 3 "Radiance Programs"
.pdfhref M "radiance-programs-7"
.LP
Table 8 shows \f[I]Radiance\f[R] programs that read and write ambient
files.
The program \f[B]lookamb\f[R] is especially useful for examining the
contents of ambient files and debugging problems in the calculation.
.PP
.na
.TS
delim(@@) tab(	);
l l l l.
T{
Program
T}	T{
Read
T}	T{
Write
T}	T{
Function
T}
_
T{
\f[B]getinfo\f[R]
T}	T{
X
T}	T{
T}	T{
Print information header from binary file
T}
T{
\f[B]lookamb\f[R]
T}	T{
X
T}	T{
X
T}	T{
Convert \f[I]Radiance\f[R] ambient file
T}
T{
\f[B]rad\f[R]
T}	T{
X
T}	T{
X
T}	T{
Render \f[I]Radiance\f[R] scene
T}
T{
\f[B]rpict\f[R]
T}	T{
X
T}	T{
X
T}	T{
Batch rendering program
T}
T{
\f[B]rpiece\f[R]
T}	T{
X
T}	T{
X
T}	T{
Parallel batch rendering program
T}
T{
\f[B]rtrace\f[R]
T}	T{
X
T}	T{
X
T}	T{
Customizable ray-tracer
T}
T{
\f[B]rview\f[R]
T}	T{
X
T}	T{
X
T}	T{
Interactive renderer
T}
.TE
.ad
.LP
\f[B]Table 8.\f[R] Programs in the \f[I]Radiance\f[R] distribution that
read and write ambient files.
.SH 3
\f[BI]Radiance\f[B] C Library
.pdfhref O 3 "Radiance C Library"
.pdfhref M "radiance-c-library-7"
.LP
The \f[C]src/rt/ambient.h\f[R] file contains definitions of the
\f[C]AMBVAL\f[R] structure and certain details of the ambient file
format.
The \f[C]src/rt/ambio.c\f[R] module contains the specialized routines
for reading and writing ambient files, and these functions in turn
access routines in \f[C]src/common/portio.c\f[R] for reading and writing
portable binary data.
The information header is handled by the routines in
s\f[C]rc/common/header.c\f[R], similar to the method described for
\f[I]Radiance\f[R] picture files.
Here are the main calls from \f[C]src/rt/ambio.c\f[R]:
.IP "\f[B]\f[CB]putambmagic(FILE *fp);\f[B]\f[R]" 3
Put out the appropriate two-byte magic number for a \f[I]Radiance\f[R]
ambient file to the stream \f[C]fp\f[R].
.RS 3
.RE
.IP "\f[B]\f[CB]int hasambmagic(FILE *fp);\f[B]\f[R]" 3
Read the next two bytes from the stream \f[C]fp\f[R] and return non-zero
if they match an ambient file\[cq]s magic number.
.RS 3
.RE
.IP "\f[B]\f[CB]int writeambval(AMBVAL *av, FILE *fp);\f[B]\f[R]" 3
Write out the ambient value structure \f[C]av\f[R] to the stream
\f[C]fp\f[R], returning -1 if a file error occurred, or 0 normally.
.RS 3
.RE
.IP "\f[B]\f[CB]int readambval(AMBVAL *av, FILE *fp);\f[B]\f[R]" 3
Read in the next ambient value structure from the stream \f[C]fp\f[R]
and put the result in \f[C]av\f[R].
Return 1 if the read was successful, 0 if the end of file was reached or
there was an error.
The \f[C]ambvalOK\f[R] function is used to check the consistency of the
value read.
.RS 3
.RE
.IP "\f[B]\f[CB]int ambvalOK(AMBVAL *av);\f[B]\f[R]" 3
Return non-zero if the member values of the \f[C]av\f[R] structure are
not too outlandish.
This is handy as insurance against a corrupted ambient file.
.RS 3
.RE
.SH 2
Conclusion
.pdfhref O 2 "Conclusion"
.pdfhref M "conclusion"
.LP
We have described the main file formats native to \f[I]Radiance\f[R] and
shown how even the binary formats can be reliably shared in
heterogeneous computing environments.
This corresponds to one of the basic philosophies of UNIX software,
which is system independence.
A more detailed understanding of the formats may still require some use
of binary dump programs and delving into the \f[I]Radiance\f[R] source
code.
