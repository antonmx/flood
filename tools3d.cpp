
#include "tools3d.h"
#include <string.h>

using namespace blitz;
using namespace std;







clargs::clargs(int argc, char *argv[]) :
  proc_table(),
  save_table(),
  table("3D advanced flood fill.",
        "Advanced flood fill algorithm implemented for the 3D volume.")
{

  proc_table
    .add(poptmx::NOTE, "Processing options:")
    .add(poptmx::OPTION, &g.start, 's', "start",
         "Starting point(s) of the fill procedure.", "")
    .add(poptmx::OPTION, &g.stop, 'S', "stop",
         "Stop sphere(s).", "Three coordinates and radius.")
    .add(poptmx::OPTION, &g.radius, 'r', "test-radius",
         "Radius of the test sphere.", "")
    .add(poptmx::OPTION, &g.radiusM, 'R', "mark-radius",
         "Radius of the fill sphere.", "")
    .add(poptmx::OPTION, &g.minval, 'm', "minval",
         "Minimum value.", "")
    .add(poptmx::OPTION, &g.maxval, 'M', "maxval",
         "Maximum value.", "")
    .add(poptmx::OPTION, &g.minggrad, 'g', "mingrad",
         "Minimum absolute value of the gradient.", "")
    .add(poptmx::OPTION, &g.maxggrad, 'G', "maxgrad",
         "Maximum absolute value of the gradient.", "")    ;
//    .add(poptmx::OPTION, &nextdim, 'n', "next",
//         "Dimension for finding neigbours.",
//         "Can be 1 - lines (6 neigbours), 2 - planes (14) and 3 - cube (26)")


  color_table
    .add(poptmx::OPTION, &g.color, 'c', "color",
         "Color of the fill volume.",
         "Instead of masking the volume, paints it with this color.")  ;

  save_table
    .add(poptmx::NOTE, "Saving options:")
    .add(color_table)
    .add(poptmx::OPTION, &g.out_filled, 'f', "outfilled",
         "Prefix to output filled.", "", g.out_filled)
    .add(poptmx::OPTION, &g.out_inverted, 'e', "outinvert",
         "Prefix to output not filled.", "", g.out_inverted)
    .add(poptmx::OPTION, &g.out_mask, 'b', "outmask",
         "Prefix to output mask.", "Outputs bit mask.", g.out_mask);


  table
    .add(poptmx::NOTE, "ARGUMENTS:")
    .add(poptmx::ARGUMENT, &g.inlist, "list", "List of the input images.", "")

    .add(poptmx::NOTE, "OPTIONS:")
    .add(poptmx::OPTION, &g.interactive, 'i', "interactive",
         "Accepts commands from the standard input.", "")
    .add(proc_table)
    .add(save_table)
    .add_standard_options(&g.beverbose);


  if ( ! table.parse(argc,argv) )
    exit(0);
  if ( ! table.count() ) {
    table.usage();
    exit(0);
  }

  g.command = table.name();

  // <list> : required argument.
  if ( ! table.count(&g.inlist) )
    exit_on_error(g.command, "Missing required argument: "+table.desc(&g.inlist)+".");


  if (g.interactive)
    return;

  check_proc(table);
  check_save(table);

}


void clargs::check_proc( poptmx::OptionTable & tab ) {

  if ( g.radius < 0 )
    throw_error(g.command, "Impossible parameter value : " + tab.desc(&g.radius) + "."
                             " Cannot be negative.");
  if ( ! tab.count(&g.radiusM) )
    g.radiusM = g.radius;
  if ( tab.count(&g.minval) && ( g.minval < 0 || g.minval > 256 ) )
    throw_error(g.command, "Impossible parameter value : " + tab.desc(&g.minval) + "."
                           " Must be in the range [0,256].");
  if ( tab.count(&g.maxval) && ( g.maxval < 0 || g.maxval > 256 ) )
    throw_error(g.command, "Impossible parameter value : " + tab.desc(&g.maxval) + "."
                           " Must be in the range [0,256].");
  if ( tab.count(&g.minggrad) && ( g.minggrad < 0 || g.minggrad > 256 ) )
    throw_error(g.command, "Impossible parameter value : " + tab.desc(&g.minggrad) + "."
                           " Must be in the range [0,256].");
  if ( tab.count(&g.minggrad) )
    g.minggrad *= g.minggrad;
  if ( tab.count(&g.maxggrad) && ( g.maxggrad < 0 || g.maxggrad > 256 ) )
    throw_error(g.command, "Impossible parameter value : " + tab.desc(&g.maxggrad) + "."
                           " Must be in the range [0,256].");
  if ( tab.count(&g.maxggrad) )
    g.maxggrad *= g.maxggrad;

    // point : required argument.
  if ( ! tab.count(&g.start) )
    throw_error(g.command, "Missing required argument: "+tab.desc(&g.start)+".");

//  if ( nextdim < 1 || nextdim > 3 )
//    throw_error(g.command, "Impossible parameter value : " + tab.desc(&nextdim) + ".");


}

void clargs::check_save( poptmx::OptionTable & tab ) {

    if ( ! ( tab.count(& g.out_filled) +
           tab.count(& g.out_inverted) +
           tab.count(& g.out_mask) ) )
      throw_error(g.command, "At least one of the two following arguments is required: "
                             +tab.desc(&g.out_mask)+ ", "
                             +tab.desc(&g.out_filled)+ ", "
                             +tab.desc(&g.out_inverted)+ ".");

}




void sig1handler(int signum) {
  g.stopProcessing=true;
}

void sig2handler(int signum) {

  Point3D pnt;
  char line[256];

  FILE *funcf = fopen( ".sliceme.txt", "r" );

  if ( ! funcf ||
       1 != fscanf ( funcf,  "%254[^\n]\n", line) ||
       1 != _conversion (&pnt, line) ||
       ! pnt.inVolume(g.wvol.shape()) )
    pnt = g.cross_point;

  if (funcf)
    fclose(funcf);

  Map8U yz( g.wvol.shape()(2), g.wvol.shape()(1) );
  yz = g.wvol( pnt.x(), Range::all(), Range::all() ).transpose(secondDim, firstDim) ;
  SaveImage(".yz.tif", yz);
  Map8U xz( g.wvol.shape()(2), g.wvol.shape()(0) );
  xz = g.wvol( Range::all(), pnt.y(), Range::all() ).transpose(secondDim, firstDim) ;
  SaveImage(".xz.tif", xz);
  Map8U xy = g.wvol( Range::all(), Range::all(), pnt.z() ).copy();
  SaveImage(".xy.tif", xy);

  printf("Written slices through (%i, %i, %i) point."
         " Files are named \".xy.tif\", \".xz.tif\" and \".yz.tif\".\n",
         pnt.x(), pnt.y(), pnt.z() );

}


void * in_zeroing_thread (void * _thread_args) {

  int threadnum = (long int) _thread_args; // dirty hack

  const long int chunk = 1 + ( g.ivol.size() / g.run_threads );
  const long int start = threadnum * chunk;

  uint8_t * wdat = g.wvol.data() + start;
  const uint8_t * wend = g.wvol.data() + min( g.wvol.size(), start + chunk );
  while ( wdat != wend )
    *wdat++ = 0;


}





void * in_wiping_thread (void * _thread_args) {

  int threadnum = (long int) _thread_args; // dirty hack

  const long int chunk = 1 + ( g.ivol.size() / g.run_threads );
  const long int start = threadnum * chunk;
  const uint8_t mask = ( ~ISBAD ) & ( ~ISGOOD );

  uint8_t * wdat = g.wvol.data() + start;
  const uint8_t * wend = g.wvol.data() + min( g.wvol.size(), start + chunk );
  while ( wdat != wend )
    *wdat++ &= mask;

}


void * in_apply_thread (void * _thread_args) {

  ApplyArgs *  args = (ApplyArgs*) _thread_args;
  if (!args)
    throw_error("apply mask thread", "Inappropriate thread function arguments.");

  const long int chunk = 1 + ( g.ivol.size() / g.run_threads );
  const long int start = args->threadnum * chunk;

  uint8_t * idat = g.ivol.data() + start,
          * wdat = g.wvol.data() + start;
  const uint8_t * iend = g.ivol.data() + min( g.ivol.size(), start + chunk );
  while ( idat != iend ) {
    if ( ( ! args->invert  &&  ( ! ( *wdat & FILLED ) ) )  ||
         (   args->invert  &&  (   *wdat & FILLED ) )  )
      *idat = g.color;
    wdat++;
    idat++;
  }

}


pthread_mutex_t ReadWriteDistributor::iolock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t ReadWriteDistributor::proglock = PTHREAD_MUTEX_INITIALIZER;



void * in_read_thread (void * _thread_args) {

  ReadWriteDistributor *  dist = (ReadWriteDistributor*) _thread_args;
  if (!dist)
    throw_error("read thread", "Inappropriate thread function arguments.");

  Shape2D imgsh = Shape2D(g.ivol.extent(secondDim) , g.ivol.extent(thirdDim) ) ;
  Map8U img(imgsh);

  long int idx;
  while ( dist->distribute(&idx) ) {
    ReadImage(g.inlist[idx], img, imgsh);
    g.ivol( idx, Range::all(), Range::all() ) = img;
    dist->updateProg();
  }

}





void * in_write_thread (void * _thread_args) {

  ReadWriteDistributor *  dist = (ReadWriteDistributor*) _thread_args;
  if (!dist)
    throw_error("write thread", "Inappropriate thread function arguments.");

  const Shape2D imgsh(g.ivol.extent(secondDim) , g.ivol.extent(thirdDim));
  Map8U img(imgsh), wmg(imgsh);
  Path outPath;

  long int idx;
  while ( dist->distribute(&idx) ) {

    wmg = g.wvol(idx, Range::all(), Range::all());

    if ( ! g.out_filled.empty() ) {
      for ( long int sh0 = 0 ; sh0 < imgsh(0) ; sh0++)
        for ( long int sh1 = 0 ; sh1 < imgsh(1) ; sh1++)
          img(sh0, sh1) =  ( wmg(sh0, sh1) & FILLED ) ? g.ivol( idx, sh0, sh1 ) : g.color ;
      outPath = g.out_filled + g.inlist[idx].name();
      SaveImage(outPath, img);
    }

    if ( ! g.out_inverted.empty() ) {
      for ( long int sh0 = 0 ; sh0 < imgsh(0) ; sh0++)
        for ( long int sh1 = 0 ; sh1 < imgsh(1) ; sh1++)
          img(sh0, sh1) =  ( wmg(sh0, sh1) & FILLED ) ? g.color : g.ivol( idx, sh0, sh1 ) ;
      outPath = g.out_inverted + g.inlist[idx].name();
      SaveImage(outPath, img);
    }

    if ( ! g.out_mask.empty() ) {
      outPath = g.out_mask + g.inlist[idx].name();
      SaveImage(outPath, wmg);
    }

    dist->updateProg();

  }


}













#ifndef NUL
#define NUL '\0'
#endif

static char const nomem[] = "no memory for %lu byte allocation\n";

static bool
copy_raw_string(char ** dest_p, char ** src_p);

static bool
copy_cooked_string(char ** dest_p, char ** src_p);

static inline void *
Xmalloc(size_t sz)
{
    void * res = malloc(sz);
    if (res == NULL) {
        fprintf(stderr, nomem, sz);
        exit(EXIT_FAILURE);
    }
    return res;
}

static inline void *
Xrealloc(void * ptr, size_t sz)
{
    void * res = realloc(ptr, sz);
    if (res == NULL) {
        fprintf(stderr, nomem, sz);
        exit(EXIT_FAILURE);
    }
    return res;
}

bool
string_to_argv(char const * str, int * argc_p, char *** argv_p)
{
    int     argc = 0;
    int     act  = 10;
    char ** res  = (char**) Xmalloc( sizeof(char *) * 10 );
    char ** argv = res;
    char *  scan;
    char *  dest;
    bool err;

    while (isspace((unsigned char)*str))  str++;
    str = scan = strdup(str);

    for (;;) {
        while (isspace((unsigned char)*scan))  scan++;
        if (*scan == NUL)
            break;

        if (++argc >= act) {
            act += act / 2;
            res  = (char**) Xrealloc(res, act * sizeof(char *));
            argv = res + (argc - 1);
        }

        *(argv++) = dest = scan;

        for (;;) {
            char ch = *(scan++);
            switch (ch) {
            case NUL:
                goto done;

            case '\\':
                if ( (*(dest++) = *(scan++)) == NUL)
                    goto done;
                break;

            case '\'':
                err = copy_raw_string(&dest, &scan);
                if (err)
                    goto error_leave;
                break;

            case '"':
                err = copy_cooked_string(&dest, &scan);
                if (err)
                    goto error_leave;
                break;

            case ' ':
            case '\t':
            case '\n':
            case '\f':
            case '\r':
            case '\v':
            case '\b':
                goto token_done;

            default:
                *(dest++) = ch;
            }
        }

    token_done:
        *dest = NUL;
    }

done:

    *argv_p = res;
    *argc_p = argc;
    *argv   = NULL;
    if (argc == 0)
        free((void *)str);

    return false;

error_leave:

    free(res);
    free((void *)str);
    return err;
}

static bool
copy_raw_string(char ** dest_p, char ** src_p)
{
    for (;;) {
        char ch = *((*src_p)++);

        switch (ch) {
        case NUL: return true;
        case '\'':
            *(*dest_p) = NUL;
            return false;

        case '\\':
            ch = *((*src_p)++);
            switch (ch) {
            case NUL:
                return true;

            default:
                /*
                 * unknown/invalid escape.  Copy escape character.
                 */
                *((*dest_p)++) = '\\';
                break;

            case '\\':
            case '\'':
                break;
            }
            /* FALLTHROUGH */

        default:
            *((*dest_p)++) = ch;
            break;
        }
    }
}

static char
escape_convt(char ** src_p)
{
    char ch = *((*src_p)++);

    /*
     *  Escape character is always eaten.  The next character is sometimes
     *  treated specially.
     */
    switch (ch) {
    case 'a': ch = '\a'; break;
    case 'b': ch = '\b'; break;
    case 't': ch = '\t'; break;
    case 'n': ch = '\n'; break;
    case 'v': ch = '\v'; break;
    case 'f': ch = '\f'; break;
    case 'r': ch = '\r'; break;
    }

    return ch;
}


static bool
copy_cooked_string(char ** dest_p, char ** src_p)
{
    for (;;) {
        char ch = *((*src_p)++);
        switch (ch) {
        case NUL: return true;
        case '"':
            *(*dest_p) = NUL;
            return false;

        case '\\':
            ch = escape_convt(src_p);
            if (ch == NUL)
                return true;
            /* FALLTHROUGH */

        default:
            *((*dest_p)++) = ch;
            break;
        }
    }
}
