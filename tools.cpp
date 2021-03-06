/******************************************************************************
 *   Copyright (C) 2007 by Anton Maksimenko                                   *
 *   antonmx@gmail.com
 *                                                                            *
 *   This program is free software; you can redistribute it and/or modify     *
 *   it under the terms of the GNU General Public License as published by     *
 *   the Free Software Foundation; either version 2 of the License, or        *
 *   (at your option) any later version.                                      *
 *                                                                            *
 *   This program is distributed in the hope that it will be useful,          *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *   GNU General Public License for more details.                             *
 *                                                                            *
 *   You should have received a copy of the GNU General Public License        *
 *   along with this program; if not, write to the                            *
 *   Free Software Foundation, Inc.,                                          *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.                *
 ******************************************************************************/


///
/// @file
/// @author antonmx <antonmx@gmail.com>
/// @date   Mon Jul 21 10:09:31 2008
///
/// @brief %Accessoires for the 3D flood fill algorithm.
///


#include "tools.h"

#include <unistd.h>
#include <vector>
#include <tiffio.h>
#include <fcntl.h>
#include <string.h>

#define MAGICKCORE_QUANTUM_DEPTH 16
#define MAGICKCORE_HDRI_ENABLE false

#ifdef _WIN32
#  define STATIC_MAGICK
#  define MAGICK_STATIC_LINK
#endif
#include<Magick++.h>

using namespace std;
using namespace blitz;






const Time::time_point startTV = Time::now();
Time::time_point prevTV = startTV;

void prdn( int a ) {
  const Time::time_point nowTV = Time::now();
  double start_elapsed = chrono::nanoseconds(nowTV - startTV).count()/1000000000.0;
  double prev_elapsed  = chrono::nanoseconds(nowTV - prevTV ).count()/1000000000.0;
  printf("DONE %i:  %f  %f\n", a, prev_elapsed, start_elapsed);
  fflush(stdout);
  prevTV=nowTV;
}


void prdn( const std::string & msg ) {
  const Time::time_point nowTV = Time::now();
  double start_elapsed = chrono::nanoseconds(nowTV - startTV).count()/1000000000.0;
  double prev_elapsed  = chrono::nanoseconds(nowTV - prevTV ).count()/1000000000.0;
  printf("DONE \"%s\":  %f  %f\n", msg.c_str(), prev_elapsed, start_elapsed);
  fflush(stdout);
  prevTV=nowTV;
}

void prdnc( int a ) {
  printf("%i ", a);
  fflush(stdout);
}


void prdnc( const std::string & msg ) {
  printf("%s ", msg.c_str());
  fflush(stdout);
}




/// \brief Prints the error into the standard error stream.
void
FloodErr::report() const {
  switch ( terr ) {
  case WARN:
    cerr << "WARNING!";
    break;
  case ERR:
    cerr << "ERROR!";
    break;
  case FATAL:
    cerr << "FATAL ERROR!";
    break;
  }
  cerr << " In module \'" << module << "\'. " << message << endl;
}


FloodErr::ErrTp
FloodErr::type() const {
  return terr;
}


void
throw_error(const string & mod, const string & msg) {
  FloodErr err(mod, msg, FloodErr::ERR);
  err.report();
  throw err;
}

void
throw_bug(const string & mod){
  throw_error(mod, "Must never happen. This is the bug, please report to developers.");
}

FloodErr
warn(const string & mod, const string & msg){
  FloodErr err(mod, msg, FloodErr::WARN);
  err.report();
  return err;
}

void
exit_on_error(const string & mod, const string & msg){
  FloodErr err(mod, msg, FloodErr::FATAL);
  err.report();
  cerr << "!!! EXITING ON ERROR !!!" << endl;
  exit(1);
}



#ifdef _WIN32
const string Path::DIRSEPARATOR = "\\";
#else
const string Path::DIRSEPARATOR = "/";
#endif

const Path Path::emptypath = Path();


/// Constructs the error which reports the unextractable element.
///
/// @param element The element whose extraction has failed.
///
void
Path::throw_unextracted(const string & element) const {
  throw_error
  ("class Path",
   "Could not extract " + element + " from path \"" + *this + "\".");
}

string
Path::drive () const {
#ifdef _WIN32
  char _drive[FILENAME_MAX];
  if ( _splitpath_s( c_str(), _drive, FILENAME_MAX, 0,0,0,0,0,0 ) )
  throw_unextracted("drive");
  return _drive;
#else
  return "";
#endif
}

string
Path::dir () const {
#ifdef _WIN32
  char _dir[FILENAME_MAX];
  if ( _splitpath_s( c_str(),0,0,_dir,FILENAME_MAX,0,0,0,0 ) )
  throw_unextracted("directory");
  return drive() + _dir;
#else
  string::size_type idx=this->rfind("/");
  return  idx == string::npos  ?  string("")  :  string(*this, 0, idx+1);
#endif
}

string
Path::name () const {
#ifdef _WIN32
  return title() + extension();
#else
  size_t found = this->rfind('/');
  if (found==string::npos)
    return *this;
  else
    return this->substr(found+1);
#endif
};


string
Path::title () const {
#ifdef _WIN32
  char _title[FILENAME_MAX];
  if ( _splitpath_s( c_str(), 0,0,0,0, _title, FILENAME_MAX, 0,0 ) )
  throw_unextracted("title");
  return _title;
#else
  string::size_type idx = name().rfind(".");
  if (idx==string::npos || idx == 0 )
  return name();
  else
  return name().substr(0, idx);
#endif
}

string
Path::dtitle () const {
  return dir() + title();
}


string
Path::extension () const {
#ifdef _WIN32
  char _ext[FILENAME_MAX];
  if ( _splitpath_s( c_str(), 0,0,0,0,0,0, _ext, FILENAME_MAX ) )
  throw_unextracted("extension");
  return _ext;
#else
  string::size_type
  dotidx = this->rfind("."),
  diridx = this->rfind("/");
  if ( dotidx != string::npos && ( diridx == string::npos || diridx+1 < dotidx ) )
  return this->substr(dotidx);
  else
  return "";
#endif
}



bool
Path::isdir() const {
  return ! empty() && name().empty();
}

Path &
Path::bedir() {
  if ( ! empty() && ! isdir() )
  *this += DIRSEPARATOR;
  return *this;
}


bool
Path::isabsolute() const {
#ifdef _WIN32
  return ! drive().empty();
#else
  return ! empty() && (*this)[0] == '/';
#endif
}


const Path Path::home() {
  char * hm;
  #ifdef _WIN32
  //return getenv("USERPROFILE");
  hm = getenv("APPDATA");
  #else
  hm = getenv("HOME");
  #endif
  return Path(hm).bedir();
}


Path
upgrade(const Path & path, const string & addthis) {
  return path.dir() + Path(addthis) + path.name();
}


Path
cdpath(const Path & dir, const Path & file){
  Path ddir = dir;
  ddir.bedir();

  if ( file.empty() )
  return ddir;
  else if ( file.isabsolute() || ddir.empty() )
  return file;
  else
      // Windows version is likely to fail here. See Path operator+.
  return ddir + file;
}


string
type_desc (Path*) {
  return "PATH";
}

int
_conversion (Path* _val, const string & in) {
  *_val = Path(in);
  return 1;
}










long
nof_threads(long _threads) {
  if (_threads)
    return _threads;

#ifdef _WIN32
#ifndef _SC_NPROCESSORS_ONLN
  SYSTEM_INFO info;
  GetSystemInfo(&info);
#define sysconf(a) info.dwNumberOfProcessors
#define _SC_NPROCESSORS_ONLN
#endif
#endif

  long nProcessorsOnline = sysconf(_SC_NPROCESSORS_ONLN);
  if (nProcessorsOnline == -1) {
    warn ("thread number",
          "Unable to read online processor count.");
    return 1;
  } else {
    return nProcessorsOnline;
  }
}






/// \cond
struct ToLower{ char operator() (char c) const { return tolower(c); } };
struct ToUpper{ char operator() (char c) const { return toupper(c); } };
/// \endcond

string
upper(string str){
  transform (str.begin(), str.end(), str.begin(), ToUpper());
  return str;
}

string
lower(string str){
  transform (str.begin(), str.end(), str.begin(), ToLower());
  return str;
}


string
toString(const string fmt, ...){
  va_list ap;
  va_start(ap, fmt);
  const int MAXLINE = 4096;
  char buf[MAXLINE];
  vsnprintf(buf, MAXLINE, fmt.c_str(), ap);
  va_end(ap);
  return string(buf);
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







/// \brief Initializing constructor.
///
/// @param _showme Tells if the progress bar should be shown.
/// @param _message The description of the progress.
/// @param _steps Number of steps in the progress. If 0 then unknown.
///
ProgressBar::ProgressBar(bool _showme, const string & _message, uint64_t _steps) :
  showme(_showme),
  message(_message),
  steps(_steps)
{

  step = 0;
  waswidth = 0;
  reservedChs = 0;

  int nums = toString(steps).length();
  reservedChs = 14 + 2*nums;

  fmt = steps ?
    "%" + toString(nums) + "u of " + toString(steps) + " [%s] %4s" :
    string( "progress: %u" );

  progPrevTV = Time::now();

  if ( showme ) {;
    cout << "Starting process";
    if (steps) cout << " (" + toString(steps) + " steps)";
    cout << ": " << message << "." << endl;
    fflush(stdout);
  }


}


string ProgressBar::print_line() {

  int progln = getwidth() - reservedChs;
  if ( progln <= 3 )  return ""; // if we have a very narrow terminal

  if ( steps && step >= steps ) {
    done();
    return "";
  }

  string outS;
  if (steps) {
    string eqs = string(progln*step/steps, '=') + string(progln, ' ') ;
    eqs.erase(progln);
    string prc = toString("%5.1f%% ", (100.0*step)/steps);
    outS = toString(fmt, step, eqs.c_str(), prc.c_str() );
  } else {
    outS = toString(fmt, step);
  }

  return outS;

}



/// \brief Updates the progress bar.
///
/// @param curstep Current step. Advances +1 if zero.
///
void
ProgressBar::update(uint64_t curstep){

  step = (curstep ? curstep : step) + 1;

  if ( !showme || !reservedChs ) return; // Uninitialized progress bar.

  Time::time_point progNowTV = Time::now();
  if ( progNowTV - progPrevTV > Time::duration(std::chrono::milliseconds(100)) )
    progPrevTV = progNowTV;
  else
    return;

  string outS = print_line();
  cout << string(waswidth+1, '\b') << outS ;
  fflush(stdout);
  waswidth = outS.length();

}


void
ProgressBar::done(){

  if ( !showme || ! reservedChs ) return;

  int progln = getwidth() - reservedChs;
  if ( progln < 0 )  progln = 0; // if we have a very narrow terminal
  string eqs(progln, '=');

  cout << string(waswidth+1, '\r')
     << ( steps ?
      toString(fmt, steps, eqs.c_str(), "DONE. ") :
      toString(fmt, step) + " steps. DONE." )
     << endl
     << "Successfully finished " << message << "." << endl;
  fflush(stdout);

  reservedChs = 0;

}


#ifdef _WIN32
#  include <windows.h>
#else
#  include<termios.h>
#  include<sys/ioctl.h>
#endif

int
ProgressBar::getwidth(){
#ifdef _WIN32
  CONSOLE_SCREEN_BUFFER_INFO info;
  HANDLE fd = GetStdHandle(STD_OUTPUT_HANDLE);
  if (fd == INVALID_HANDLE_VALUE)
    return 0;
  return GetConsoleScreenBufferInfo (fd, &info) ? info.dwSize.X - 1 : 0;
  //return (info.srWindow.Right - info.srWindow.Left + 1);
#else
  winsize size;
  return ( ioctl (STDOUT_FILENO, TIOCGWINSZ, &size ) < 0 ) ?  0 : size.ws_col - 1;
#endif
}

















Shape2D
ImageSizes(const Path & filename){
  Magick::Image imag;
  try {
    imag.ping(filename);
  }
  catch ( Magick::WarningCoder err ) {}
  catch ( Magick::Exception & error) {
    throw_error("get image size", "Could not read image file\""+filename+"\"."
                        " Caught Magick++ exception: \""+error.what()+"\".");
  }
  return Shape2D( imag.rows(), imag.columns() );
}




void
BadShape(const Path & filename, const Shape2D & shp){
  Shape2D ashp = ImageSizes(filename);
  if ( ashp != shp )
    throw_error("load image", "The shape of the image"
                "\"" + filename + "\"  (" + toString(ashp) + ") is not equal"
                " to the requested shape (" + toString(shp)  + ").");
}




/// Loads an image (lines) using TIFF library.
///
/// @param filename Name of the image
/// @param storage The array to store the image.
/// @param idxs The indexes of the line to read.
///        if empty then reads whole image.
///
static void
ReadImageLine_TIFF (const Path & filename, Map8U & storage,
                    const vector<int> & idxs ) {

  static const string modname = "load 8-bit image tiff";

  // BUG in libtiff
  // On platforms (f.e. CentOS) the TIFFOpen function fails,
  // while TIFFFdOpen works well. On the MS Windows the
  // TIFFFdOpen does not work, while TIFFOpen does.

  int fd=0;
  #ifdef _WIN32
  TIFF *tif = TIFFOpen(filename.c_str(), "r");
  #else
  fd = open (filename.c_str(), O_RDONLY);
  if (fd < 1)
    throw_error(modname,
                "Could not open file \"" + filename + "\" for reading.");
  TIFF *tif = TIFFFdOpen(fd, filename.c_str(), "r");
  #endif

  if( ! tif ) {
    if (fd) close(fd);
    throw FloodErr(modname, "Could not read tif from file\"" + filename + "\".", FloodErr::WARN);
  }

  uint32 width = 0, height = 0;
  uint16 spp = 0, bps = 0, fmt = 0, photo;

  if ( ! TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width) ) {
    TIFFClose(tif);
    throw warn(modname,
               "Image \"" + filename + "\" has undefined width.");
  }

  if ( ! TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height) ) {
    TIFFClose(tif);
    throw warn(modname,
               "Image \"" + filename + "\" has undefined height.");
  }

  if ( ! TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &spp) ) {
    TIFFClose(tif);
    throw_error(modname,
               "Image \"" + filename + "\" has undefined samples per pixel.");
  }
  if ( spp != 1 ) {
    TIFFClose(tif);
    throw_error(modname,
                "Image \"" + filename + "\" is not grayscale.");
  }

  if ( ! TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps) ) {
    TIFFClose(tif);
    throw_error(modname,
                "Image \"" + filename + "\" has undefined bits per sample.");
  }
  if ( bps != 8 ) {
    TIFFClose(tif);
    throw_error(modname,
               "Image \"" + filename + "\" is not 8 bits per sample.");
  }

  if ( ! TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &fmt) )
    fmt = SAMPLEFORMAT_UINT;
  if ( fmt != SAMPLEFORMAT_UINT ) {
    TIFFClose(tif);
    throw_error(modname,
               "Image \"" + filename + "\" has unsupported sample format.");
  }

  const int readheight = idxs.size() ? idxs.size() : height;
  storage.resize(readheight,width);

  tdata_t buf = _TIFFmalloc(TIFFScanlineSize(tif));

  for (uint curidx = 0; curidx < readheight; curidx++) {

    uint32 row = idxs.size() ? idxs[curidx] : curidx;

    if ( row >= height || row < 0 ) {

      warn("load imagelines tiff",
      "The index of the line to be read (" + toString(row) + ")"
      " is outside the image boundaries (" + toString(height) + ").");

      storage(curidx, blitz::Range::all()) = 0;

    } else {

      if ( TIFFReadScanline(tif, buf, row) < 0 ) {
        _TIFFfree(buf);
        TIFFClose(tif);
        throw_error(modname,
                   "Failed to read line " + toString(row) +
                   " in image \"" + filename + "\".");
      }

      storage(curidx, blitz::Range::all()) =  blitz::Array<unsigned char,1> (
        (unsigned char *) buf, blitz::shape(width), blitz::neverDeleteData);

    }

  }

  _TIFFfree(buf);
  TIFFClose(tif);

}


/// Loads an image using ImageMagick library.
///
/// @param filename Name of the image
/// @param storage The array to store the image.
///
inline static void
ReadImage_TIFF (const Path & filename, Map8U & storage) {
  ReadImageLine_TIFF(filename, storage, vector<int>() );
}

void
ReadImage_IM (const Path & filename, Map8U & storage ) {

  Magick::Image imag;
  try { imag.read(filename); }
  catch ( Magick::WarningCoder err ) {}
  catch ( Magick::Exception & error) {
    throw_error("load image IM", "Could not read image file\""+filename+"\"."
                " Caught Magick++ exception: \""+error.what()+"\".");
  }
  if ( imag.type() != Magick::GrayscaleType )
    throw_error("load image IM", "Input image \""+filename+"\" is not grayscale.");
  if ( imag.depth() != 8 )
    throw_error("load image IM", "Input image \""+filename+"\" is not 8-bit depth.");

  const int
    width = imag.columns(),
    hight = imag.rows();
  storage.resize( hight, width );

  for (blitz::MyIndexType curw = 0 ; curw < width ; curw++)
    for (blitz::MyIndexType curh = 0 ; curh < hight ; curh++)
      storage(curh,curw) = Magick::ColorGray(imag.pixelColor(curw, curh)).shade();

}

void
ReadImage (const Path & filename, Map8U & storage ){
  try { ReadImage_TIFF(filename, storage); }
  catch (FloodErr err) {
    if (err.type() != FloodErr::WARN)
      throw;
    ReadImage_IM(filename, storage);
  }
}


void
ReadImage(const Path & filename, Map8U & storage, const Shape2D & shp) {
  BadShape(filename, shp);
  ReadImage(filename, storage);
}




/// Saves image in integer format using ImageMagick library.
///
/// @param filename file to save image into.
/// @param storage array with the image.
///
void
SaveImage (const Path & filename, const Map8U & storage) {

  if ( ! storage.size() ) {
    warn("save image", "Zero-sized array for image \"" + filename + "\" - will not be saved.");
    return;
  }

  const int
    width = storage.columns(),
    hight = storage.rows();

  Map8U _storage;
  if ( storage.isStorageContiguous()  &&  storage.stride() == Shape2D(width,1) )
    _storage.reference(storage);
  else {
    _storage.resize(storage.shape());
    _storage = storage;
  }

  const string extension = lower( filename.extension() );
  if ( extension.empty() || extension == ".tif" || extension == ".tiff" ) {

    TIFF *image = TIFFOpen(filename.c_str(), "w");
    if( ! image )
      throw_error("save tiff image", "Could create file\"" + filename + "\".");

    // We need to set some values for basic tags before we can add any data
    TIFFSetField(image, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(image, TIFFTAG_IMAGELENGTH, hight);
    TIFFSetField(image, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(image, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(image, TIFFTAG_ROWSPERSTRIP, hight);
    TIFFSetField(image, TIFFTAG_SAMPLEFORMAT,SAMPLEFORMAT_UINT);
    TIFFSetField(image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);

    int wret = TIFFWriteRawStrip(image, 0, (void*) _storage.data(), width*hight);
    TIFFClose(image);
    if ( -1 == wret )
      throw_error("save tiff image", "Could not save image to file \"" + filename + "\".");

  } else {

    Magick::Image imag( Magick::Geometry(width, hight), "black" );
    imag.classType(Magick::DirectClass);
    imag.type( Magick::GrayscaleType );
    imag.depth(8);
    imag.magick("TIFF"); // saves to tif if not overwritten by the extension.

    for (blitz::MyIndexType curh = 0 ; curh < hight ; curh++)
      for (blitz::MyIndexType curw = 0 ; curw < width ; curw++)
        imag.pixelColor(curw, curh, Magick::ColorGray(_storage(curh,curw)));

    try { imag.write(filename); }
    catch ( Magick::Exception & error) {
      throw_error("save image IM", "Could not write image file\""+filename+"\"."
        " Caught Magick++ exception: \""+error.what()+"\".");
    }

  }

}


















