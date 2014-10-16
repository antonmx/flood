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
#include <sys/mman.h>
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






static timespec nowTime() {
  timespec nt;
  clock_gettime(CLOCK_MONOTONIC, &nt);
  return nt;
}

const timespec startTV = nowTime();
timespec prevTV = startTV;

inline double sec( const timespec & tm ) {
  return tm.tv_sec + tm.tv_nsec / 1000000000.0; 
}

void prdn( int a ) {
  const timespec nowTV = nowTime();
  double start_elapsed = sec(nowTV) - sec(startTV);
  double prev_elapsed  = sec(nowTV) - sec(prevTV);
  printf("DONE %i:  %f  %f\n", a, prev_elapsed, start_elapsed);
  fflush(stdout);
  prevTV=nowTV;
}


void prdn( const std::string & msg ) {
  timespec nowTV = nowTime();
  double start_elapsed = ( nowTV.tv_nsec - startTV.tv_nsec ) / 1000000000.0;
  double prev_elapsed = ( nowTV.tv_nsec - prevTV.tv_nsec ) / 1000000000.0;
  printf("DONE \"%s\":  %f  %f\n", msg.c_str(), prev_elapsed, start_elapsed);
  fflush(stdout);
  prevTV=nowTV;
}


/// \brief Constructor.
///
/// @param _terr Sets CtasErr::terr
/// @param mod   Sets CtasErr::module
/// @param msg   Sets CtasErr::message
///
FloodErr::FloodErr(ErrTp _terr, const string & mod, const string & msg){
  terr = _terr;
  module = mod;
  message = msg;
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
  FloodErr err(FloodErr::ERR, mod, msg);
  err.report();
  throw err;
}

void
throw_bug(const string & mod){
  throw_error(mod, "Must never happen. This is the bug, please report to developers.");
}

FloodErr
warn(const string & mod, const string & msg){
  FloodErr err(FloodErr::WARN, mod, msg);
  err.report();
  return err;
}

void
exit_on_error(const string & mod, const string & msg){
  FloodErr err(FloodErr::FATAL, mod, msg);
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









/// \brief Initializing constructor.
///
/// @param _showme Tells if the progress bar should be shown.
/// @param _message The description of the progress.
/// @param _steps Number of steps in the progress. If 0 then unknown.
///
ProgressBar::ProgressBar(bool _showme, const string & _message, long int _steps) :
  showme(_showme),
  message(_message),
  steps(_steps)
{

  if ( ! showme ) return;

  step = 0;
  waswidth = 0;
  reservedChs = 0;

  cout << "Starting process";
  if (steps) cout << " (" + toString(steps) + " steps)";
  cout << ": " << message << "." << endl;
  fflush(stdout);

  int nums = toString(steps).length();
  reservedChs = 14 + 2*nums;

  fmt = steps ?
    "%" + toString(nums) + "u of " + toString(steps) + " [%s] %4s" :
    string( "progress: %u" );
  
  
  progPrevTV = sec(nowTime());

}

/// \brief Updates the progress bar.
///
/// @param curstep Current step. Advances +1 if zero.
///
void
ProgressBar::update(long curstep){

  if ( !showme || !reservedChs ) return; // Uninitialized progress bar.
  
  step = curstep ? curstep+1 : step + 1;
  
  const double progNowTV = sec(nowTime());
  if ( progNowTV - progPrevTV > 0.1 )
    progPrevTV = progNowTV;
  else
    return;
  
  int progln = getwidth() - reservedChs;
  if ( progln <= 3 )  return; // if we have a very narrow terminal

  if ( steps && step >= steps ) {
    done();
    return;
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
  
  cout << string(waswidth+1, '\b') << outS;
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













int
allocateBigVolume( const Shape3D & shape, Volume8U & data ) {
  
  const std::string modname="Big volume allocation";
  
  int mapfile = 0;
  
  try { data.resize(shape); }
  catch(std::bad_alloc err) {

    warn (modname,
          "Failed to allocate memory. Trying to play with the memory mapped to file.");

    const char * tmpenv =  "TMPDIR";
    const std::string tempdesc = "You can control the path to the temporary file via"
                            " the environment variable \"" + (std::string) tmpenv + "\".";
    const char * tmpdir = getenv (tmpenv);
    if ( ! tmpdir ) tmpdir = "./" ;
    if ( strlen(tmpdir) >= 248 ) { // little under the maximum length (256)
      warn("temp file", (std::string)
        "The environment variable \"" + tmpenv + "\" is too long:\n"
        + tmpdir + "\n"
        "Switching to the default one \"./\". " + tempdesc);
      tmpdir = "./";
    }
    char tmpfilename[256]; // maximum file name length
    strcpy(tmpfilename, tmpdir);
    if ( tmpfilename[ strlen(tmpfilename)-1 ] != '/' )
      strcpy(tmpfilename+strlen(tmpfilename), "/");
    strcpy(tmpfilename+strlen(tmpfilename), "XXXXXX");

    mapfile = mkstemp(tmpfilename);
    if ( mapfile < 0 )
      throw_error("temp file", "Could not open temporary file"
                  " \"" + std::string(tmpfilename) + "\". " + tempdesc);
    if ( unlink(tmpfilename) < 0 )
    warn("temp file", "IMPORTANT! Could not unlink temporary file"
                      " \"" + std::string(tmpfilename) + "\". Don't forget"
                      " to delete it after the program has finished.");

    const off_t size = shape(0) * shape(1) * shape(2) * sizeof(uint8_t);

    if ( lseek(mapfile, size, SEEK_SET) != size ||
         write(mapfile, "END", 3) != 3 ) {
      // Should I pad the holes ???
      close (mapfile);
      throw_error("temp file",
                "Could not set size of the  temporary file"
                " \"" + std::string(tmpfilename) + "\". You need " + toString((long unsigned)size) +
                "bytes on the hard disk. " + tempdesc);
    }

    uint8_t * datap = (uint8_t*) mmap( 0, size, PROT_WRITE|PROT_READ, MAP_SHARED, mapfile, 0);
    if ( ! datap || datap == MAP_FAILED ) {
      close (mapfile);
      throw_error(modname, "Failed to mmap memory.");
    }

    data.reference( Volume8U(datap, shape, blitz::neverDeleteData) );
    
  }
  
  return mapfile;
  
}





Shape2D
ImageSizes(const Path & filename){
  Magick::Image imag;
  try { imag.ping(filename); }
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
    throw FloodErr(FloodErr::WARN, modname,
                  "Could not read tif from file\"" + filename + "\".");
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
  const Magick::PixelPacket * pixels = imag.getConstPixels(0,0,width,hight);

  storage.resize( hight, width );
  unsigned char * data = storage.data();

  for ( int k = 0 ; k < hight*width ; k++ )
    *data++ = (pixels++)->red;

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

    int wret = TIFFWriteRawStrip(image, 0, (void*) storage.data(), width*hight);
    TIFFClose(image);
    if ( -1 == wret )
      throw_error("save tiff image", "Could not save image to file \"" + filename + "\".");
    
  } else { 
    
    Magick::Image imag( Magick::Geometry(width, hight), "black" );
    imag.classType(Magick::DirectClass);
    imag.type( Magick::GrayscaleType );
    imag.depth(8);
    imag.magick("TIFF"); // saves to tif if not overwritten by the extension.

    const unsigned char *data = storage.data();
    Magick::PixelPacket * pixels = imag.getPixels(0,0,width,hight);
    for ( int k = 0 ; k < hight*width ; k++ )
      *pixels++ = Magick::PixelPacket( Magick::ColorGray( *data++ / 256 ) );

    imag.syncPixels();
    try { imag.write(filename); }
    catch ( Magick::Exception & error) {
      throw_error("save image IM", "Could not write image file\""+filename+"\"."
        " Caught Magick++ exception: \""+error.what()+"\".");
    }
      
  }

}


















