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


#ifndef _H_TOOLS_H_
#define _H_TOOLS_H_



#include <poptmx.h>
#include <pthread.h>
#include <queue>
#include <stdint.h>
#include <vector>

#include "blitz-long/blitz/array.h"






const uint8_t FILLED    = 0b10000000;
const uint8_t CHECKED   = 0b01000000;
const uint8_t SCHEDULED = 0b00100000;
const uint8_t ISGOOD    = 0b00010000;
const uint8_t ISBAD     = 0b00001000;






void prdn( int a );
void prdn( const std::string & msg );
void prdnc( int a );
void prdnc( const std::string & msg );



/// Error type.
class FloodErr{
public:
  /// Error severity
  typedef enum {
    WARN,                        ///< Warning
    ERR,                         ///< Error
    FATAL,                        ///< Fatal error
    NONE=0
  } ErrTp ;

private:
  std::string message;               ///< The message which describes the error
  std::string module;                ///< Name of the module where the error happened
  ErrTp terr;                   ///< Error severity

public:
  inline FloodErr(const std::string & mod = "", const std::string & msg = "", ErrTp _terr=NONE ) :
      terr(_terr), module(mod), message(msg) {}
  ErrTp type() const;     ///< Returns error type.
  void report() const;        ///< Reports the error to the ctderr.
};


/// \brief Constructs the error, reports and then throws it.
///
/// @param mod Module where the error happened.
/// @param msg Description.
///
void
throw_error(const std::string & mod, const std::string & msg);

/// \brief Reports and throws the bug.
///
/// Constructs the error with the message that it must never happen,
/// reports and then throws it.
///
/// @param mod Module where the error happened.
///
void
throw_bug(const std::string & mod);


/// \brief Constructs the warning and reports it.
///
/// @param mod Module where the error happened
/// @param msg Description
///
/// @return the error generated in warning.
///
FloodErr
warn(const std::string & mod, const std::string & msg);


/// \brief Constructs the error, reports and then exits.
///
/// @param mod Module where the error happened
/// @param msg Description
///
void
exit_on_error(const std::string & mod, const std::string & msg);



/// \brief %Path in the directory structure.
///
/// This class just a wrapper around std::string. It has some additional
/// useful functions which have sense only if the string represents a path.
/// Also a lot of platform specific code is hidden here.
class Path : public std::string {

private:

  const static std::string DIRSEPARATOR; ///< Directory separator
  void throw_unextracted(const std::string & element) const;
  ///< throws error in extracting an element.

public:

  /// \brief Constructor from string
  /// @param str the initializing string.
  inline Path(const std::string & str = std::string()) : std::string(str) {};

  /// \brief Constructor from C-string
  /// @param str the initializing string.
  inline Path(const char *str) : std::string(str) {};

  const static Path emptypath;

  std::string drive () const; ///< Extracts drive letter (has sense on MS WIN only).
  std::string dir () const;   ///< Extracts directory.
  std::string title () const; ///< Extracts title (filename without the extension).
  std::string dtitle () const;  ///< Extracts title with the preceding full path.
  std::string extension () const; ///< Extracts the extension.
  std::string name () const;  ///< Extracts file name.

  bool isdir() const;     ///< Tells if the path definitely represents the a directory.
  bool isabsolute() const;    ///< Tells if the path is absolute.

  Path & bedir();       ///< Makes the path to be the directory (adds DIRSEPARATOR).

  static const Path home();

};

/// \brief Prints type name.
/// To be used in the CLI parsing via "poptmx" library
/// @return type name.
std::string
type_desc (Path*);

/// \brief Converts the string "in" into the Path _val.
/// To be used in the CLI parsing via "poptmx" library
/// @param _val value to be updated.
/// @param in string to be parsed.
///
/// @return \c 1 if success, \c -1 otherwise.
int
_conversion (Path* _val, const std::string & in);


/// \brief Upgrades the filename with an addition.
///
/// Prefixes the base of the filename with the addition while the
/// path to the file remains unchanged.
///
/// @param path Path to upgrade.
/// @param addthis the upgrade string.
///
/// @return itself
///
Path
upgrade(const Path & path, const std::string & addthis);


/// \brief Adds the path to the filename.
///
/// If the filename is absolute then just returns it as is. If the
/// filename is relative then prefixes it with the "dir" and returns.
///
/// @param dir Path to be added.
/// @param file File name to be updated with the path.
///
/// @return Updated filename.
///
Path
cdpath(const Path & dir, const Path & file);




/// \brief 2D Array with 8bit data.
///
/// Two dimensional array of the ::uchar elements.
/// Used for sinograms, input and output images etc.
typedef blitz::Array<uint8_t,2> Map8U;

/// \brief 3D Array with 8bit data.
///
/// Three dimensional array of the ::uchar elements.
typedef blitz::Array<uint8_t,3> Volume8U;


/// \brief Shape of an 2D array.
typedef blitz::TinyVector<long int,2> Shape2D;

/// \brief Compare shapes.
///
/// @param sh1 first shape.
/// @param sh2 second shape.
///
/// @return \c true if the shapes are equal, \c false otherwise.
///
inline bool
operator==( const Shape2D & sh1, const Shape2D & sh2){
  return ( sh1(0)==sh2(0)  &&  sh1(1)==sh2(1) );
}

/// \brief Compare shapes.
///
/// @param sh1 first shape.
/// @param sh2 second shape.
///
/// @return \c false if the shapes are equal, \c true otherwise.
///
inline bool
operator!=( const Shape2D & sh1, const Shape2D & sh2){
  return ( sh1(0)!=sh2(0)  ||  sh1(1)!=sh2(1) );
}


/// \brief Shape of an 3D array.
typedef blitz::TinyVector<long int,3> Shape3D;

/// \brief Compare shapes.
///
/// @param sh1 first shape.
/// @param sh2 second shape.
/// @param sh3 third shape.
///
/// @return \c true if the shapes are equal, \c false otherwise.
///
inline bool
operator==( const Shape3D & sh1, const Shape3D & sh2){
  return ( sh1(0)==sh2(0)  &&  sh1(1)==sh2(1)  &&  sh1(2)==sh2(2) );
}

/// \brief Compare shapes.
///
/// @param sh1 first shape.
/// @param sh2 second shape.
///
/// @return \c false if the shapes are equal, \c true otherwise.
///
inline bool
operator!=( const Shape3D & sh1, const Shape3D & sh2){
  return ( sh1(0)!=sh2(0)  ||  sh1(1)!=sh2(1)  ||  sh1(2)!=sh2(2)  );
}




/// \brief Number of threads for the process.
///
/// @param _threads Requested number of threads (0 for auto).
///
/// @return Number of threads for the architecture where the process is running
/// if automatic number of threads was requested and just _threads if set in stone.
///
long
nof_threads(long _threads=0);





/// \brief Convert string to upper case
///
/// @param str Input string.
///
/// @return New string which represents input string converted into the upper case.
///
std::string
upper(std::string str);

/// \brief Convert string to lower case
///
/// @param str Input string.
///
/// @return New string which represents input string converted into the lower case.
///
std::string
lower(std::string str);


/// \brief Prints formatted message to string.
///
/// Like 'sprintf', but prints into the STL std::string
///
/// @param fmt format string.
/// @param ... whatever goes into the format string.
///
/// @return New STL string with the printed expression.
///
// Don't use the reference type "const string &" here: will
// not work on Windows
std::string
toString(const std::string fmt, ...);


/// \brief Prints the value to string.
///
/// @param n number to be printed
///
/// @return string with printed number
///
inline std::string toString (long double   n)   { return toString("%g", n); }
/// \cond
inline std::string toString (double        n)   { return toString("%g", n); }
inline std::string toString (float         n)   { return toString("%g", n); }
inline std::string toString (long int      n)   { return toString("%li", n); }
inline std::string toString (int           n)   { return toString("%i", n); }
inline std::string toString (long unsigned n)   { return toString("%lu", n); }
inline std::string toString (unsigned      n)   { return toString("%u", n); }
inline std::string toString (const Shape2D & shp) { return toString("%u, %u", shp(0), shp(1));}
inline std::string toString (const Shape3D & shp) { return toString("%u, %u, %u", shp(0), shp(1), shp(2));}
/// \endcond



bool string_to_argv(char const * str, int * argc_p, char *** argv_p);





/// \brief CLI progress bar (PB).
///
/// This class is used to visualize a progress.
class  ProgressBar {

private:
  bool showme;      ///< Tells if the bar should be shown.
  std::string message;    ///< Description of the process.
  uint64_t steps;     ///< Total number of steps in the progress.
  uint64_t step;        ///< Current step.
  int waswidth;     ///< The width of the terminal in the previous step
  int reservedChs;    ///< Symbols for the constant part of the PB.
  std::string fmt;      ///< Format string used to print the numbers.
  int getwidth();   ///< Returns current terminal width.
  double progPrevTV;   ///< time of the last update

public:

  ProgressBar(bool _showme=false, const std::string & _message="", uint64_t _steps=0);

  void update(uint64_t curstep=0);  ///< Updates the progress bar
  inline uint64_t progress() const {return step;}
  void done();          ///< Finalizes the progress bar.
  std::string print_line();

};


/// \brief Image sizes.
///
/// Reads image size and returns them as the Shape vector.
///
/// @param filename Image filename.
///
/// @return Image sizes.
///
Shape2D
ImageSizes(const Path & filename);


/// \brief Load image from preopened object.
///
/// Similar to ReadImage(), but reads into 8-bit array.
///
/// @param filename Image filename.
/// @param storage  Array to load the image into.
///
void
ReadImage(const Path & filename, Map8U & storage );

/// \brief Load image checking the shape.
///
/// Same as ReadImage(const Path & , Map &), but checks the image shape
/// before loading and throws the error if the shape is different from :::shp.
///
/// @param filename Image filename.
/// @param storage  Array to load the image into.
/// @param shp The expected shape.
///
void
ReadImage(const Path & filename, Map8U & storage, const Shape2D & shp);

/// \brief Save the array into integer image.
///
/// Stores the array in the integer-based image. If minval is equal to maxval
/// then the minimum and maximum values of the array data corresponds to black
/// and white respectively.
///
/// @param filename the name of the image file image.
/// @param storage the array to be written to the image.
/// @param minval the value corresponding to black.
/// @param maxval the value corresponding to white.
///
void
SaveImage(const Path & filename, const Map8U & storage);




class Point3D {

private:

  blitz::TinyVector<int,3> point;

public:

  inline Point3D(long int z=0, long int y=0, long int x=0) : point(z,y,x) {};

  inline int x() const { return point(2); }

  inline int y() const { return point(1); }

  inline int z() const { return point(0); }

  inline bool inVolume(const Shape3D & shape) const {
     return x() >= 0 && x() < shape(2) &&
            y() >= 0 && y() < shape(1) &&
            z() >= 0 && z() < shape(0);
  }

  inline int r2() const { return x()*x() + y()*y() + z()*z() ; }

  inline operator blitz::TinyVector<long int,3> () const { return point; }

  inline const Point3D operator+(const Point3D &other) const {
    return Point3D( point(0)+other.point(0), point(1)+other.point(1), point(2)+other.point(2) );
  }

  inline const Point3D operator-(const Point3D &other) const {
    return Point3D( point(0)-other.point(0), point(1)-other.point(1), point(2)-other.point(2) );
  }

  inline const bool operator==(const Point3D &other) const {
    return point(0) == other.point(0) &&
           point(1) == other.point(1) &&
           point(2) == other.point(2) ;
  }

};




inline std::string type_desc (Point3D*) {
  return "UINT:UINT:UINT";
}

inline int _conversion (Point3D* _val, const std::string & in) {
  int x, y, z;
  int scanres = sscanf( in.c_str(), "%i:%i:%i", &x, &y, &z);
  if (scanres != 3) // try , instead of :
    scanres = sscanf( in.c_str(), "%i,%i,%i", &x, &y, &z);
  if ( 3 != scanres || x<0 || y<0 || z<0 )
    return -1;
  *_val = Point3D(z,y,x); // it's not a mistake!
  return 1;
}


struct Sphere {
  Point3D center;
  int radius;
};


inline std::string type_desc (Sphere*) {
  return "UINT:UINT:UINT[:UINT]";
}

inline int _conversion (Sphere* _val, const std::string & in) {
  int x, y, z, r=0;
  int scanres = sscanf( in.c_str(), "%i:%i:%i:%i", &x, &y, &z, &r);
  if (scanres < 3) // try , instead of :
    scanres = sscanf( in.c_str(), "%i,%i,%i,%i", &x, &y, &z, &r);
  if ( scanres < 3 || r<0 )
    return -1;
  _val->center = Point3D(z,y,x);
  _val->radius = r;
  return 1;
}














class Point2D {

private:

  blitz::TinyVector<int,2> point;

public:

  inline Point2D(long int x=0, long int y=0) : point(x,y) {};

  inline int x() const { return point(0); }

  inline int y() const { return point(1); }

  inline bool inMap(const Shape2D & shape) const {
     return x() >= 0 && x() < shape(0) &&
            y() >= 0 && y() < shape(1) ;
  }

  inline int r2() const { return x()*x() + y()*y() ; }

  inline operator blitz::TinyVector<long int,2> () const { return point; }

  inline const Point2D operator+(const Point2D &other) const {
    return Point2D( point(0)+other.point(0), point(1)+other.point(1) );
  }

  inline const Point2D operator-(const Point2D &other) const {
    return Point2D( point(0)-other.point(0), point(1)-other.point(1));
  }

  inline const bool operator==(const Point2D &other) const {
    return point(0) == other.point(0) &&
           point(1) == other.point(1) ;
  }

};




inline std::string type_desc (Point2D*) {
  return "UINT:UINT";
}

inline int _conversion (Point2D* _val, const std::string & in) {
  int x, y;
  int scanres = sscanf( in.c_str(), "%i:%i", &x, &y);
  if (scanres != 2) // try , instead of :
    scanres = sscanf( in.c_str(), "%i,%i", &x, &y);
  if ( 2 != scanres || x<0 || y<0 )
    return -1;
  *_val = Point2D(y,x); // it's not a mistake!
  return 1;
}


struct Disk {
  Point2D center;
  int radius;
};


inline std::string type_desc (Disk*) {
  return "UINT:UINT[:UINT]";
}

inline int _conversion (Disk* _val, const std::string & in) {
  int x, y, r=0;
  int scanres = sscanf( in.c_str(), "%i:%i:%i", &x, &y, &r);
  if (scanres < 2) // try , instead of :
    scanres = sscanf( in.c_str(), "%i,%i,%i", &x, &y, &r);
  if ( scanres < 2 || r<0 )
    return -1;
  _val->center = Point2D(y,x);
  _val->radius = r;
  return 1;
}

















#endif // _H_TOOLS_H_

