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
/// @brief %3D flood fill algorithm.
///

#include "tools.h"
//#include "tools3d.h"
#include "algorithm3d.h"

#include <pthread.h>
#include <queue>
#include <algorithm>
#include <unistd.h>

#include<readline/readline.h>
#include<readline/history.h>

using namespace std;
using namespace blitz;


static const int run_threads( nof_threads() );
static Volume8U ivol;
static Volume8U wvol;
static Shape3D  volsh;




















struct clargs {

  Path command;               ///< Command name as it was invoked.
  bool interactive;
  vector<Path> inlist;        ///< List of all input images.
  Path out_filled;       ///< The mask for the output file names.
  Path out_inverted;
  Path out_mask;
  bool beverbose;       ///< Be verbose flag
  vector<Point3D> start;          ///< Point to start the fill at.
  vector<Sphere> stop;           ///< Point to stop at.

  // Parameters to use in the check function.
  uint radius;
  uint radiusM;
  uint color;
  uint minval;
  uint maxval;

  // sub-tables
  poptmx::OptionTable proc_table;
  poptmx::OptionTable save_table;
  poptmx::OptionTable color_table;

  poptmx::OptionTable table;

  clargs(int argc, char *argv[])
    : interactive(false)
    , beverbose(false)
    , radius(0)
    , radiusM(0)
    , color(0)
    , minval(-1)
    , maxval(-1)
    , proc_table()
    , save_table()
    , table("3D advanced flood fill.",
           "Advanced flood fill algorithm implemented for the 3D volume.")
  {

    proc_table
      .add(poptmx::NOTE, "Processing options:")
      .add(poptmx::OPTION, &start, 's', "start",
           "Starting point(s) of the fill procedure.", "")
      .add(poptmx::OPTION, &stop, 'S', "stop",
           "Stop sphere(s).", "Three coordinates and radius.")
      .add(poptmx::OPTION, &radius, 'r', "test-radius",
           "Radius of the test sphere.", "")
      .add(poptmx::OPTION, &radiusM, 'R', "mark-radius",
           "Radius of the fill sphere.", "")
      .add(poptmx::OPTION, &minval, 'm', "minval",
           "Minimum value.", "")
      .add(poptmx::OPTION, &maxval, 'M', "maxval",
           "Maximum value.", "") ;

    color_table
      .add(poptmx::OPTION, &color, 'c', "color",
           "Color of the fill volume.",
           "Instead of masking the volume, paints it with this color.")  ;

    save_table
      .add(poptmx::NOTE, "Saving options:")
      .add(color_table)
      .add(poptmx::OPTION, &out_filled, 'f', "outfilled",
           "Prefix to output filled.", "", out_filled)
      .add(poptmx::OPTION, &out_inverted, 'e', "outinvert",
           "Prefix to output not filled.", "", out_inverted)
      .add(poptmx::OPTION, &out_mask, 'b', "outmask",
           "Prefix to output mask.", "Outputs bit mask.", out_mask);


    table
      .add(poptmx::NOTE, "ARGUMENTS:")
      .add(poptmx::ARGUMENT, &inlist, "list", "List of the input images.", "")

      .add(poptmx::NOTE, "OPTIONS:")
      .add(poptmx::OPTION, &interactive, 'i', "interactive",
           "Accepts commands from the standard input.", "")
      .add(proc_table)
      .add(save_table)
      .add_standard_options(&beverbose);


    if ( ! table.parse(argc,argv) )
      exit(0);
    if ( ! table.count() ) {
      table.usage();
      exit(0);
    }

    command = table.name();

    // <list> : required argument.
    if ( ! table.count(&inlist) )
      exit_on_error(command, "Missing required argument: "+table.desc(&inlist)+".");

    if (interactive)
      return;

    check_proc(table);
    check_save(table);

  }


  void check_proc( poptmx::OptionTable & tab ) {

    if ( radius < 0 )
      throw_error(command, "Impossible parameter value : " + tab.desc(&radius) + "."
      " Cannot be negative.");
    if ( ! tab.count(&radiusM) )
      radiusM = radius;
    if ( tab.count(&minval) && ( minval < 0 || minval > 256 ) )
      throw_error(command, "Impossible parameter value : " + tab.desc(&minval) + "."
      " Must be in the range [0,256].");
    if ( tab.count(&maxval) && ( maxval < 0 || maxval > 256 ) )
      throw_error(command, "Impossible parameter value : " + tab.desc(&maxval) + "."
      " Must be in the range [0,256].");

    // point : required argument.
    if ( ! tab.count(&start) )
      throw_error(command, "Missing required argument: "+tab.desc(&start)+".");

  }


  void check_save( poptmx::OptionTable & tab ) {
    if ( ! ( tab.count(& out_filled) +
      tab.count(& out_inverted) +
      tab.count(& out_mask) ) )
      throw_error(command, "At least one of the two following arguments is required: "
      +tab.desc(&out_mask)+ ", "
      +tab.desc(&out_filled)+ ", "
      +tab.desc(&out_inverted)+ ".");
  }

};

























class ThreadDistributor {

private :

  static pthread_mutex_t lock;
  long int currentidx;
  vector<pthread_t> threads;

  void * arg;
  bool (*sub_routine0)  ();
  bool (*sub_routine1) (long int);
  bool (*sub_routine2) (void *, long int);
  bool sub_routine (long int idx) {
    if (sub_routine0) return sub_routine0();
    if (sub_routine1) return sub_routine1(idx);
    if (sub_routine2) return sub_routine2(arg, idx);
    return false;
  }

  ThreadDistributor( bool (*_sub_routine0)(),
                     bool (*_sub_routine1)(long int),
                     bool (*_sub_routine2) (void *, long int),
                     void * _arg)
    : currentidx(0)
    , arg(_arg)
    , sub_routine0(_sub_routine0)
    , sub_routine1(_sub_routine1)
    , sub_routine2(_sub_routine2)
  {}

  ThreadDistributor( bool (*_sub_routine)() )
    : ThreadDistributor(_sub_routine, 0, 0, 0)
  {}


  ThreadDistributor( bool (*_sub_routine)(long int) )
    : ThreadDistributor(0, _sub_routine, 0, 0)
  {}

  ThreadDistributor( bool (*_sub_routine) (void *, long int), void * _arg )
    : ThreadDistributor(0, 0, _sub_routine, _arg)
  {}


  long int distribute() {
    long int idx;
    pthread_mutex_lock( & lock );
    idx = currentidx++;
    pthread_mutex_unlock( & lock );
    return idx;
  }


  static void * in_thread (void * vdist) {
    ThreadDistributor * dist = (ThreadDistributor*) vdist;
    while ( dist->sub_routine(dist->distribute()) ) {}
    return 0;
  }


  void start() {
    pthread_t thread;
    for (int ith = 0 ; ith < run_threads ; ith++)
      if ( pthread_create( & thread, NULL, in_thread, this ) )
        warn("Thread operation", "Can't create thread.");
      else
        threads.push_back(thread);
  }

  void finish() {
    for (int ith = 0 ; ith < threads.size() ; ith++)
      pthread_join( threads[ith], 0);
  }

public:

  static void execute( bool (*_thread_routine) (void *, long int), void * _arg ) {
    ThreadDistributor dist(_thread_routine, _arg);
    dist.start();
    dist.finish();
  }

  static void execute( bool (*_thread_routine) (long int)) {
    ThreadDistributor dist(_thread_routine);
    dist.start();
    dist.finish();
  }

  static void execute( bool (*_thread_routine)()) {
    ThreadDistributor dist(_thread_routine);
    dist.start();
    dist.finish();
  }


};

pthread_mutex_t ThreadDistributor::lock = PTHREAD_MUTEX_INITIALIZER;




struct InitArgs {
  Volume8U & vol;
  uint crop;
  static const uint8_t bounds = CHECKED | ISGOOD | SCHEDULED;
  const long bg0;
  const long bg1;
  const long bg2;
  // const Range bgn;
  // const Range bg1;
  // const Range bg2;
  InitArgs(Volume8U & _vol, uint _crop)
    : vol(_vol)
    , crop(_crop)
    , bg0(vol.extent(0) - crop)
    , bg1(vol.extent(1) - crop)
    , bg2(vol.extent(2) - crop)
  {}
};

bool inThread_initXvol (void * _thread_args, long int idx) {

  // variant with range objects causes segfaults in here
  // have to use macro instead
  #define forSubVol( m1, M1, m2, M2 ) \
    for (int idx1 = m1 ; idx1 < M1 ; idx1++ ) \
      for (int idx2 = m2 ; idx2 < M2 ; idx2++ ) \
        args->vol(Point3D(idx, idx1, idx2)) = InitArgs::bounds;

  const InitArgs* args = (InitArgs*) _thread_args;
  if (idx >= args->vol.extent(0) )
    return false;
  if ( idx < args->crop  ||  idx >= args->vol.extent(0) - args->crop ) {
    //args->vol(idx, Range::all()   , Range::all()   ) = InitArgs::bounds;
    forSubVol(0, args->vol.extent(1), 0, args->vol.extent(2));
  } else {
    //args->vol(idx, args->bgn, args->bgn) = InitArgs::bounds;
    forSubVol(0, args->crop, 0, args->crop);
    //args->vol(idx, args->bg1, args->bgn) = InitArgs::bounds;
    forSubVol(args->bg1, args->vol.extent(1), 0, args->crop);
    //args->vol(idx, args->bgn, args->bg2) = InitArgs::bounds;
    forSubVol(0, args->crop, args->bg2, args->vol.extent(2));
    //args->vol(idx, args->bg1, args->bg2) = InitArgs::bounds;
    forSubVol(args->bg1, args->vol.extent(1), args->bg2, args->vol.extent(2));
  }

  #undef forSubVol
  return true;
}


bool inThread_zerowvol (long int idx) {
  if ( idx >= volsh(0) )
    return false;
  wvol(idx, Range::all(), Range::all()) = 0;
  return true;
}


bool inThread_wipewvol (long int idx) {
  if ( idx >= volsh(0) )
    return false;
  wvol(idx, Range::all(), Range::all()) &= ( ~ISBAD ) & ( ~ISGOOD );
  return true;
}



struct ApplyArgs {
  bool invert;
  uint color;
};

bool inThread_apply(void * _thread_args, long int idx) {
  if ( idx >= volsh(0))
    return false;
  ApplyArgs *  args = (ApplyArgs*) _thread_args;
  if (!args)
    throw_error("apply mask thread", "Inappropriate thread function arguments.");
  for ( long int sh1 = 0 ; sh1 < volsh(1) ; sh1++)
    for ( long int sh2 = 0 ; sh2 < volsh(2) ; sh2++)
      if ( ( ! args->invert  &&  ( ! ( wvol(idx, sh1, sh2) & FILLED ) ) )  ||
           (   args->invert  &&  (     wvol(idx, sh1, sh2) & FILLED   ) )  )
        ivol(idx, sh1, sh2) = args->color;
  return true;
}



struct ReadVolArgs {

  const vector<Path> & names;
  const Shape2D imgsh;
  ProgressBar bar;
  pthread_mutex_t proglock;

  ReadVolArgs(const vector<Path> & _names, bool verbose=false)
    : names(_names)
    , imgsh(volsh(1),volsh(2))
    , bar(verbose , "reading volume", volsh(0))
    , proglock(PTHREAD_MUTEX_INITIALIZER)
  {
    if (names.size() != volsh(0))
      throw_error("Reading volume", "Number of files not matching array size.");
  }

};

bool inThread_readVol (void * _thread_args, long int idx) {
  ReadVolArgs *  dist = (ReadVolArgs*) _thread_args;
  if (!dist)
    throw_error("read thread", "Inappropriate thread function arguments.");
  if ( idx >= dist->names.size())
    return false;
  Map8U img(dist->imgsh);
  ReadImage(dist->names[idx], img, dist->imgsh);
  ivol ( idx, Range::all(), Range::all() ) = img;
  pthread_mutex_lock(&dist->proglock);
  dist->bar.update();
  pthread_mutex_unlock(&dist->proglock);
  return true;
}



struct SaveVolArgs {

  const vector<Path> & names;
  const Shape2D imgsh;
  const Path out_filled;
  const Path out_inverted;
  const Path out_mask;
  const uint8_t color;
  ProgressBar bar;
  pthread_mutex_t proglock;

  SaveVolArgs(const vector<Path> & _names, bool verbose=false, uint8_t _color=0,
              const Path & _out_filled=Path(), const Path & _out_inverted=Path(), const Path & _out_mask=Path() )
    : names(_names)
    , imgsh(volsh(1),volsh(2))
    , out_filled(_out_filled)
    , out_inverted(_out_inverted)
    , out_mask(_out_mask)
    , color(_color)
    , bar(verbose , "saving volume", volsh(0))
    , proglock(PTHREAD_MUTEX_INITIALIZER)
  {
    if (names.size() != volsh(0))
      throw_error("Saving volume", "Number of files not matching array size.");
    if ( out_filled.empty() && out_inverted.empty() && out_mask.empty() )
      throw_error("Saving volume", "Nothing to output.");
  }

};


bool inThread_saveVol (void * _thread_args, long int idx) {

  SaveVolArgs *  dist = (SaveVolArgs*) _thread_args;
  if (!dist)
    throw_error("read thread", "Inappropriate thread function arguments.");
  if ( idx >= dist->names.size() )
    return false;

  const Map8U wmg = wvol (idx, Range::all(), Range::all());
  Map8U img(dist->imgsh);

  if ( ! dist->out_filled.empty() ) {
    for ( long int sh0 = 0 ; sh0 < dist->imgsh(0) ; sh0++)
      for ( long int sh1 = 0 ; sh1 < dist->imgsh(1) ; sh1++)
        img(sh0, sh1) =  ( wmg(sh0, sh1) & FILLED ) ? ivol( idx, sh0, sh1 ) : dist->color ;
    SaveImage(dist->out_filled + dist->names[idx].name(), img);
  }
  if ( ! dist->out_inverted.empty() ) {
    for ( long int sh0 = 0 ; sh0 < dist->imgsh(0) ; sh0++)
      for ( long int sh1 = 0 ; sh1 < dist->imgsh(1) ; sh1++)
        img(sh0, sh1) =  ( wmg(sh0, sh1) & FILLED ) ? dist->color : ivol( idx, sh0, sh1 ) ;
    SaveImage(dist->out_inverted + dist->names[idx].name(), img);
  }
  if ( ! dist->out_mask.empty() ) {
    SaveImage(dist->out_mask + dist->names[idx].name(), wmg);
  }

  pthread_mutex_lock(&dist->proglock);
  dist->bar.update();
  pthread_mutex_unlock(&dist->proglock);

  return true;

}



























void saveSlice(const Point3D & cross, bool original=false, const string & prefix = string()) {

  if ( ! cross.inVolume(volsh) )
    throw_error("slice", "Point outside volume.");
  const Volume8U & voltosave = original ? ivol : wvol;

  Map8U yz = voltosave( Range::all(), Range::all(), cross.x() ).copy();
  SaveImage(prefix + ".yz.tif", yz);
  Map8U xz = voltosave( Range::all(), cross.y(), Range::all() ).copy();
  SaveImage(prefix + ".xz.tif", xz);
  Map8U xy = voltosave( cross.z(), Range::all(), Range::all() ).copy();
  SaveImage(prefix + ".xy.tif", xy);

}


#include <signal.h>

void sig2handler(int signum) {

  Point3D pnt;
  char line[256];

  FILE *funcf = fopen( ".sliceme.txt", "r" );

  if ( ! funcf ||
       1 != fscanf ( funcf,  "%254[^\n]\n", line) ||
       1 != _conversion (&pnt, line) ||
       ! pnt.inVolume(volsh))
    pnt = Point3D(volsh(2)/2, volsh(1)/2, volsh(0)/2);

  if (funcf)
    fclose(funcf);

  saveSlice(pnt);

  printf("Written slices through (%i, %i, %i) point."
         " Files are named \".xy.tif\", \".xz.tif\" and \".yz.tif\".\n",
         pnt.x(), pnt.y(), pnt.z() );

}
























int main(int argc, char *argv[]) {

  clargs args(argc, argv);

  struct sigaction act2;
  act2.sa_handler = sig2handler;
  sigaction(SIGUSR2, &act2, NULL);

  if (args.beverbose)
    printf("My PID is %i. You can send me USR2 signal to save current"
           " mask slices through the first start point.\n", getpid() );

  const Shape2D imgsh = ImageSizes(args.inlist[0]);

  const uint mrad = max(args.radius, args.radiusM);
  volsh=Shape3D(args.inlist.size(), imgsh(0), imgsh(1));
  Volume8U xivol( Shape3D(volsh(0)+2*mrad, volsh(1)+2*mrad, volsh(2)+2*mrad) );
  Volume8U xwvol( Shape3D(volsh(0)+2*mrad, volsh(1)+2*mrad, volsh(2)+2*mrad) );
  ivol.reference(xivol(Range(mrad, mrad+volsh(0)-1), Range(mrad, mrad+volsh(1)-1), Range(mrad, mrad+volsh(2)-1)));
  wvol.reference(xwvol(Range(mrad, mrad+volsh(0)-1), Range(mrad, mrad+volsh(1)-1), Range(mrad, mrad+volsh(2)-1)));

  InitArgs initArgs(xwvol, mrad);
  ThreadDistributor::execute(inThread_initXvol, &initArgs);
  ThreadDistributor::execute(inThread_zerowvol);
  ReadVolArgs readArgs(args.inlist, args.beverbose);
  ThreadDistributor::execute(inThread_readVol, &readArgs);


  if ( args.interactive ) {

    char* flns;
    char** aargv;
    int aargc;
    string sflns;
    ProcDistributor * dist = 0;

    while ( (flns = readline("enter a command > ")) ) {

      if (flns && *flns) {
        add_history (flns);
        sflns = string(flns);
      } else {
        sflns = " ";
      }

      string lns;
      istringstream tokenStream(sflns);

      while (getline(tokenStream, lns, ';')) {

        try {


          if ( string_to_argv(lns.c_str(), &aargc, &aargv )  || ! aargc ) {

            if (dist)
              dist->update_process();
            else
              printf("  No process is running.\n");

            printf ( "  Possible commands are: start, stop, apply, save, slice.\n" );


          } else if ( string(aargv[0]) == "save" ) {

            if (dist)
              throw_error("save", "A process is running. Stop it first.");

            poptmx::OptionTable table = args.save_table;
            table.parse(aargc, aargv);
            args.check_save(table);

            SaveVolArgs saveArgs(args.inlist, args.beverbose, args.color,
                                 args.out_filled, args.out_inverted, args.out_mask);
            ThreadDistributor::execute(inThread_saveVol, &saveArgs);

          } else if ( string(aargv[0]) == "apply" ) {

            if (dist)
              throw_error("apply", "A process is running. Stop it first.");

            args.color=0;
            bool invert(false);
            poptmx::OptionTable table = args.color_table;
            table.add(poptmx::OPTION, &invert, 'i', "invert", "Inverts mask", "");
            table.parse(aargc, aargv);

            ApplyArgs apply = {invert, args.color};
            ThreadDistributor::execute(inThread_apply, &apply);


          } else if ( string(aargv[0]) == "reset" ) {

            if (dist)
              throw_error("save", "A process is running. Stop it first.");
            ThreadDistributor::execute(inThread_zerowvol);


          } else if ( string(aargv[0]) == "start" ) {

            if (dist)
              throw_error("save", "A process is already running. Stop it first.");

            args.radius = 0;
            args.radiusM = 0;
            args.color = 0;
            args.minval = -1;
            args.maxval = -1;
            args.start.clear();
            args.stop.clear();

            poptmx::OptionTable table = args.proc_table;
            table.parse(aargc, aargv);
            args.check_proc(table);

            ThreadDistributor::execute(inThread_wipewvol);
            dist = new ProcDistributor(ivol, wvol, args.start, args.stop,
                                       args.radius, args.radiusM, args.minval, args.maxval,
                                       run_threads, args.beverbose);
            dist->start_process();

          } else if ( string(aargv[0]) == "stop" ) {

            if (!dist)
              throw_error("stop", "No process running.");

            dist->stopProcessing=true;
            dist->finish_process();

          } else if ( string(aargv[0]) == "slice" ) {

            Point3D cross;
            bool original=false;
            string prefix;
            poptmx::OptionTable table;
            table
              .add(poptmx::ARGUMENT, &cross, "point", "slices through this point", "")
              .add(poptmx::OPTION, &original, 'o', "original",
                  "Save slices through the data volume, not mask.", "")
              .add(poptmx::OPTION, &original, 'p', "prefix",
                  "Prefix of output files.", "");
            table.parse(aargc, aargv);

            saveSlice(cross, original, prefix);

          } else {

            printf ( "  Unknown command \"%s\". Possible commands are: start, stop, apply, save, slice.\n", aargv[0] );

          }

        }

        catch (FloodErr err) {}
        catch (poptmx::Err err) {}

      }

      free(flns);
      free(aargv);

    }

  } else {

    ProcDistributor(ivol, wvol, args.start, args.stop,
                    args.radius, args.radiusM, args.minval, args.maxval,
                    run_threads, args.beverbose)
      .process();

    SaveVolArgs saveArgs(args.inlist, args.beverbose, args.color,
                         args.out_filled, args.out_inverted, args.out_mask);
    ThreadDistributor::execute(inThread_saveVol, &saveArgs);

  }

}
