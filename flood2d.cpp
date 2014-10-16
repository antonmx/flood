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

#include <pthread.h>
#include <queue>
#include <algorithm>
#include <unistd.h>

using namespace std;
using namespace blitz;

struct GlobalVariables {
  
  Path command;               ///< Command name as it was invoked.
  Path input;        ///< List of all input images.
  Path out_filled;       ///< The mask for the output file names.
  Path out_inverted;
  Path out_mask;
  bool beverbose;       ///< Be verbose flag
  std::vector<Point2D> start;          ///< Point to start the fill at.
  std::vector<Disk> stop;           ///< Point to stop at.
//  int nextdim;          ///< dimensions for the neibour.
  
  
  const int run_threads;
  
  Map8U ivol;
  Map8U wvol;
  
  // Parameters to use in the check function. 
  int radius;
  int radiusM;
  unsigned int color;
  int minval;
  int maxval;
  int minggrad; // square of the grad!
  int maxggrad; // square of the grad!
  
  blitz::Array < std::vector<Point2D>, 8 > newPoints;    
  blitz::Array < std::vector<Point2D>, 8> markPoints;
  
  GlobalVariables() :
    beverbose(false),
//  nextdim(1),
    run_threads(nof_threads()),
    radius(0),
    radiusM(0),
    color(0),
    minval(-1),
    maxval(-1),
    minggrad(-1), // square of the grad!
    maxggrad(-1), // square of the grad!    
    newPoints (2,2,2,2),
    markPoints(2,2,2,2)
  {}
  
} g;


/// \CLARGS
struct clargs {
  
  poptmx::OptionTable table;
  
  /// \CLARGSF
  clargs(int argc, char *argv[]);
  
};



clargs::clargs(int argc, char *argv[]) :
  table("2D advanced flood fill.",
        "Advanced flood fill algorithm implemented for the 2D volume.")
{

  table
    .add(poptmx::NOTE, "ARGUMENTS:")
    .add(poptmx::ARGUMENT, &g.input, "image", "Input image.", "")

    .add(poptmx::NOTE, "OPTIONS:")
    .add(poptmx::OPTION, &g.out_filled, 'f', "outfilled",
         "Prefix to output filled.", "", g.out_filled)
    .add(poptmx::OPTION, &g.out_inverted, 'e', "outinvert",
         "Prefix to output not filled.", "", g.out_inverted)
    .add(poptmx::OPTION, &g.out_mask, 'b', "outmask",
         "Prefix to output mask.", "Outputs bit mask.", g.out_mask)
    .add(poptmx::OPTION, &g.start, 's', "start",
         "Starting point(s) of the fill procedure.", "")
    .add(poptmx::OPTION, &g.stop, 'S', "stop",
         "Stop sphere(s).", "Three coordinates and radius.")
    .add(poptmx::OPTION, &g.radius, 'r', "test-radius",
         "Radius of the test sphere.", "")    
    .add(poptmx::OPTION, &g.radiusM, 'R', "mark-radius",
         "Radius of the fill sphere.", "")    
    .add(poptmx::OPTION, &g.color, 'c', "color",
         "Color of the fill volume.",
         "Instead of masking the volume, paints it with this color.")        
    .add(poptmx::OPTION, &g.minval, 'm', "minval",
         "Minimum value.", "")
    .add(poptmx::OPTION, &g.maxval, 'M', "maxval",
         "Maximum value.", "")
    .add(poptmx::OPTION, &g.minggrad, 'g', "mingrad",
         "Minimum absolute value of the gradient.", "")
    .add(poptmx::OPTION, &g.maxggrad, 'G', "maxgrad",
         "Maximum absolute value of the gradient.", "")    
//    .add(poptmx::OPTION, &nextdim, 'n', "next",
//         "Dimension for finding neigbours.",
//         "Can be 1 - lines (6 neigbours), 2 - planes (14) and 3 - cube (26)")
    .add_standard_options(&g.beverbose);
    
    
    if ( ! table.parse(argc,argv) )
      exit(0);
    if ( ! table.count() ) {
      table.usage();
      exit(0);
    }
  
  g.command = table.name();
      
  // <list> : required argument.
  if ( ! table.count(&g.input) )
    exit_on_error(g.command, "Missing required argument: "+table.desc(&g.input)+".");
  
  // point : required argument.
  if ( ! table.count(&g.start) )
    exit_on_error(g.command, "Missing required argument: "+table.desc(&g.start)+".");
  
//  if ( nextdim < 1 || nextdim > 3 )
//    exit_on_error(command, "Impossible parameter value : " + table.desc(&nextdim) + ".");    
  
  if ( ! ( table.count(& g.out_filled) +
           table.count(& g.out_inverted) +
           table.count(& g.out_mask) ) )
    exit_on_error(g.command, "At least one of the two following arguments is required: "
                            +table.desc(&g.out_mask)+ ", "
                            +table.desc(&g.out_filled)+ ", "
                            +table.desc(&g.out_inverted)+ ".");

  if ( g.radius < 0 )
    exit_on_error(g.command, "Impossible parameter value : " + table.desc(&g.radius) + "."
                           " Cannot be negative.");  
  if ( ! table.count(&g.radiusM) )
    g.radiusM = g.radius;
  if ( table.count(&g.minval) && ( g.minval < 0 || g.minval > 256 ) )
    exit_on_error(g.command, "Impossible parameter value : " + table.desc(&g.minval) + "."
                           " Must be in the range [0,256].");  
  if ( table.count(&g.maxval) && ( g.maxval < 0 || g.maxval > 256 ) )
    exit_on_error(g.command, "Impossible parameter value : " + table.desc(&g.maxval) + "."
                           " Must be in the range [0,256].");  
  if ( table.count(&g.minggrad) && ( g.minggrad < 0 || g.minggrad > 256 ) )
    exit_on_error(g.command, "Impossible parameter value : " + table.desc(&g.minggrad) + "."
                           " Must be in the range [0,256].");
  if ( table.count(&g.minggrad) )
    g.minggrad *= g.minggrad;
  if ( table.count(&g.maxggrad) && ( g.maxggrad < 0 || g.maxggrad > 256 ) )
    exit_on_error(g.command, "Impossible parameter value : " + table.desc(&g.maxggrad) + "."
                           " Must be in the range [0,256].");  
  if ( table.count(&g.maxggrad) )
    g.maxggrad *= g.maxggrad;
  
}









int ggradient( const Point2D & pnt, const Map8U & map ) {

  int sum=0, num=0;
  Point2D tpnt;
  
  #define adspn(pnt, x, y) \
    tpnt = pnt + Point2D(x, y); \
    if( tpnt.inMap(map.shape()) ) { \
      int gg = map(pnt) - map (tpnt) ; \
      sum+=gg*gg; \
      num++; \
    } 
  
  adspn(pnt, 1, 0);
  adspn(pnt,-1, 0);
  adspn(pnt, 0, 1);
  adspn(pnt, 0,-1);
  
  return num ? sum/num : 0 ;
      
}


  
bool subcheck( const Point2D & pnt) {

  if ( g.wvol(pnt) & ISBAD )
    return false;
  if ( g.wvol(pnt) & ISGOOD )
    return true;

  const int ggrad = ( g.minggrad >= 0 || g.maxggrad >= 0 ) ? ggradient(pnt, g.ivol) : -1 ;

  bool pass =
    ( g.minval < 0  ||  g.ivol(pnt) >= g.minval ) &&
    ( g.maxval < 0  ||  g.ivol(pnt) <= g.maxval ) &&
    ( g.minggrad < 0  ||  ggrad >= g.minggrad ) &&
    ( g.maxggrad < 0  ||  ggrad <= g.maxggrad ) ;

  g.wvol(pnt) |= pass ? ISGOOD : ISBAD ; // can I do it in multithreaded?

  return pass;

}



inline long int spn(const Point2D & pnt, int xx, int yy) {
  Point2D pntinhere( pnt+Point2D(xx,yy) );
  return
    pntinhere.inMap(g.wvol.shape())  &&  g.wvol(pntinhere) & CHECKED
      ? 1l : 0l ;
}

// THIS IS THE KEY FUNCTION
int checkMe(const Point2D & pnt) {
  
  static const int radius2 = g.radius * g.radius;
  
  const vector<Point2D> & tocheck = g.newPoints(
    spn(pnt, 1, 0), spn(pnt, -1, 0), spn(pnt, 0, 1), spn(pnt, 0, -1) );
  vector<Point2D>::const_iterator it = tocheck.begin();  
  while ( it != tocheck.end() )  {
    Point2D tpnt( (*it++) + pnt );
    if ( tpnt.inMap(g.ivol.shape()) && ! subcheck(tpnt) )
      return 0;    
  }

  return 1;
  
}


int nearestBad2(const Point2D & pnt) {
  
  static const int radius2 = g.radius * g.radius;
  int closest2 = radius2+1;

  const vector<Point2D> & tocheck = g.newPoints(
    spn(pnt, 1, 0), spn(pnt, -1, 0), spn(pnt, 0, 1), spn(pnt, 0, -1) );
  vector<Point2D>::const_iterator it = tocheck.begin();  
  while ( it != tocheck.end() )  {
    const Point2D tpnt( *it + pnt );
    const int tr2=it->r2();
    if ( tpnt.inMap(g.ivol.shape()) && tr2 < closest2 && ! subcheck(tpnt) )
      closest2 = tr2 ;
    it++;    
  }

  return closest2;

}











class ProcDistributor {
    
private:
  
  list<Point2D> schedule;
  list<Point2D> inwork;
  
  static pthread_mutex_t picklock;
  static pthread_cond_t check_again;
  static pthread_mutex_t proglock;
  
  ProgressBar bar;
  
public:
  
  ProcDistributor() {
      bar=ProgressBar( g.beverbose , "processing volume", g.ivol.size() );
      for ( int cpnt=0 ; cpnt < g.start.size() ; cpnt++ )
        schedule.push_back( g.start[cpnt] );  
  }
      
  bool distribute( Point2D * pnt ) {
    
    bool ret = true;
    
    pthread_mutex_lock( & picklock );
    
    while ( schedule.empty() && ! inwork.empty() )
      pthread_cond_wait(&check_again, &picklock);
    
    if ( schedule.empty() ) {
      ret = false;
    } else {
      *pnt = schedule.front();
      schedule.pop_front();
      inwork.push_back(*pnt);
    }
   
    pthread_mutex_unlock( & picklock );
    pthread_cond_signal( & check_again ) ; // do i need it here?
    
    return ret;
    
  }
  
  void collect( const Point2D & pnt, int ffrad ) {
    
    list<Point2D> add_to_schedule;
    
    if ( ffrad ) {
      
      Point2D pntt;
      
      #define scheduleMe( shift1, shift2 ) \
        pntt = pnt + Point2D(shift1, shift2); \
        if ( pntt.inMap(g.wvol.shape() ) && \
           ! ( g.wvol(pntt) & SCHEDULED ) ) { \
        add_to_schedule.push_back(pntt); \
        g.wvol(pntt) |= SCHEDULED; \
      }      
   
      scheduleMe( 1, 0);
      scheduleMe(-1, 0);
      scheduleMe( 0, 1);
      scheduleMe( 0,-1);

/*      
      if (g.nextdim > 1) {
      
        scheduleMe( 0, 1, 0);
        scheduleMe( 0,-1, 0);
        scheduleMe( 0, 1, 1);
        scheduleMe( 0,-1,-1);
        
        scheduleMe( 1, 0, 0);
        scheduleMe(-1, 0, 0);
        scheduleMe( 1, 0, 1);
        scheduleMe(-1, 0,-1);
        
        scheduleMe( 1, 0, 0);
        scheduleMe(-1, 0, 0);
        scheduleMe( 1, 1, 0);
        scheduleMe(-1,-1, 0);
        
      }
      
      if (g.nextdim > 2) {
              
        scheduleMe( 1, 1, 1);
        scheduleMe(-1, 1, 1);
        scheduleMe( 1,-1, 1);
        scheduleMe(-1,-1, 1);
        scheduleMe( 1, 1,-1);
        scheduleMe(-1, 1,-1);
        scheduleMe( 1,-1,-1);
        scheduleMe(-1,-1,-1);
      
      }
*/


      const vector<Point2D> & tomark = g.markPoints(
        spn(pnt, 1, 0), spn(pnt, -1, 0), spn(pnt, 0, 1), spn(pnt, 0, -1) );
      vector<Point2D>::const_iterator it = tomark.begin();  
      while ( it != tomark.end() ) {
        Point2D tpnt( (*it++) + pnt );
        if ( tpnt.inMap(g.ivol.shape()) )
          g.wvol(tpnt) |= FILLED;  // it works MUCH faster if done without thread locking. Has seen no difference
      }       
      
    } else {
 
      // mark all on the boundary to the nearest.
      
    }
    
    
    pthread_mutex_lock( & picklock );    
    schedule.splice(schedule.end(), add_to_schedule);
    inwork.remove(pnt);
    if ( ffrad )      
      g.wvol(pnt) |= CHECKED;
    pthread_mutex_unlock( & picklock );
    pthread_cond_signal( & check_again ) ;    
    
    
    pthread_mutex_lock( & proglock );
    bar.update();
    pthread_mutex_unlock( & proglock );

    
  }
  
  void finilizeProgressBar() {
    bar.update();
  }
  

  
};

pthread_mutex_t ProcDistributor::picklock = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t ProcDistributor::check_again = PTHREAD_COND_INITIALIZER;
pthread_mutex_t ProcDistributor::proglock = PTHREAD_MUTEX_INITIALIZER;



void * in_proc_thread (void * _thread_args) {
  
  ProcDistributor *  dist = (ProcDistributor*) _thread_args;
  if (!dist)
    throw_error("process thread", "Inappropriate thread function arguments.");
  
  Point2D pnt; 
  while ( dist->distribute(&pnt)  )    
    dist->collect( pnt, checkMe(pnt) );    
    
}






void sig2handler(int signum) {  
  SaveImage(".tmp.tif", g.wvol);
}



#include <signal.h>


int main(int argc, char *argv[]) {

  clargs args(argc, argv);
  
  struct sigaction act2;
  act2.sa_handler = sig2handler;
  sigaction(SIGUSR2, &act2, NULL);

  if (g.beverbose) 
    printf("My PID is %i. You can send me USR2 signal to save current"
    " mask slices through the first start point.\n", getpid() );
  
  ReadImage( g.input, g.ivol );
  Shape2D shape = g.ivol.shape();
  g.wvol.resize(shape);
  g.wvol = 0;
  
  for ( int cpnt=0 ; cpnt < g.start.size() ; cpnt++ )
    if ( ! g.start[cpnt].inMap(shape) )
      exit_on_error(g.command, "At least one starting point is outside the volume.");
  
  // prepare the shift-volumes  
  
  const int maxrad = max(g.radius, g.radiusM);
  const int radius2 = g.radius * g.radius;
  const int radiusM2 = g.radiusM * g.radiusM;
  
  if ( g.radius == g.radiusM )
    g.markPoints.reference( g.newPoints );
    
  for ( int x = -maxrad ; x <= maxrad ; x++ )
    for ( int y = -maxrad ; y <= maxrad ; y++ ) {
        
        const Point2D pnt(x,y);
        
        if ( pnt.r2() <= radius2 )

          for ( long _P0 = 0 ; _P0 <= 1 ; _P0++)
            if ( ! _P0  || ( pnt - Point2D( 1, 0) ).r2() > radius2 )
              for ( long _M0 = 0 ; _M0 <= 1 ; _M0++)
                if ( ! _M0  || ( pnt - Point2D(-1, 0) ).r2() > radius2 )
                  for ( long _0P = 0 ; _0P <= 1 ; _0P++)
                    if ( ! _0P  || ( pnt - Point2D( 0, 1) ).r2() > radius2 )                            
                      for ( long _0M = 0 ; _0M <= 1 ; _0M++)
                        if ( ! _0M  || ( pnt - Point2D( 0,-1) ).r2() > radius2 )
                          g.newPoints( _P0, _M0, _0P, _0M ).push_back(pnt);
                    
        if ( g.radius != g.radiusM  &&  pnt.r2() <= radiusM2 )

          for ( long _P0 = 0 ; _P0 <= 1 ; _P0++)
            if ( ! _P0  || ( pnt - Point2D( 1, 0) ).r2() > radiusM2 )
              for ( long _M0 = 0 ; _M0 <= 1 ; _M0++)
                if ( ! _M0  || ( pnt - Point2D(-1, 0) ).r2() > radiusM2 )
                  for ( long _0P = 0 ; _0P <= 1 ; _0P++)
                    if ( ! _0P  || ( pnt - Point2D( 0, 1) ).r2() > radiusM2 )                            
                      for ( long _0M = 0 ; _0M <= 1 ; _0M++)
                        if ( ! _0M  || ( pnt - Point2D( 0,-1) ).r2() > radiusM2 )
                          g.markPoints( _P0, _M0, _0P, _0M ).push_back(pnt);

    }     

    
  for ( int cpnt=0 ; cpnt < g.stop.size() ; cpnt++ ) {
    int rs2 = g.stop[cpnt].radius * g.stop[cpnt].radius;
    for ( int x = - g.stop[cpnt].radius ; x <= g.stop[cpnt].radius ; x++ )
      for ( int y = - g.stop[cpnt].radius ; y <= g.stop[cpnt].radius ; y++ ) { 
        const Point2D pnt = Point2D(x,y) + g.stop[cpnt].center;
        if ( Point2D(x,y).r2() <= rs2  &&  pnt.inMap(shape) )
          g.wvol(pnt) |= ISBAD | SCHEDULED;
      }
  }          
   
  // process
  vector<pthread_t> threads(g.run_threads);
  ProcDistributor procdist;
  for (int ith = 0 ; ith < g.run_threads ; ith++)
    if ( pthread_create( & threads[ith], NULL, in_proc_thread, &procdist ) )
      throw_error("Process data", "Can't create thread.");
  for (int ith = 0 ; ith < threads.size() ; ith++)
    pthread_join( threads[ith], 0);
  procdist.finilizeProgressBar();
  
  Map8U img(g.ivol.shape());
    
  if ( ! g.out_filled.empty() ) {
    for ( long int sh0 = 0 ; sh0 < shape(0) ; sh0++)
      for ( long int sh1 = 0 ; sh1 < shape(1) ; sh1++)
        img(sh0, sh1) =  ( g.wvol(sh0, sh1) & FILLED ) ? g.ivol( sh0, sh1 ) : g.color ;
    SaveImage(g.out_filled, img);
  }
      
  if ( ! g.out_inverted.empty() ) {
    for ( long int sh0 = 0 ; sh0 < shape(0) ; sh0++)
      for ( long int sh1 = 0 ; sh1 < shape(1) ; sh1++)
        img(sh0, sh1) =  ( g.wvol(sh0, sh1) & FILLED ) ? g.color : g.ivol( sh0, sh1 ) ;
      SaveImage(g.out_inverted, img);
  }
  
    if ( ! g.out_mask.empty() )
      SaveImage(g.out_mask, g.wvol);

    
}
