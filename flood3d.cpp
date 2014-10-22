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
#include "tools3d.h"
#include <blitz/rand-tt800.h>

#include <pthread.h>
#include <queue>
#include <algorithm>
#include <unistd.h>

#include<readline/readline.h>
#include<readline/history.h>

using namespace std;
using namespace blitz;




GlobalVariables g;


int ggradient( const Point3D & pnt, const Volume8U & vol ) {

  int sum=0, num=0;
  Point3D tpnt;
  
  #define adspn(pnt, x, y, z) \
    tpnt = pnt + Point3D(x, y, z); \
    if( tpnt.inVolume(vol.shape()) ) { \
      int gg = vol(pnt) - vol (tpnt) ; \
      sum+=gg*gg; \
      num++; \
    } 
  
  adspn(pnt, 1, 0, 0);
  adspn(pnt,-1, 0, 0);
  adspn(pnt, 0, 1, 0);
  adspn(pnt, 0,-1, 0);
  adspn(pnt, 0, 0, 1);
  adspn(pnt, 0, 0,-1);
  
  return num ? sum/num : 0 ;
      
}


  
bool subcheck( const Point3D & pnt) {

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



inline long int spn(const Point3D & pnt, int xx, int yy, int zz) {
  Point3D pntinhere( pnt+Point3D( xx,yy,zz) );
  return
    pntinhere.inVolume(g.wvol.shape())  &&  g.wvol(pntinhere) & CHECKED
      ? 1l : 0l ;
}

// THIS IS THE KEY FUNCTION
int checkMe(const Point3D & pnt) {
  
  static const int radius2 = g.radius * g.radius;
  
  const vector<Point3D> & tocheck = g.newPoints( spn(pnt, 1, 0, 0), spn(pnt,-1, 0, 0),
                                                 spn(pnt, 0, 1, 0), spn(pnt, 0,-1, 0),
                                                 spn(pnt, 0, 0, 1), spn(pnt, 0, 0,-1) );
  vector<Point3D>::const_iterator it = tocheck.begin();  
  while ( it != tocheck.end() )  {
    Point3D tpnt( (*it++) + pnt );
    if ( tpnt.inVolume(g.ivol.shape()) && ! subcheck(tpnt) )
      return 0;    
  }

  return 1;
  
}


int nearestBad2(const Point3D & pnt) {
  
  static const int radius2 = g.radius * g.radius;
  int closest2 = radius2+1;
  
  const vector<Point3D> & tocheck = g.newPoints( spn(pnt, 1, 0, 0), spn(pnt,-1, 0, 0),
                                                 spn(pnt, 0, 1, 0), spn(pnt, 0,-1, 0),
                                                 spn(pnt, 0, 0, 1), spn(pnt, 0, 0,-1) );
  vector<Point3D>::const_iterator it = tocheck.begin();  
  while ( it != tocheck.end() )  {
    const Point3D tpnt( *it + pnt );
    const int tr2=it->r2();
    if ( tpnt.inVolume(g.ivol.shape()) && tr2 < closest2 && ! subcheck(tpnt) )
      closest2 = tr2 ;
    it++;    
  }

  return closest2;

}











class ProcDistributor {
    
private:
  
  list<Point3D> schedule;
  list<Point3D> inwork;
  
  static pthread_mutex_t picklock;
  static pthread_cond_t check_again;
  static pthread_mutex_t proglock;
  
public:
  
  ProgressBar bar;
  
  ProcDistributor() {
      bar=ProgressBar( g.beverbose , "processing volume", g.ivol.size() );
      for ( int cpnt=0 ; cpnt < g.start.size() ; cpnt++ )
        schedule.push_back( g.start[cpnt] );  
  }
      
  bool distribute( Point3D * pnt ) {
    
    pthread_mutex_lock( & picklock );
    
    bool ret = true;
    
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
  
  void collect( const Point3D & pnt, int ffrad ) {
    
    list<Point3D> add_to_schedule;
    
    if ( ffrad ) {
      
      Point3D pntt;
      
      #define scheduleMe( shift1, shift2, shift3 ) \
        pntt = pnt + Point3D(shift1, shift2, shift3); \
        if ( pntt.inVolume(g.wvol.shape() ) && \
           ! ( g.wvol(pntt) & SCHEDULED ) ) { \
        add_to_schedule.push_back(pntt); \
        g.wvol(pntt) |= SCHEDULED; \
      }      
   
      scheduleMe( 1, 0, 0);
      scheduleMe(-1, 0, 0);
      scheduleMe( 0, 1, 0);
      scheduleMe( 0,-1, 0);
      scheduleMe( 0, 0, 1);
      scheduleMe( 0, 0,-1);

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


      const vector<Point3D> & tomark = g.markPoints( spn(pnt, 1, 0, 0), spn(pnt,-1, 0, 0),
                                                   spn(pnt, 0, 1, 0), spn(pnt, 0,-1, 0),
                                                   spn(pnt, 0, 0, 1), spn(pnt, 0, 0,-1) );  
      vector<Point3D>::const_iterator it = tomark.begin();  
      while ( it != tomark.end() ) {
        Point3D tpnt( (*it++) + pnt );
        if ( tpnt.inVolume(g.ivol.shape()) )
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
  
  Point3D pnt; 
  while ( dist->distribute(&pnt) && ! g.stopProcessing )    
    dist->collect( pnt, checkMe(pnt) ); 
  
  if (g.beverbose)
    prdn("Processing thread finished.");
    
}


void operate_wvol( void *(*start_routine) (void *) ) {
  
  //g.wvol = 0; // can be time consuming especially if mmaped. 
  // instead I do it in multithreading:
  
  vector<pthread_t> threads(g.run_threads);
  
  for (int ith = 0 ; ith < g.run_threads ; ith++)
    if ( pthread_create( & threads[ith], NULL, start_routine, (void*) ith ) ) // dirty hack with (void*)
      throw_error("Zeroing volume", "Can't create thread.");
  for (int ith = 0 ; ith < threads.size() ; ith++)
    pthread_join( threads[ith], 0);
  
}


void apply_mask(bool invert) {
  
  vector<pthread_t> threads(g.run_threads);
  vector<ApplyArgs> args(g.run_threads);
  
  for (int ith = 0 ; ith < g.run_threads ; ith++) {
    args[ith].invert = invert;
    args[ith].threadnum = ith;
    if ( pthread_create( & threads[ith], NULL, in_apply_thread, &(args[ith]) ) ) // dirty hack with (void*)
      throw_error("Apply mask", "Can't create thread.");
  }
  for (int ith = 0 ; ith < threads.size() ; ith++)
    pthread_join( threads[ith], 0);
   
}


void prepare_shift_volumes() {
  
    // prepare the shift-volumes  
  
  const int maxrad = max(g.radius, g.radiusM);
  const int radius2 = g.radius * g.radius;
  const int radiusM2 = g.radiusM * g.radiusM;
  
  g.newPoints = vector<Point3D>();
  if ( g.radius == g.radiusM )
    g.markPoints.reference( g.newPoints );
  else 
    g.markPoints = vector<Point3D>();
  
  for ( int x = -maxrad ; x <= maxrad ; x++ )
    for ( int y = -maxrad ; y <= maxrad ; y++ )
      for ( int z = -maxrad ; z <= maxrad ; z++ ) {
        
        const Point3D pnt(x,y,z);
        
        if ( pnt.r2() <= radius2 )

          for ( long _P00 = 0 ; _P00 <= 1 ; _P00++)
            if ( ! _P00  || ( pnt - Point3D( 1, 0, 0) ).r2() > radius2 )
              for ( long _M00 = 0 ; _M00 <= 1 ; _M00++)
                if ( ! _M00  || ( pnt - Point3D(-1, 0, 0) ).r2() > radius2 )
                  for ( long _0P0 = 0 ; _0P0 <= 1 ; _0P0++)
                    if ( ! _0P0  || ( pnt - Point3D( 0, 1, 0) ).r2() > radius2 )
                      for ( long _0M0 = 0 ; _0M0 <= 1 ; _0M0++)
                        if ( ! _0M0  || ( pnt - Point3D( 0,-1, 0) ).r2() > radius2 )
                          for ( long _00P = 0 ; _00P <= 1 ; _00P++)
                            if ( ! _00P  || ( pnt - Point3D( 0, 0, 1) ).r2() > radius2 )                            
                              for ( long _00M = 0 ; _00M <= 1 ; _00M++)
                                if ( ! _00M  || ( pnt - Point3D( 0, 0,-1) ).r2() > radius2 )
                                  g.newPoints( _P00, _M00, _0P0, _0M0, _00P, _00M ).push_back(pnt);
                    
        if ( g.radius != g.radiusM  &&  pnt.r2() <= radiusM2 )

          for ( long _P00 = 0 ; _P00 <= 1 ; _P00++)
            if ( ! _P00  || ( pnt - Point3D( 1, 0, 0) ).r2() > radiusM2 )
              for ( long _M00 = 0 ; _M00 <= 1 ; _M00++)
                if ( ! _M00  || ( pnt - Point3D(-1, 0, 0) ).r2() > radiusM2 )
                  for ( long _0P0 = 0 ; _0P0 <= 1 ; _0P0++)
                    if ( ! _0P0  || ( pnt - Point3D( 0, 1, 0) ).r2() > radiusM2 )
                      for ( long _0M0 = 0 ; _0M0 <= 1 ; _0M0++)
                        if ( ! _0M0  || ( pnt - Point3D( 0,-1, 0) ).r2() > radiusM2 )
                          for ( long _00P = 0 ; _00P <= 1 ; _00P++)
                            if ( ! _00P  || ( pnt - Point3D( 0, 0, 1) ).r2() > radiusM2 )                            
                              for ( long _00M = 0 ; _00M <= 1 ; _00M++)
                                if ( ! _00M  || ( pnt - Point3D( 0, 0,-1) ).r2() > radiusM2 )
                                  g.markPoints( _P00, _M00, _0P0, _0M0, _00P, _00M ).push_back(pnt);

      }     
      
  
}



void read_input() {
  
  vector<pthread_t> threads(g.run_threads);  
  
  // Read images to the array.
  ReadWriteDistributor rwdist;
  rwdist.PrepareForRead();
  for (int ith = 0 ; ith < g.run_threads ; ith++)
    if ( pthread_create( & threads[ith], NULL, in_read_thread, &rwdist ) )
      throw_error("Read volume", "Can't create thread.");
  for (int ith = 0 ; ith < threads.size() ; ith++)
    pthread_join( threads[ith], 0);
  
}


void save_results() {
  
  vector<pthread_t> threads(g.run_threads);   
  ReadWriteDistributor rwdist;
  rwdist.PrepareForWrite();
  
  for (int ith = 0 ; ith < g.run_threads ; ith++)
    if ( pthread_create( & threads[ith], NULL, in_write_thread, &rwdist ) )
      throw_error("Write volume", "Can't create thread.");
    for (int ith = 0 ; ith < threads.size() ; ith++)
      pthread_join( threads[ith], 0);
  
}


ProcDistributor procdist;
vector<pthread_t> proc_threads;


void start_process() {
  
  if (proc_threads.size())
    throw_error("Process data", "Another process is already running.");
  proc_threads.resize(g.run_threads);
  
  g.cross_point = g.start.front();
  
  for ( int cpnt=0 ; cpnt < g.start.size() ; cpnt++ )
    if ( ! g.start[cpnt].inVolume( g.ivol.shape() ) )
      throw_error("Process data", "At least one starting point is outside the volume.");
  
  for ( int cpnt=0 ; cpnt < g.stop.size() ; cpnt++ ) {
    int rs2 = g.stop[cpnt].radius * g.stop[cpnt].radius;
    for ( int x = - g.stop[cpnt].radius ; x <= g.stop[cpnt].radius ; x++ )
      for ( int y = - g.stop[cpnt].radius ; y <= g.stop[cpnt].radius ; y++ )
        for ( int z = -g.stop[cpnt].radius ; z <= g.stop[cpnt].radius ; z++ ) { 
          const Point3D pnt = Point3D(x,y,z) + g.stop[cpnt].center;
          if ( Point3D(x,y,z).r2() <= rs2  &&  pnt.inVolume(g.ivol.shape()) )
            g.wvol(pnt) |= ISBAD | SCHEDULED;
        }
  }   
  
  prdn(0);
  // process
  procdist = ProcDistributor();
  for (int ith = 0 ; ith < g.run_threads ; ith++)
    if ( pthread_create( & proc_threads[ith], NULL, in_proc_thread, &procdist ) )
      throw_error("Process data", "Can't create thread.");
  
}

void update_process() {
  if ( proc_threads.empty() )
    throw_error("Process data", "No process are running to finish.");
  printf( "%s", procdist.bar.print_line().c_str() );
}

void finish_process() {
  if ( proc_threads.empty() )
    throw_error("Process data", "No process are running to finish.");
  for (int ith = 0 ; ith < proc_threads.size() ; ith++)
    pthread_join( proc_threads[ith], 0);
  procdist.finilizeProgressBar();
  prdn(1);
  proc_threads.clear();
}

void process() {
  start_process();
  finish_process();  
}






#include <signal.h>


int main(int argc, char *argv[]) {

  clargs args(argc, argv);
  
  struct sigaction act1;
  act1.sa_handler = sig1handler;
  sigaction(SIGUSR1, &act1, NULL);
  
  struct sigaction act2;
  act2.sa_handler = sig2handler;
  sigaction(SIGUSR2, &act2, NULL);

  if (g.beverbose) 
    printf("My PID is %i. You can send me USR1 signal to stop filling and output whatever"
           " has been done so far; or USR2 to save current"
           " mask slices through the first start point.\n", getpid() );

  const Shape2D imgsh = ImageSizes(g.inlist[0]);
  const int nx=imgsh(0), ny=imgsh(1), nz=g.inlist.size();
  const Shape3D vshape(nx, ny, nz);
  
  prepare_shift_volumes();
      
  const int imapfile = allocateBigVolume( vshape, g.ivol);
  const int wmapfile = allocateBigVolume( vshape, g.wvol);
  
  
  read_input();
  
  operate_wvol( in_zeroing_thread ); 

  
  if ( g.interactive ) {
    
    const bool beverbose = g.beverbose;
    char* lns;
    char** aargv;
    int aargc;

    while ( lns = readline("enter a command > ") ) { 
      
      if (lns && *lns)
        add_history (lns);
      
      try {
      
        if ( string_to_argv(lns, &aargc, &aargv )  || ! aargc ) {
        
          if ( ! proc_threads.size() )
            printf("  No process is running.\n");
          else
            update_process();
            
          printf ( "  Possible commands are: start, stop, apply, save, slice.\n" );
      
        } else if ( string(aargv[0]) == "save" ) {
          
          if ( proc_threads.size() )
            throw_error("save", "A process is already running. Stop it first.");
        
          poptmx::OptionTable table = args.save_table; 
          if ( table.parse(aargc, aargv)  &&  clargs::check_save(table) )
            save_results();
          else
            table.usage();        
        
        } else if ( string(aargv[0]) == "apply" ) {
          
          if ( proc_threads.size() )
            throw_error("apply", "A process is already running. Stop it first.");
        
          g.color=0;
          bool invert(false);
          poptmx::OptionTable table = args.color_table; 
          table.add(poptmx::OPTION, &invert, 'i', "invert", "Inverts mask", "");
          table.parse(aargc, aargv);
          
          apply_mask(invert);
          
            
        } else if ( string(aargv[0]) == "reset" ) { 
          
          if ( proc_threads.size() )
            throw_error("reset", "A process is already running. Stop it first.");
         
          operate_wvol( in_zeroing_thread );          
          
        } else if ( string(aargv[0]) == "start" ) { 
          
          if ( proc_threads.size() )
            throw_error("start", "A process is already running. Stop it first.");
                    
          g.radius = 0;
          g.radiusM = 0;
          g.color = 0;
          g.minval = -1;
          g.maxval = -1;
          g.minggrad = -1; 
          g.maxggrad = -1;
          g.stopProcessing = false;
          g.start.clear();
          g.stop.clear();        
          
          poptmx::OptionTable table = args.proc_table; 
          table.parse(aargc, aargv);   
          clargs::check_proc(table);
          
          operate_wvol( in_wiping_thread ); 
          prepare_shift_volumes();
          
          g.beverbose = false;
          start_process();  
        
        } else if ( string(aargv[0]) == "stop" ) { 

          g.beverbose = beverbose;
          g.stopProcessing = true;
          finish_process();
        
        } else if ( string(aargv[0]) == "slice" ) {        
        
          Point3D cross;
          bool original=false;
          poptmx::OptionTable table; 
          table
            .add(poptmx::ARGUMENT, &cross, "point", "slices through this point", "")
            .add(poptmx::OPTION, &original, 'o', "original",
                "Save slices through the data volume, not mask.", "");
          
          if ( ! cross.inVolume(vshape) )
            cross = g.cross_point;
          
          table.parse(aargc, aargv);
        
          const Volume8U & voltosave = original ? g.ivol : g.wvol;
          Map8U yz( voltosave.shape()(2), voltosave.shape()(1) );
          yz = voltosave( cross.x(), Range::all(), Range::all() ).transpose(secondDim, firstDim) ;
          SaveImage(".yz.tif", yz);
          Map8U xz( voltosave.shape()(2), voltosave.shape()(0) );
          xz = voltosave( Range::all(), cross.y(), Range::all() ).transpose(secondDim, firstDim) ;
          SaveImage(".xz.tif", xz);
          Map8U xy = voltosave( Range::all(), Range::all(), cross.z() ).copy();
          SaveImage(".xy.tif", xy);
          

        } else {
        
          printf ( "  Unknown command \"%s\". Possible commands are: start, stop, apply, save, slice.\n", aargv[0] );
        
        }
        
      }
      catch (FloodErr err) {}
      catch (poptmx::Err err) {}
      
      free(lns); 
      free(aargv);
      
    }
     
  } else {
      
    prepare_shift_volumes();
    process();  
    save_results();
     
  }
      
  if (imapfile)
    close(imapfile);
  if (wmapfile)
    close(wmapfile);
    
}
