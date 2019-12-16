
#include <pthread.h>
#include <algorithm>
#include <unistd.h>

#include "algorithm3d.h"

using namespace std;
using namespace blitz;

static const array < Point3D, 6 > crossPoints = {
  Point3D( 1, 0, 0),  Point3D(-1, 0, 0),
  Point3D( 0, 1, 0),  Point3D( 0,-1, 0),
  Point3D( 0, 0, 1),  Point3D( 0, 0,-1) };


/*
int ggradient( const Point3D & pnt, const Volume8U & vol ) {
  int sum=0, num=0;
  for( CrsIt itcs = crsBgn ; itcs < crsEnd ; itcs++ ) {
    Point3D tpnt = pnt + *itcs++;
    if( tpnt.inVolume(vol.shape()) ) {
      int gg = vol(pnt) - vol (tpnt) ;
      sum+=gg*gg;
      num++;
    }
  }
  return num ? sum/num : 0 ;
}
*/









bool ProcDistributor::stopProcessing = false;


ProcDistributor::ProcDistributor(Volume8U & _ivol, Volume8U & _wvol,
                  const std::vector< Point3D > & _schedule,
                  const std::vector<Sphere> & fence,
                  int _radius, int _radiusM, int _minval, int _maxval,
                  ulong _run_threads, bool verbose)
  : ivol(_ivol)
  , wvol(_wvol)
  , volsh(_ivol.shape())
  , schedule()
  , thinwork(0)
  , run_threads(_run_threads)
  , radius(_radius)
  , radiusM(_radiusM)
  , minval(_minval)
  , maxval(_maxval)
  , newPoints (2,2,2,2,2,2)
  , markPoints(2,2,2,2,2,2)
  , proglock(PTHREAD_MUTEX_INITIALIZER)
  , picklock(PTHREAD_MUTEX_INITIALIZER)
  , check_again(PTHREAD_COND_INITIALIZER)
  , proc_threads()
  , tMon(0)
  , tMon2(0)
  , bar( verbose , "processing volume", _ivol.size() )
{

  proc_threads.reserve(run_threads);

  if (_schedule.empty())
    throw_error("Process data", "No starting points.");
  vector<Point3D>::const_iterator it = _schedule.begin();
  while ( it != _schedule.end() ) {
    if ( it->inVolume(volsh) )
      schedule.push(*it);
    else
      warn("Process data", "At least one starting point is outside the volume.");
    it++;
  }

  for ( int cpnt=0 ; cpnt < fence.size() ; cpnt++ ) {
    int rs2 = fence[cpnt].radius * fence[cpnt].radius;
    for ( int x = - fence[cpnt].radius ; x <= fence[cpnt].radius ; x++ )
      for ( int y = - fence[cpnt].radius ; y <= fence[cpnt].radius ; y++ )
        for ( int z = -fence[cpnt].radius ; z <= fence[cpnt].radius ; z++ ) {
          const Point3D pnt = Point3D(z,y,x) + fence[cpnt].center;
          if ( Point3D(z,y,x).r2() <= rs2  &&  pnt.inVolume(ivol.shape()) )
            wvol(pnt) |= ISBAD | SCHEDULED;
        }
  }

  const int maxrad = max(radius, radiusM);
  const int radius2 = radius * radius;
  const int radiusM2 = radiusM * radiusM;

  for ( int x = -maxrad ; x <= maxrad ; x++ )
    for ( int y = -maxrad ; y <= maxrad ; y++ )
      for ( int z = -maxrad ; z <= maxrad ; z++ ) {

        const Point3D pnt(z,y,x);

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
                                  newPoints( _P00, _M00, _0P0, _0M0, _00P, _00M ).push_back(pnt);

        if ( pnt.r2() <= radiusM2 )

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
                                  markPoints( _P00, _M00, _0P0, _0M0, _00P, _00M ).push_back(pnt);

      }

}


bool ProcDistributor::distribute( Fqueue<Point3D> & pnts ) {

  pthread_mutex_lock( & picklock );
  const Time::time_point nowSt = Time::now();
  while ( thinwork && ! schedule.size() )
    pthread_cond_wait(&check_again, &picklock);
  if (schedule.size()) {
    pnts = schedule.pop( std::max(schedule.size()/run_threads, ulong(1)) );
    thinwork++;
  }
  tMon += Time::now() - nowSt;
  pthread_mutex_unlock( & picklock );
  pthread_cond_signal( & check_again ) ; // do i need it here?

  return pnts.size();

}


void ProcDistributor::collect( Fqueue<Point3D> & add_to_schedule) {
  pthread_mutex_lock( & picklock );
  schedule.push(add_to_schedule);
  thinwork--;
  pthread_mutex_unlock( & picklock );
  pthread_cond_signal( & check_again ) ;
}


int ProcDistributor::checkMe(const Point3D & pnt, const blitz::TinyVector<long int, 6> & spnv) {

  const vector<Point3D> & tocheck = newPoints(spnv);
  vector<Point3D>::const_iterator it = tocheck.begin();
  while ( it != tocheck.end() )  {

    Point3D tpnt( (*it++) + pnt );
    uint8_t & wal=wvol(tpnt);

    if ( wal & ISBAD )
      return 0;

    if ( ! ( wal & ISGOOD ) ) {
      const uint8_t & ial=ivol(tpnt);
      const bool pass = ( minval < 0  ||  ial >= minval ) &&
                        ( maxval < 0  ||  ial <= maxval ) ;
      wal |= pass ? ISGOOD : ISBAD ; // can I do it in multithreaded?
      if (!pass)
        return 0;
    }

  }
  return 1;
}


void * ProcDistributor::in_proc_thread (void * _thread_args) {

  ProcDistributor *  dist = (ProcDistributor*) _thread_args;
  if (!dist)
    throw_error("process thread", "Inappropriate thread function arguments.");

  Fqueue<Point3D> pnts;
  while ( dist->distribute(pnts)  &&  ! ProcDistributor::stopProcessing ) {

    size_t cnt = 0;
    do {

      Point3D pnt = pnts.pop();
      #define spn(a1, a2, a3) ( dist->wvol(pnt + Point3D(a1, a2, a3) ) & CHECKED ?  1l  :  0l )
      const TinyVector<long int, 6> spnv( spn(1, 0, 0), spn(-1, 0, 0),
                                          spn(0, 1, 0), spn( 0,-1, 0),
                                          spn(0, 0, 1), spn( 0, 0,-1) );
      #undef spn

      if ( dist->checkMe(pnt, spnv) ) {

        array<Point3D,6>::const_iterator itcs = crossPoints.begin();
        while(itcs < crossPoints.end()) {
          const Point3D pntt = pnt + *itcs++;
          uint8_t & wal = dist->wvol(pntt);
          if ( ! ( wal & SCHEDULED ) ) {
            pnts.push(pntt);
            cnt++;
            wal |= SCHEDULED;
           }
         }

        const vector<Point3D> & tomark = dist->markPoints(spnv);
        vector<Point3D>::const_iterator it = tomark.begin();
        while ( it != tomark.end() )
          dist->wvol(pnt + *it++) |= FILLED ;

        dist->wvol(pnt) |= CHECKED;

      }

      if (cnt * 10000 > dist->ivol.size() ) {
        pthread_mutex_lock( & dist->proglock );
        dist->bar.update( dist->bar.progress() + cnt );
        pthread_mutex_unlock( & dist->proglock );
        cnt = 0;
      }

    } while( pnts.size() &&
             pnts.size() < dist->schedule.size() );

    dist->collect(pnts);

  }

  return 0;

}

void ProcDistributor::start_process() {

  if (proc_threads.size())
    throw_error("Process data", "Another process is already running.");

  for (int ith = 0 ; ith < run_threads ; ith++) {
    pthread_t curthread;
    if ( pthread_create( & curthread, NULL, in_proc_thread, this ) )
      throw_error("Process data", "Can't create thread.");
    else
      proc_threads.push_back(curthread);
  }

}

void ProcDistributor::update_process() {
  if ( proc_threads.empty() )
    throw_error("Process data", "No process are running to finish.");
  printf( "%s", bar.print_line().c_str() );
}

void ProcDistributor::finish_process() {
  if ( proc_threads.empty() )
    throw_error("Process data", "No process are running to finish.");
  for (int ith = 0 ; ith < proc_threads.size() ; ith++)
    pthread_join( proc_threads[ith], 0);
  prdn(toString(chrono::nanoseconds(tMon).count()/1000000000.0));
  prdn(toString(chrono::nanoseconds(tMon2).count()/1000000000.0));
  bar.update();
  proc_threads.clear();
}

void ProcDistributor::process() {
  start_process();
  finish_process();
}
















