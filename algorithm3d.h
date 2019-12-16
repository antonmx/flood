
#ifndef _H_ALG3D_H_
#define _H_ALG3D_H_

#include "tools.h"
#include <queue>


class ProcDistributor {

private:

  Volume8U & ivol;
  Volume8U & wvol;
  const Shape3D volsh;

  const int radius;
  const int radiusM;
  const uint8_t minval;
  const uint8_t maxval;

  blitz::Array < std::vector<Point3D>, 6 > newPoints;
  blitz::Array < std::vector<Point3D>, 6 > markPoints;

  std::queue<Point3D> schedule;
  int thinwork;
  ulong run_threads;
  pthread_mutex_t proglock;
  pthread_mutex_t picklock;
  pthread_cond_t check_again;
  std::vector<pthread_t> proc_threads;
  const size_t updateBlock;
  //Time::duration tMon;
  // // Time::time_point nowSt = Time::now();
  // // tMon += Time::now() - nowSt;
  // // prdn(toString(chrono::nanoseconds(tMon ).count()/1000000000.0));


  int checkMe(const Point3D & pnt, const blitz::TinyVector<long int, 6> & spnv);

  static void * in_proc_thread (void * _thread_args);
  bool distribute( std::queue<Point3D> & pnts );
  void collect( std::queue<Point3D> & pnts);

public:

  ProgressBar bar;
  static bool stopProcessing;

  ProcDistributor(Volume8U & _ivol, Volume8U & _wvol,
                  const std::vector< Point3D > & _schedule,
                  const std::vector<Sphere> & fence,
                  int _radius, int _radiusM, int _minval, int _maxval,
                  ulong _run_threads, bool verbose);

  void start_process();
  void update_process();
  void finish_process();
  void process();
  int progress() const {return proc_threads.size() ?  bar.progress() + 1 : 0 ;}

};


#endif // _H_ALG3D_H_

