
#ifndef _H_ALG3D_H_
#define _H_ALG3D_H_

#include "tools.h"
#include <queue>


template<class ItemType> class Fqueue {

private:

  struct Fitem {
    ItemType val;
    Fitem * next;
    Fitem( ItemType _val, Fitem * _next) : val(_val), next(_next) {}
  };

  size_t sz;
  Fitem * first;
  Fitem * last;

  Fqueue(size_t _sz, Fitem * _first, Fitem * _last)
    : sz(_sz)
    , first(_first)
    , last(_last)
  {
    if (sz)
      last->next=0;
  }

public:

  Fqueue()
   : Fqueue(0, 0, 0)
  {}

  ~Fqueue() {
    while(sz--) {
      Fitem * del = first;
      first = del->next;
      delete del;
    }
  }

  const size_t & size() const { return sz; }

  Fqueue & push(const ItemType& el) {
    Fitem * ni = new Fitem(el, 0);
    if (sz)
      last->next = ni;
    last = ni;
    if (!sz)
      first = last;
    sz++;
    return * this;
  }

  Fqueue & push(Fqueue & addMe) {
    if ( ! addMe.sz )
      return * this;
    if(sz)
      last->next = addMe.first;
    last = addMe.last;
    if (!sz)
      first = addMe.first;
    sz += addMe.sz;
    addMe = Fqueue();
    return * this;
  }

  ItemType pop() {
    if (!sz)
      return ItemType();
    Fitem * del = first;
        ItemType ret = del->val;
    first = del->next;
    last->next = first;
    delete del;
    sz--;
    return ret;
  }

  Fqueue & pop(const size_t nit) {

    if ( ! nit  ||  ! sz )
      return *this;
    if (nit>=sz) {
      Fqueue * nq =new Fqueue(sz, first, last);
      *this = Fqueue();
      return *nq;
    }

    Fitem * cur = first;
    for (int idx=1 ; idx<nit ; idx++)
      cur = cur->next;
    last->next = cur->next;
    Fqueue * nq = new Fqueue(nit, first, cur);
    first = last->next;
    last->next = 0;
    sz -= nit;
    return *nq;

  }

};



class ProcDistributor {

private:

  Volume8U & ivol;
  Volume8U & wvol;
  const Shape3D volsh;

  Fqueue<Point3D> schedule;
  int thinwork;
  ulong run_threads;

  const int radius;
  const int radiusM;
  const uint8_t minval;
  const uint8_t maxval;

  blitz::Array < std::vector<Point3D>, 6 > newPoints;
  blitz::Array < std::vector<Point3D>, 6 > markPoints;

  pthread_mutex_t proglock;
  pthread_mutex_t picklock;
  pthread_cond_t check_again;
  std::vector<pthread_t> proc_threads;
  Time::duration tMon;
  Time::duration tMon2;
  // // Time::time_point nowSt = Time::now();
  // // tMon += Time::now() - nowSt;
  // // prdn(toString(chrono::nanoseconds(tMon ).count()/1000000000.0));


  int checkMe(const Point3D & pnt, const blitz::TinyVector<long int, 6> & spnv);

  static void * in_proc_thread (void * _thread_args);
  bool distribute( Fqueue<Point3D> & pnts );
  void collect( Fqueue<Point3D> & pnts);

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

