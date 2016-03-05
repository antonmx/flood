#include "tools.h"



struct GlobalVariables {
  
  Path command;               ///< Command name as it was invoked.
  bool interactive;
  std::vector<Path> inlist;        ///< List of all input images.
  Path out_filled;       ///< The mask for the output file names.
  Path out_inverted;
  Path out_mask;
  bool beverbose;       ///< Be verbose flag
  std::vector<Point3D> start;          ///< Point to start the fill at.
  std::vector<Sphere> stop;           ///< Point to stop at.
//  int nextdim;          ///< dimensions for the neibour.
  
  const int run_threads;
  
  Volume8U ivol;
  Volume8U wvol;
  
  // Parameters to use in the check function. 
  int radius;
  int radiusM;
  unsigned int color;
  int minval;
  int maxval;
  int minggrad; // square of the grad!
  int maxggrad; // square of the grad!
  
  blitz::Array < std::vector<Point3D>, 6 > newPoints;
  blitz::Array < std::vector<Point3D>, 6 > markPoints;
  
  Point3D cross_point;
  
  bool stopProcessing;
  
  GlobalVariables() :
    interactive(false),
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
    newPoints (2,2,2,2,2,2),
    markPoints(2,2,2,2,2,2),
    stopProcessing(false)
  {}
  
};

extern GlobalVariables g;







class ReadWriteDistributor {
  
private:
  
  long int currentidx;
  static pthread_mutex_t iolock;
  static pthread_mutex_t proglock;
  ProgressBar bar;
  
public:
  
  ReadWriteDistributor() {}
    
  void PrepareForRead(){
    bar = ProgressBar( g.beverbose , "reading volume", g.ivol.extent(blitz::firstDim) );
    currentidx=0;
  }
  
  void PrepareForWrite(){
    bar = ProgressBar( g.beverbose, "writing results", g.ivol.extent(blitz::firstDim) );
    currentidx=0;
  }

    
  bool distribute( long int * idx ) {
    pthread_mutex_lock( & iolock );
    *idx = currentidx++;
    pthread_mutex_unlock( & iolock );
    return *idx < g.ivol.extent(blitz::firstDim) ;
  }
  
  void updateProg() {
    pthread_mutex_lock( & proglock );
    bar.update();
    pthread_mutex_unlock( & proglock );    
  }
  
};


void * in_read_thread (void * _thread_args);

void * in_write_thread (void * _thread_args);


int allocateBigVolume( const Shape3D & shape, Volume8U & data );

void * in_zeroing_thread (void * _thread_args);
void * in_wiping_thread (void * _thread_args);


struct ApplyArgs {
  bool invert;
  int threadnum;  
};

void * in_apply_thread (void * _thread_args);







/// \CLARGS
struct clargs {
  
  poptmx::OptionTable proc_table;
  poptmx::OptionTable save_table;
  poptmx::OptionTable color_table;
  
  poptmx::OptionTable table;
  
  /// \CLARGSF
  clargs(int argc, char *argv[]);
  
  static void check_proc( poptmx::OptionTable & tab );
  static void check_save( poptmx::OptionTable & tab );
  
};



bool string_to_argv(char const * str, int * argc_p, char *** argv_p);


void sig1handler(int signum);

void sig2handler(int signum);






