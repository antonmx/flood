
#include "tools3d.h"

using namespace blitz;
using namespace std;







clargs::clargs(int argc, char *argv[]) :
  table("3D advanced flood fill.",
        "Advanced flood fill algorithm implemented for the 3D volume.")
{

  table
    .add(poptmx::NOTE, "ARGUMENTS:")
    .add(poptmx::ARGUMENT, &g.inlist, "list", "List of the input images.", "")

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
  if ( ! table.count(&g.inlist) )
    exit_on_error(g.command, "Missing required argument: "+table.desc(&g.inlist)+".");
  
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
  
  g.cross_point = g.start.front();
      
}



#include <signal.h>

void sig1handler(int signum) {
  g.stopProcessing=true;
}

void sig2handler(int signum) {
  
  Point3D pnt;  
  char line[256];
  
  FILE *funcf = fopen( ".sliceme.txt", "r" );
  
  if ( ! funcf || 
       1 != fscanf ( funcf,  "%254[^\n]\n", line) ||
       1 != _conversion (&pnt, line) ||
       ! pnt.inVolume(g.wvol.shape()) ) 
    pnt = g.cross_point;
  
  if (funcf)
    fclose(funcf);
  
  Map8U yz( g.wvol.shape()(2), g.wvol.shape()(1) );
  yz = g.wvol( pnt.x(), Range::all(), Range::all() ).transpose(secondDim, firstDim) ;
  SaveImage(".yz.tif", yz);
  Map8U xz( g.wvol.shape()(2), g.wvol.shape()(0) );
  xz = g.wvol( Range::all(), pnt.y(), Range::all() ).transpose(secondDim, firstDim) ;
  SaveImage(".xz.tif", xz);
  Map8U xy = g.wvol( Range::all(), Range::all(), pnt.z() ).copy();
  SaveImage(".xy.tif", xy);
  
  printf("Written slices through (%i, %i, %i) point."
         " Files are named \".xy.tif\", \".xz.tif\" and \".yz.tif\".\n",
         pnt.x(), pnt.y(), pnt.z() );
  
}


void * in_zeroing_thread (void * _thread_args) {
  
  int threadnum = (long int) _thread_args; // dirty hack
    
  const int thichness =  g.run_threads ?
      1 + g.wvol.shape()(2) / g.run_threads :
      g.wvol.shape()(2);
  Range zRange( threadnum*thichness,
                min( (threadnum+1)*thichness - 1 , int( g.wvol.shape()(2) ) ) );
  g.wvol( Range::all(), Range::all(),  zRange) = 0;
  
}




pthread_mutex_t ReadWriteDistributor::iolock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t ReadWriteDistributor::proglock = PTHREAD_MUTEX_INITIALIZER;



void * in_read_thread (void * _thread_args) {
  
  ReadWriteDistributor *  dist = (ReadWriteDistributor*) _thread_args;
  if (!dist)
    throw_error("read thread", "Inappropriate thread function arguments.");
  
  Shape2D imgsh = Shape2D(g.ivol.extent(firstDim) , g.ivol.extent(secondDim) ) ;
  Map8U img(imgsh);
  
  long int idx;
  while ( dist->distribute(&idx) ) {
    ReadImage(g.inlist[idx], img, imgsh);
    g.ivol( Range::all(), Range::all(), idx ) = img;
    dist->updateProg();
  }
    
}





void * in_write_thread (void * _thread_args) {
  
  ReadWriteDistributor *  dist = (ReadWriteDistributor*) _thread_args;
  if (!dist)
    throw_error("write thread", "Inappropriate thread function arguments.");
  
  const Shape2D imgsh(g.ivol.extent(firstDim) , g.ivol.extent(secondDim));
  Map8U img(imgsh), wmg(imgsh);
  Path outPath;
  
  long int idx;
  while ( dist->distribute(&idx) ) {
    
    wmg = g.wvol( Range::all(), Range::all(), idx );
    
    if ( ! g.out_filled.empty() ) {
      for ( long int sh0 = 0 ; sh0 < imgsh(0) ; sh0++)
        for ( long int sh1 = 0 ; sh1 < imgsh(1) ; sh1++)
          img(sh0, sh1) =  ( wmg(sh0, sh1) & FILLED ) ? g.ivol( sh0, sh1, idx ) : g.color ;
      outPath = g.out_filled + g.inlist[idx].name();
      SaveImage(outPath, img);
    }
    
    if ( ! g.out_inverted.empty() ) {
      for ( long int sh0 = 0 ; sh0 < imgsh(0) ; sh0++)
        for ( long int sh1 = 0 ; sh1 < imgsh(1) ; sh1++)
          img(sh0, sh1) =  ( wmg(sh0, sh1) & FILLED ) ? g.color : g.ivol( sh0, sh1, idx ) ;
      outPath = g.out_inverted + g.inlist[idx].name();
      SaveImage(outPath, img);
    }
    
    if ( ! g.out_mask.empty() ) {
      outPath = g.out_mask + g.inlist[idx].name();
      SaveImage(outPath, wmg);
    }
        
    dist->updateProg();
    
  }
  
  
}






