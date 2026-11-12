#include "local_checks.h"
int main(int argc,char** argv) {
  int rank=0;
#ifdef LOCAL_VORONOI_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
  int success=1;
  std::ofstream output(rank==0 ? "output" : "/dev/null");
  output<<std::unitbuf;
  try {
    runChecks(output);
  } catch(const std::exception& error) {
    output<<"FAIL "<<error.what()<<"\n";
    std::cerr<<error.what()<<"\n";
    success=0;
  }
#ifdef LOCAL_VORONOI_MPI
  int globalSuccess=0;
  MPI_Allreduce(&success,&globalSuccess,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  success=globalSuccess;
  MPI_Finalize();
#endif
  return success ? 0 : 1;
}
