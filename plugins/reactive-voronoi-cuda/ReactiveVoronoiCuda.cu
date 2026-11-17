// FP64 streamed soft assignment: O(Ncenters+Nassigned) device storage.
// Build: nvcc -O3 -std=c++17 -arch=sm_70 -Xcompiler=-fPIC -shared ReactiveVoronoiCuda.cu -o libReactiveVoronoiCuda.so
#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <new>

namespace {
constexpr int blockSize=128;
struct State {
  int nc,nh;
  double *x=nullptr,*occupancy=nullptr,*maximum=nullptr,*normalization=nullptr;
  double *sensitivity=nullptr,*derivatives=nullptr,*virial=nullptr;
  int* invalid=nullptr;
  double box[3]{},kappa=0;
  bool pbc=false;
  cudaStream_t stream=nullptr;
};
__device__ double reduce(double value,bool maximum,double* scratch) {
  const int t=threadIdx.x;
  scratch[t]=value; __syncthreads();
  for(int offset=blockSize/2;offset;offset/=2) {
    if(t<offset) scratch[t]=maximum ? fmax(scratch[t],scratch[t+offset]) : scratch[t]+scratch[t+offset];
    __syncthreads();
  }
  const double result=scratch[0]; __syncthreads(); return result;
}
__device__ double distance(const State s,int i,int j,double* d) {
  double r2=0;
  for(int k=0;k<3;++k) {
    d[k]=s.x[3*(s.nc+j)+k]-s.x[3*i+k];
    if(s.pbc) d[k]-=nearbyint(d[k]/s.box[k])*s.box[k];
    r2+=d[k]*d[k];
  }
  return sqrt(r2);
}
__global__ void forwardKernel(State s) {
  __shared__ double scratch[blockSize];
  const int j=blockIdx.x,t=threadIdx.x;
  double maximum=-INFINITY;
  for(int i=t;i<s.nc;i+=blockSize) {
    double d[3]; const double r=distance(s,i,j,d), score=-s.kappa*r;
    if(!isfinite(r)||r<=2.220446049250313e-16||!isfinite(score)) atomicExch(s.invalid,1);
    maximum=fmax(maximum,score);
  }
  maximum=reduce(maximum,true,scratch);
  double norm=0;
  for(int i=t;i<s.nc;i+=blockSize) {
    double d[3]; norm+=exp(-s.kappa*distance(s,i,j,d)-maximum);
  }
  norm=reduce(norm,false,scratch);
  if(t==0) {s.maximum[j]=maximum;s.normalization[j]=norm;}
  for(int i=t;i<s.nc;i+=blockSize) {
    double d[3]; const double w=exp(-s.kappa*distance(s,i,j,d)-maximum)/norm;
    atomicAdd(s.occupancy+i,w);
  }
}
__global__ void backwardKernel(State s) {
  __shared__ double scratch[blockSize];
  const int j=blockIdx.x,t=threadIdx.x;
  const double maximum=s.maximum[j],norm=s.normalization[j];
  double mean=0;
  for(int i=t;i<s.nc;i+=blockSize) {
    double d[3]; const double w=exp(-s.kappa*distance(s,i,j,d)-maximum)/norm;
    mean+=w*s.sensitivity[i];
  }
  mean=reduce(mean,false,scratch);
  double assigned[3]={0,0,0},virial[9]={};
  for(int i=t;i<s.nc;i+=blockSize) {
    double d[3]; const double r=distance(s,i,j,d);
    const double w=exp(-s.kappa*r-maximum)/norm;
    const double radial=-s.kappa*w*(s.sensitivity[i]-mean)/r;
    for(int k=0;k<3;++k) {
      const double grad=radial*d[k];
      atomicAdd(s.derivatives+3*i+k,-grad);
      assigned[k]+=grad;
      for(int l=0;l<3;++l) virial[3*k+l]-=grad*d[l];
    }
  }
  for(int k=0;k<3;++k) {
    const double value=reduce(assigned[k],false,scratch);
    if(t==0) s.derivatives[3*(s.nc+j)+k]=value;
  }
  for(int k=0;k<9;++k) {
    const double value=reduce(virial[k],false,scratch);
    if(t==0) atomicAdd(s.virial+k,value);
  }
}
bool ok(cudaError_t status) {
  if(status==cudaSuccess) return true;
  std::fprintf(stderr,"ReactiveVoronoi CUDA: %s\n",cudaGetErrorString(status)); return false;
}
void release(State* s) {
  if(!s) return;
  cudaFree(s->x);cudaFree(s->occupancy);cudaFree(s->maximum);cudaFree(s->normalization);
  cudaFree(s->sensitivity);cudaFree(s->derivatives);cudaFree(s->virial);cudaFree(s->invalid);
  if(s->stream) cudaStreamDestroy(s->stream);
  delete s;
}
}
extern "C" void* reactive_voronoi_create_v1(int nc,int nh) {
  if(nc<=0||nh<=0) return nullptr;
  State* s=new(std::nothrow) State; if(!s) return nullptr; s->nc=nc;s->nh=nh;
  if(!ok(cudaStreamCreateWithFlags(&s->stream,cudaStreamNonBlocking)) ||
     !ok(cudaMalloc(&s->x,3ULL*(nc+nh)*sizeof(double))) ||
     !ok(cudaMalloc(&s->occupancy,nc*sizeof(double))) ||
     !ok(cudaMalloc(&s->maximum,nh*sizeof(double))) ||
     !ok(cudaMalloc(&s->normalization,nh*sizeof(double))) ||
     !ok(cudaMalloc(&s->sensitivity,nc*sizeof(double))) ||
     !ok(cudaMalloc(&s->derivatives,3ULL*(nc+nh)*sizeof(double))) ||
     !ok(cudaMalloc(&s->virial,9*sizeof(double))) ||
     !ok(cudaMalloc(&s->invalid,sizeof(int)))) {release(s);return nullptr;}
  return s;
}
extern "C" void reactive_voronoi_destroy_v1(void* p) {release(static_cast<State*>(p));}
extern "C" int reactive_voronoi_forward_v1(void* p,const double* x,const double* box,int pbc,double kappa,double* occupancy) {
  State& s=*static_cast<State*>(p); s.pbc=pbc;s.kappa=kappa;
  for(int k=0;k<3;++k) s.box[k]=box[k];
  if(!ok(cudaMemcpyAsync(s.x,x,3ULL*(s.nc+s.nh)*sizeof(double),cudaMemcpyHostToDevice,s.stream)) ||
     !ok(cudaMemsetAsync(s.occupancy,0,s.nc*sizeof(double),s.stream)) ||
     !ok(cudaMemsetAsync(s.invalid,0,sizeof(int),s.stream))) return 1;
  forwardKernel<<<s.nh,blockSize,0,s.stream>>>(s);
  int invalid=0;
  if(!ok(cudaGetLastError()) ||
     !ok(cudaMemcpyAsync(occupancy,s.occupancy,s.nc*sizeof(double),cudaMemcpyDeviceToHost,s.stream)) ||
     !ok(cudaMemcpyAsync(&invalid,s.invalid,sizeof(int),cudaMemcpyDeviceToHost,s.stream)) ||
     !ok(cudaStreamSynchronize(s.stream))) return 1;
  return invalid;
}
extern "C" int reactive_voronoi_backward_v1(void* p,const double* sensitivity,double* derivatives,double* virial) {
  State& s=*static_cast<State*>(p);
  if(!ok(cudaMemcpyAsync(s.sensitivity,sensitivity,s.nc*sizeof(double),cudaMemcpyHostToDevice,s.stream)) ||
     !ok(cudaMemsetAsync(s.derivatives,0,3ULL*(s.nc+s.nh)*sizeof(double),s.stream)) ||
     !ok(cudaMemsetAsync(s.virial,0,9*sizeof(double),s.stream))) return 1;
  backwardKernel<<<s.nh,blockSize,0,s.stream>>>(s);
  if(!ok(cudaGetLastError()) ||
     !ok(cudaMemcpyAsync(derivatives,s.derivatives,3ULL*(s.nc+s.nh)*sizeof(double),cudaMemcpyDeviceToHost,s.stream)) ||
     !ok(cudaMemcpyAsync(virial,s.virial,9*sizeof(double),cudaMemcpyDeviceToHost,s.stream)) ||
     !ok(cudaStreamSynchronize(s.stream))) return 1;
  return 0;
}
