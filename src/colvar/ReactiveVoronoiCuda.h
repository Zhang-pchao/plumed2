// Optional, explicitly selected CUDA plugin. CPU builds require no CUDA headers.
#ifndef PLUMED_REACTIVE_VORONOI_CUDA_H
#define PLUMED_REACTIVE_VORONOI_CUDA_H
#include <dlfcn.h>
#include <stdexcept>
#include <string>

namespace PLMD {
namespace colvar {
class ReactiveVoronoiCuda {
  void* library_=nullptr;
  void* state_=nullptr;
  using Create=void*(*)(int,int);
  using Destroy=void(*)(void*);
  using Forward=int(*)(void*,const double*,const double*,int,double,double*);
  using Backward=int(*)(void*,const double*,double*,double*);
  Destroy destroy_=nullptr;
  Forward forward_=nullptr;
  Backward backward_=nullptr;
  template<class T> T symbol(const char* name) {
    void* value=dlsym(library_,name);
    if(!value) throw std::runtime_error(std::string("Missing CUDA backend symbol: ")+name);
    return reinterpret_cast<T>(value);
  }
public:
  ReactiveVoronoiCuda(const std::string& path,int centers,int assigned) {
    library_=dlopen(path.c_str(),RTLD_NOW|RTLD_LOCAL);
    if(!library_) throw std::runtime_error(std::string("Cannot load CUDA backend: ")+dlerror());
    try {
      auto create=symbol<Create>("reactive_voronoi_create_v1");
      destroy_=symbol<Destroy>("reactive_voronoi_destroy_v1");
      forward_=symbol<Forward>("reactive_voronoi_forward_v1");
      backward_=symbol<Backward>("reactive_voronoi_backward_v1");
      state_=create(centers,assigned);
      if(!state_) throw std::runtime_error("CUDA backend allocation failed");
    } catch(...) { dlclose(library_); library_=nullptr; throw; }
  }
  ~ReactiveVoronoiCuda() { if(state_) destroy_(state_); if(library_) dlclose(library_); }
  ReactiveVoronoiCuda(const ReactiveVoronoiCuda&)=delete;
  ReactiveVoronoiCuda& operator=(const ReactiveVoronoiCuda&)=delete;
  void forward(const double* x,const double* box,bool pbc,double kappa,double* occupancy) {
    if(forward_(state_,x,box,pbc,kappa,occupancy)!=0)
      throw std::runtime_error("CUDA reactive Voronoi forward failed (invalid geometry or CUDA error)");
  }
  void backward(const double* sensitivity,double* derivatives,double* virial) {
    if(backward_(state_,sensitivity,derivatives,virial)!=0)
      throw std::runtime_error("CUDA reactive Voronoi backward failed");
  }
};
}
}
#endif
