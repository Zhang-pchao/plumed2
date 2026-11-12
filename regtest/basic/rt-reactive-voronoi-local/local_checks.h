#include "plumed/wrapper/Plumed.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef LOCAL_VORONOI_MPI
#include <mpi.h>
#endif

namespace {
using Box=std::array<double,9>;
struct Result {
  double value=0;
  std::vector<double> force;
  Box virial{};
};
const Box normalBox{{2,0,0,0,2,0,0,0,2}};
const std::vector<double> fixture{
  0,0,0, 0.27,0.02,0, 0.54,0.05,0.01,
  0.095,0.02,0, -0.03,0.092,0.01, 0.13,0.015,0.03,
  0.30,0.10,0.02, 0.60,0.07,0.02, 0.49,0.08,-0.02,
  0.08,-0.05,0.03};
const std::string common=
  "v: VORONOI_LOCAL_COORDINATION CENTERS=1-3 ASSIGNED=4-9 "
  "REGION_CENTER=10 REFERENCE=2 KAPPA=12 CORE_RADIUS=0.12 OUTER_RADIUS=0.38 ";
const std::string compact="PAIR_SWITCH=0.25 PAIR_CUTOFF=0.45 ";

class Engine {
  PLMD::Plumed plumed;
  int step=0;
  bool biased;
  double monitoredValue=0;
public:
  explicit Engine(const std::string& input,bool applyBias=true,
                  const std::string& bias="b: BIASVALUE ARG=v",const std::string& setup=""):biased(applyBias) {
    int precision=8, natoms=10;
    double timestep=0.001;
#ifdef LOCAL_VORONOI_MPI
    MPI_Comm communicator=MPI_COMM_WORLD;
    plumed.cmd("setMPIComm",&communicator);
#endif
    plumed.cmd("setRealPrecision",&precision);
    plumed.cmd("setNatoms",&natoms);
    plumed.cmd("setTimestep",&timestep);
    plumed.cmd("setLogFile","/dev/null");
    plumed.cmd("init");
    if(!setup.empty()) {
      plumed.cmd("readInputLine",setup.c_str());
    }
    plumed.cmd("readInputLine",input.c_str());
    if(biased) {
      plumed.cmd("readInputLine",bias.c_str());
    } else {
      plumed.cmd("readInputLine","PRINT ARG=v FILE=monitor STRIDE=1");
      long dataRank=-1;
      plumed.cmd("getDataRank v",&dataRank);
      plumed.cmd("setMemoryForData v",&monitoredValue);
    }
  }
  Result evaluate(std::vector<double> x, Box box=normalBox) {
    Result result;
    result.force.resize(x.size(),0);
    std::vector<double> masses(x.size()/3,1.0);
    plumed.cmd("setStep",&step);
    plumed.cmd("setPositions",x.data());
    plumed.cmd("setMasses",masses.data());
    plumed.cmd("setBox",box.data());
    plumed.cmd("setForces",result.force.data());
    plumed.cmd("setVirial",result.virial.data());
    plumed.cmd("calc");
    plumed.cmd("getBias",&result.value);
    if(!biased) {
      if(result.value!=0) {
        throw std::runtime_error("monitor-only action applied bias");
      }
      result.value=monitoredValue;
    }
#ifdef LOCAL_VORONOI_MPI
    // DomainDecomposition::Set_comm suppresses Box force output on non-root
    // ranks. Compare the communicator-total host virial, not each local buffer.
    Box total{};
    MPI_Allreduce(result.virial.data(),total.data(),9,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    result.virial=total;
#endif
    ++step;
    return result;
  }
};

void require(bool condition,const std::string& name) {
#ifdef LOCAL_VORONOI_MPI
  int local=condition ? 1 : 0, global=0;
  MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  condition=global!=0;
#endif
  if(!condition) {
    throw std::runtime_error(name);
  }
}
void close(double a,double b,double tolerance,const std::string& name) {
  const bool ok=std::isfinite(a) && std::isfinite(b) && std::fabs(a-b)<=tolerance*(1+std::max(std::fabs(a),std::fabs(b)));
  if(!ok) {
    std::cerr<<name<<" actual="<<a<<" expected="<<b<<" tolerance="<<tolerance<<"\n";
  }
  require(ok,name);
}
void compare(const Result& a,const Result& b,const std::string& name) {
  close(a.value,b.value,1e-11,name+" value");
  for(unsigned i=0; i<a.force.size(); ++i) {
    close(a.force[i],b.force[i],1e-10,name+" force");
  }
  for(unsigned i=0; i<9; ++i) {
    close(a.virial[i],b.virial[i],1e-10,name+" virial");
  }
}

double exponentialOracle(const std::vector<double>& x) {
  std::array<double,3> occupancy{};
  for(unsigned j=3; j<9; ++j) {
    std::array<double,3> weights{};
    double norm=0;
    for(unsigned i=0; i<3; ++i) {
      double r2=0;
      for(unsigned k=0; k<3; ++k) {
        const double delta=x[3*i+k]-x[3*j+k];
        r2+=delta*delta;
      }
      weights[i]=std::exp(-12*std::sqrt(r2));
      norm+=weights[i];
    }
    for(unsigned i=0; i<3; ++i) {
      occupancy[i]+=weights[i]/norm;
    }
  }
  double value=0;
  for(unsigned i=0; i<3; ++i) {
    double r2=0;
    for(unsigned k=0; k<3; ++k) {
      const double delta=x[3*i+k]-x[27+k];
      r2+=delta*delta;
    }
    const double t=std::max(0.0,std::min(1.0,(std::sqrt(r2)-0.12)/0.26));
    const double chi=1-10*std::pow(t,3)+15*std::pow(t,4)-6*std::pow(t,5);
    value+=chi*(occupancy[i]-2)*(occupancy[i]-2);
  }
  return value;
}

void checkDerivatives(Engine& engine) {
  const Result base=engine.evaluate(fixture);
  for(double h: {
        1e-5,1e-6,1e-7
      }) {
    for(unsigned i=0; i<fixture.size(); ++i) {
      auto plus=fixture,minus=fixture;
      plus[i]+=h;
      minus[i]-=h;
      const double numerical=-(engine.evaluate(plus).value-engine.evaluate(minus).value)/(2*h);
      close(base.force[i],numerical,3e-6,"coordinate derivative");
    }
  }
  for(unsigned a=0; a<3; ++a) {
    for(unsigned b=0; b<3; ++b) {
      const double h=1e-6;
      auto plus=fixture,minus=fixture;
      auto bp=normalBox,bm=normalBox;
      for(unsigned i=0; i<10; ++i) {
        plus[3*i+a]+=h*fixture[3*i+b];
        minus[3*i+a]-=h*fixture[3*i+b];
      }
      for(unsigned i=0; i<3; ++i) {
        bp[3*i+a]+=h*normalBox[3*i+b];
        bm[3*i+a]-=h*normalBox[3*i+b];
      }
      const double numerical=(engine.evaluate(plus,bp).value-engine.evaluate(minus,bm).value)/(2*h);
      close(base.virial[3*b+a],numerical,3e-6,"strain derivative");
    }
  }
}

[[maybe_unused]] void runChecks(std::ostream& output) {
  Engine exact(common+"EVALUATION=EXACT");
  const Result exactBase=exact.evaluate(fixture);
  close(exactBase.value,exponentialOracle(fixture),1e-12,"independent exponential oracle");
  output<<"PASS independent-exponential-oracle\n";
  Engine dense(common+compact+"EVALUATION=DENSE");
  Engine sparse(common+compact+"EVALUATION=SPARSE");
  Engine cached(common+compact+"EVALUATION=SPARSE SKIN=0.12 NL_STRIDE=17");
  const Result base=dense.evaluate(fixture);
  compare(base,sparse.evaluate(fixture),"dense sparse");
  compare(base,cached.evaluate(fixture),"dense cached");
  output<<"PASS dense-sparse-cached-parity\n";
  checkDerivatives(exact);
  checkDerivatives(sparse);
  output<<"PASS atomic-anchor-and-strain-derivatives\n";
  for(unsigned frame=0; frame<120; ++frame) {
    auto x=fixture;
    // Move the anchor through and beyond all waters; cross padding/rebuild edges.
    x[27]=-0.5+0.012*frame;
    x[28]=-0.05+0.025*std::sin(0.31*frame);
    // The shared proton exchanges host without changing the atom lists.
    x[15]=0.07+0.0015*(frame%100);
    const Result reference=dense.evaluate(x);
    compare(reference,sparse.evaluate(x),"moving support");
    compare(reference,cached.evaluate(x),"moving cached support");
  }
  // Reuse the high-water allocation through empty and populated evaluations.
  auto empty=fixture;
  empty[27]=empty[28]=empty[29]=0.9;
  for(unsigned repeat=0; repeat<8; ++repeat) {
    const Result zero=cached.evaluate(empty);
    close(zero.value,0,0,"empty reused scratch value");
    for(double force:zero.force) {
      close(force,0,0,"empty reused scratch force");
    }
    const Result restored=cached.evaluate(fixture);
    require(restored.value==base.value && restored.force==base.force &&
            restored.virial==base.virial,"scratch reuse changed the deterministic result");
  }
  output<<"PASS region-entry-exit-host-exchange-and-empty-region\n";
  auto translated=fixture;
  for(unsigned i=0; i<10; ++i) {
    for(unsigned k=0; k<3; ++k) {
      translated[3*i+k]+=0.41;
    }
  }
  compare(base,sparse.evaluate(translated),"translation");
  translated=fixture;
  translated[0]+=2;
  translated[13]-=2;
  translated[27]+=2;
  compare(base,sparse.evaluate(translated),"pbc images");
  output<<"PASS translation-and-pbc-images\n";
  auto permuted=fixture;
  for(unsigned k=0; k<3; ++k) {
    std::swap(permuted[k],permuted[6+k]);
    std::swap(permuted[9+k],permuted[21+k]);
  }
  Result perm=sparse.evaluate(permuted);
  for(unsigned k=0; k<3; ++k) {
    std::swap(perm.force[k],perm.force[6+k]);
    std::swap(perm.force[9+k],perm.force[21+k]);
  }
  compare(base,perm,"permutation");
  output<<"PASS center-and-assigned-permutation\n";
  Engine nonperiodic(common+compact+"NOPBC");
  compare(base,nonperiodic.evaluate(fixture),"nopbc");
  auto changedBox=normalBox;
  changedBox[0]=2.3;
  changedBox[4]=2.2;
  compare(dense.evaluate(fixture,changedBox),cached.evaluate(fixture,changedBox),"box rebuild");
  output<<"PASS nonperiodic-and-box-rebuild\n";
  Engine serial(common+compact+"SERIAL");
  compare(base,serial.evaluate(fixture),"serial mode");
  output<<"PASS serial-mode-parity\n";
  bool caught=false;
  try {
    Engine invalid(common+"PAIR_SWITCH=0.01 PAIR_CUTOFF=0.45");
    invalid.evaluate(fixture);
  } catch(const PLMD::Plumed::Exception& error) {
    caught=std::string(error.what()).find("no CENTER inside PAIR_SWITCH")!=std::string::npos;
  }
  require(caught,"unsafe normalization did not fail");
  output<<"PASS unsafe-normalization-fails-clearly\n";
  caught=false;
  try {
    Engine invalid(common+compact+"NL_STRIDE=10");
  } catch(const PLMD::Plumed::Exception& error) {
    caught=std::string(error.what()).find("SKIN > 0")!=std::string::npos;
  }
  require(caught,"invalid skin did not fail");
  output<<"PASS invalid-buffer-input-fails-clearly\n";
  for(const std::string mode: {
        "EXACT","DENSE","SPARSE"
      }) {
    Engine monitor(common+compact+"EVALUATION="+mode,false);
    for(unsigned frame=0; frame<3; ++frame) {
      const Result result=monitor.evaluate(fixture);
      close(result.value,mode=="EXACT" ? exactBase.value : base.value,1e-11,"monitor value");
      for(double force:result.force) {
        close(force,0,1e-14,"monitor force");
      }
      for(double virial:result.virial) {
        close(virial,0,1e-14,"monitor virial");
      }
    }
  }
  output<<"PASS monitor-only-value-and-zero-force\n";
  const Box triclinic{{2,0,0,0.3,2,0,-0.2,0.25,2}};
  auto wrapped=fixture;
  for(unsigned k=0; k<3; ++k) {
    wrapped[k]+=triclinic[3+k];
    wrapped[27+k]-=triclinic[6+k];
  }
  compare(dense.evaluate(fixture,triclinic),cached.evaluate(wrapped,triclinic),"triclinic images");
  output<<"PASS triclinic-images-and-rebuild\n";
}
}
