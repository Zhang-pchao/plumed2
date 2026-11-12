#include "local_checks.h"
namespace {
void runRegionChecks(std::ostream& output) {
  const std::string global="v: VORONOI_LOCAL_COORDINATION CENTERS=1-3 ASSIGNED=4-9 REFERENCE=2 KAPPA=12 GLOBAL ";
  std::array<double,3> occupancy{};
  for(unsigned j=3; j<9; ++j) {
    std::array<double,3> raw{};
    double sum=0;
    for(unsigned i=0; i<3; ++i) {
      double r2=0;
      for(unsigned axis=0; axis<3; ++axis) {
        const double d=fixture[3*i+axis]-fixture[3*j+axis];
        r2+=d*d;
      }
      raw[i]=std::exp(-12*std::sqrt(r2));
      sum+=raw[i];
    }
    for(unsigned i=0; i<3; ++i) {
      occupancy[i]+=raw[i]/sum;
    }
  }
  double positiveOracle=0, negativeOracle=0;
  for(double n:occupancy) {
    if(n>2) {
      positiveOracle+=(n-2)*(n-2);
    } else {
      negativeOracle+=(n-2)*(n-2);
    }
  }
  Engine exactPositive(global+"EVALUATION=EXACT DEFECT=POSITIVE");
  Engine exactNegative(global+"EVALUATION=EXACT DEFECT=NEGATIVE");
  close(exactPositive.evaluate(fixture).value,positiveOracle,1e-12,"positive oracle");
  close(exactNegative.evaluate(fixture).value,negativeOracle,1e-12,"negative oracle");
  output<<"PASS independent-signed-global-oracle\n";
  for(const std::string mode: {
        "EXACT","DENSE","SPARSE"
      }) {
    for(const std::string& scope: {
          common,global
        }) {
      const std::string input=scope+compact+"EVALUATION="+mode;
      Engine total(input),positive(input+" DEFECT=POSITIVE"),negative(input+" DEFECT=NEGATIVE");
      const Result all=total.evaluate(fixture),pos=positive.evaluate(fixture),neg=negative.evaluate(fixture);
      Result sum=pos;
      sum.value+=neg.value;
      for(unsigned i=0; i<sum.force.size(); ++i) {
        sum.force[i]+=neg.force[i];
      }
      for(unsigned i=0; i<9; ++i) {
        sum.virial[i]+=neg.virial[i];
      }
      compare(all,sum,"signed decomposition");
      require(pos.value>0 && neg.value>0,"fixture must exercise both signs");
      checkDerivatives(positive);
      checkDerivatives(negative);
      for(const std::string sign: {
            "POSITIVE","NEGATIVE"
          }) {
        Engine monitor(input+" DEFECT="+sign,false);
        const Result m=monitor.evaluate(fixture);
        close(m.value,sign=="POSITIVE" ? pos.value : neg.value,1e-11,"signed monitor value");
        for(double f:m.force) {
          close(f,0,1e-14,"signed monitor force");
        }
      }
    }
  }
  output<<"PASS signed-decomposition-and-complete-derivatives\n";
  output<<"PASS signed-monitor-only-zero-force\n";
  for(const std::string sign: {
        "TOTAL","POSITIVE","NEGATIVE"
      }) {
    Engine dense(global+compact+"EVALUATION=DENSE DEFECT="+sign);
    Engine sparse(global+compact+"EVALUATION=SPARSE SKIN=0.12 NL_STRIDE=17 DEFECT="+sign);
    for(unsigned frame=0; frame<40; ++frame) {
      auto x=fixture;
      x[15]=0.07+0.004*(frame%30);
      compare(dense.evaluate(x),sparse.evaluate(x),"global sparse movement");
    }
  }
  output<<"PASS global-dense-sparse-cached-parity\n";
  Engine full(global+compact),local(common+compact);
  const Result before=full.evaluate(fixture);
  auto escaped=fixture;
  for(unsigned i=0; i<9; ++i) {
    escaped[3*i+1]+=0.8;
  }
  const Result after=full.evaluate(escaped);
  compare(before,after,"defect translation keeps global activity");
  require(local.evaluate(fixture).value>1e-4,"local initial activity");
  close(local.evaluate(escaped).value,0,1e-13,"escaped activity invisible locally");
  for(unsigned k=27; k<30; ++k) {
    close(after.force[k],0,1e-14,"unused anchor force");
  }
  output<<"PASS escape-is-not-global-recombination\n";
  for(const std::string& bad: {
        common+compact+"GLOBAL",global+compact+"CORE_RADIUS=0.1",global+compact+"DEFECT=NET"
      }) {
    bool caught=false;
    try {
      Engine invalid(bad);
    } catch(const PLMD::Plumed::Exception&) {
      caught=true;
    }
    require(caught,"invalid diagnostic configuration did not fail");
  }
  output<<"PASS invalid-global-and-sign-input\n";
  bool caught=false;
  try {
    Engine isolated(global+compact);
    auto x=fixture;
    x[9]=0.9;
    x[10]=0.9;
    x[11]=0.9;
    isolated.evaluate(x);
  } catch(const PLMD::Plumed::Exception&) {
    caught=true;
  }
  require(caught,"global unassigned atom silently ignored");
  output<<"PASS global-missing-normalization-fails\n";
  for(const std::string geometry: {
        "SPHERE","SLAB_X","SLAB_Y","SLAB_Z"
      }) {
    const std::string shape=" REGION_GEOMETRY="+geometry;
    Engine inside(common+compact+shape),outside(common+compact+shape+" REGION_SIDE=OUTSIDE");
    const Result a=inside.evaluate(fixture),b=outside.evaluate(fixture);
    Result sum=a;
    sum.value+=b.value;
    for(unsigned i=0; i<sum.force.size(); ++i) {
      sum.force[i]+=b.force[i];
    }
    for(unsigned i=0; i<9; ++i) {
      sum.virial[i]+=b.virial[i];
    }
    compare(full.evaluate(fixture),sum,"inside plus outside equals global");
    checkDerivatives(inside);
    checkDerivatives(outside);
    Engine dense(common+compact+shape+" EVALUATION=DENSE REGION_SIDE=OUTSIDE");
    Engine cached(common+compact+shape+" EVALUATION=SPARSE REGION_SIDE=OUTSIDE SKIN=0.12 NL_STRIDE=17");
    for(unsigned frame=0; frame<30; ++frame) {
      auto x=fixture;
      x[27]+=0.013*frame;
      compare(dense.evaluate(x),cached.evaluate(x),"outside moving anchor sparse parity");
    }
  }
  output<<"PASS sphere-slab-inside-outside-derivatives-and-parity\n";
  for(const std::string bad: {
        " REGION_GEOMETRY=CYLINDER"," REGION_SIDE=BOTH"
      }) {
    bool rejected=false;
    try {
      Engine invalid(common+compact+bad);
    } catch(const PLMD::Plumed::Exception&) {
      rejected=true;
    }
    require(rejected,"invalid general region input did not fail");
  }
  output<<"PASS invalid-general-region-input\n";
  const std::string wall="b: UPPER_WALLS ARG=v AT=0 KAPPA=2 EXP=2 EPS=1";
  for(const std::string geometry: {
        "SPHERE","SLAB_X"
      }) {
    const std::string input=common+compact+" REGION_SIDE=OUTSIDE REGION_GEOMETRY="+geometry;
    Engine activity(input),confined(input,true,wall);
    const double value=activity.evaluate(fixture).value;
    close(confined.evaluate(fixture).value,2*value*value,1e-12,"outside wall scalar oracle");
    checkDerivatives(confined);
    std::string mobileInput=input;
    mobileInput.replace(mobileInput.find("REGION_CENTER=10"),16,"REGION_CENTER=site");
    Engine mobile(mobileInput,true,wall,"site: COM ATOMS=1,3");
    checkDerivatives(mobile);
  }
  std::ostringstream neutralReference;
  neutralReference<<std::setprecision(17)<<"REFERENCE="<<occupancy[0]<<","<<occupancy[1]<<","<<occupancy[2];
  std::string neutral=common+" EVALUATION=EXACT REGION_SIDE=OUTSIDE";
  neutral.replace(neutral.find("REFERENCE=2"),11,neutralReference.str());
  Engine neutralWall(neutral,true,wall);
  const auto zero=neutralWall.evaluate(fixture);
  close(zero.value,0,1e-14,"neutral wall energy");
  for(double f:zero.force) {
    close(f,0,1e-14,"neutral wall force");
  }
  output<<"PASS optional-wall-mobile-com-and-neutral-derivatives\n";


}
}
