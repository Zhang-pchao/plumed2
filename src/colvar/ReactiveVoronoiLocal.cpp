/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Colvar.h"
#include "core/ActionRegister.h"
#include "tools/Communicator.h"
#include "tools/LinkCells.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR VORONOI_LOCAL_COORDINATION
/*
Measure smooth coordination-defect activity inside or outside a spatial reaction region.

The fixed CENTERS and ASSIGNED lists contain the complete chemical candidate
universe, not a fixed list of reacting molecules. REGION_CENTER is one physical
or virtual atom. The contribution of a center is one inside CORE_RADIUS, falls
smoothly to zero at OUTER_RADIUS, and is zero outside. The scalar is

\f[ S=\sum_i\chi_i(n_i-\nu_i)^2,\qquad
 n_i=\sum_j w_{ij},\qquad w_{ij}=g(r_{ij})/\sum_k g(r_{kj}). \f]

REFERENCE supplies the neutral occupancies (one number or one per CENTER).
For water oxygen use 2; other chemistries require appropriate references.
DEFECT=POSITIVE retains only positive occupancy defects, DEFECT=NEGATIVE
retains only negative defects, and the default DEFECT=TOTAL retains both.
The outputs are squared activities, not integer ion counts. Their sum equals
the TOTAL value and derivatives. The signed parts are continuously differentiable
at zero defect, but their second derivatives need not be continuous.

REGION_GEOMETRY=SLAB_X, SLAB_Y or SLAB_Z uses the absolute separation
along the named Cartesian axis. CORE_RADIUS and OUTER_RADIUS are then the
inner and outer half-widths; the default SPHERE uses radial distance. The
normal is fixed, while REGION_CENTER may be a COM or CENTER virtual atom.
All anchor and box derivatives are retained. Membership is spatial and does
not constrain the number or identity of molecules in the region.

REGION_SIDE=OUTSIDE uses one minus the usual spatial weight. It includes all
exterior defects, including defects formed independently far from the core.
It therefore requires complete assignment coverage and retains global candidate
work even in SPARSE mode. This is an optional confinement diagnostic, not a
local-cost shortcut. Combining an OUTSIDE activity with UPPER_WALLS defines a
finite conditional-ensemble restraint. Such a wall suppresses exterior activity
but cannot absolutely prohibit it. Keep positive and negative activities separate.
No division by a defect population is used, so neutral states are nonsingular.
In a contracted PIMD PLUMED instance this action sees contracted coordinates:
a wall intended for physical bead defects requires an uncontracted physical
coordinate view, not an implicitly contracted one. Bias removal does not restore
states excluded by a hard wall or recover a different quantum ensemble.

GLOBAL sets every spatial weight to one. Omit REGION_CENTER, CORE_RADIUS and
OUTER_RADIUS when using GLOBAL. The same complete normalization and optional
compact pair kernel apply. Compare low-frequency GLOBAL monitors to local
monitors to distinguish loss of local visibility from disappearance of defects.
For compact GLOBAL evaluation every ASSIGNED atom must have an unswitched
CENTER, including atoms outside the local reaction region. GLOBAL does not
provide spatial acceleration; choose its reporting stride independently.

These geometric defects are not electronic charges. This activity does not
measure ion-pair separation and does not confine reactions to the core.

By default \f$g(r)=\exp(-\kappa r)c(r)\f$. The compact factor c is one up to
PAIR_SWITCH and zero beyond PAIR_CUTOFF. Both c and the spatial factor chi use
\f$1-10t^3+15t^4-6t^5\f$ in their transition intervals. All spatial, pair-kernel,
reference-center and box derivatives are included. The compact kernel is an
explicit approximation to the exponential assignment in VORONOI_COORDINATION;
converge its range against EVALUATION=EXACT before scientific use.

EVALUATION=SPARSE uses native link cells and retains every assigned atom that
can affect the region, together with ALL its competing centers inside the pair
cutoff. Competing centers can be outside the reaction region. EVALUATION=DENSE
uses the identical compact formula with complete pair enumeration, for testing.
EVALUATION=EXACT uses all exponential pairs without truncation. These alternatives
permit separate tests of sparse implementation error and compact-kernel error.
For compact evaluation, every assigned atom contributing to the region must have
at least one center at or inside PAIR_SWITCH. Otherwise calculation stops: a
vanishing denominator must not create an artificial discontinuity at the cutoff.

SKIN and NL_STRIDE control only candidate-list reuse. A list is also rebuilt
when any requested atom (including the reference center) moves more than half
SKIN, the box changes or time is restarted. NL_STRIDE greater than one requires
positive SKIN. Lists are padded both for region entry and pair entry. Storage,
coordinate exchange and list rebuilds can still scale with the total system.
Do not infer simulation acceleration from a small reaction region alone.

Distances use current PLUMED length units. NOPBC uses direct distances and
requires consistently unwrapped coordinates. With PBC, select region and kernel
supports that avoid minimum-image switching surfaces. A moving virtual center
must itself have a continuous, PBC-consistent definition. No quantum bead counts
are changed by this action.

```plumed
UNITS LENGTH=A
site: FIXEDATOM AT=0,0,0
local: VORONOI_LOCAL_COORDINATION ...
  CENTERS=1,4,7 ASSIGNED=2,3,5,6,8,9 REFERENCE=2
  REGION_CENTER=site CORE_RADIUS=2 OUTER_RADIUS=4
  KAPPA=4 PAIR_SWITCH=2.5 PAIR_CUTOFF=4.5
...
PRINT ARG=local FILE=COLVAR
```

For separate whole-system diagnostics, with the same chemical candidate lists:

```plumed
UNITS LENGTH=A
positive: VORONOI_LOCAL_COORDINATION ...
  CENTERS=1,4,7 ASSIGNED=2,3,5,6,8,9 REFERENCE=2
  GLOBAL DEFECT=POSITIVE KAPPA=4 PAIR_SWITCH=2.5 PAIR_CUTOFF=4.5
...
negative: VORONOI_LOCAL_COORDINATION ...
  CENTERS=1,4,7 ASSIGNED=2,3,5,6,8,9 REFERENCE=2
  GLOBAL DEFECT=NEGATIVE KAPPA=4 PAIR_SWITCH=2.5 PAIR_CUTOFF=4.5
...
PRINT ARG=positive,negative STRIDE=100 FILE=DEFECTS
```

The distances above illustrate syntax, not transferable converged parameters.
Start with EVALUATION=EXACT and DENSE comparisons on representative structures.
Biasing is optional; ordinary MD, PIMD and offline driver evaluation can all use
this action. For PIMD preserve the chosen centroid/CV-mean/shared-density contract.
*/
//+ENDPLUMEDOC

class VoronoiLocalCoordination : public Colvar {
  struct Row {
    unsigned assigned;
    std::vector<unsigned> centers;
  };
  struct Pair {
    unsigned center;
    unsigned assigned;
    Vector delta;
    double distance;
    double weight;
    double logarithmicDerivative;
  };
  unsigned centersCount_, assignedCount_, anchorIndex_, listStride_;
  double kappa_, coreRadius_, outerRadius_, pairSwitch_, pairCutoff_, skin_;
  bool pbc_, serial_, haveList_, global_;
  std::string evaluation_, defect_;
  int regionAxis_=-1;
  bool outside_=false;
  std::vector<double> reference_;
  std::vector<Row> rows_;
  // Retain allocation capacity, but reset every value on every evaluation.
  std::vector<double> chiScratch_, dchiScratch_, radiiScratch_, occupancyScratch_, sensitivityScratch_;
  std::vector<Pair> pairScratch_;
  std::vector<std::pair<unsigned,unsigned>> rangeScratch_;
  std::vector<Vector> derivativeScratch_;
  std::vector<Vector> savedPositions_;
  Tensor savedBox_;
  long long savedStep_;
  static double taper(double, double, double, double&);
  Vector displacement(unsigned, unsigned) const;
  double regionDistance(const Vector&) const;
  Vector regionDirection(Vector) const;
  bool needsList();
  void buildList();
public:
  explicit VoronoiLocalCoordination(const ActionOptions&);
  static void registerKeywords(Keywords&);
  void calculate() override;
};

PLUMED_REGISTER_ACTION(VoronoiLocalCoordination,"VORONOI_LOCAL_COORDINATION")

void VoronoiLocalCoordination::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms","CENTERS","Complete list of atoms receiving the assigned atoms");
  keys.add("atoms","ASSIGNED","Complete list of atoms being smoothly assigned");
  keys.add("atoms","REGION_CENTER","One physical or virtual atom anchoring the spatial region");
  keys.add("compulsory","REFERENCE","Neutral occupancy, one number or one per CENTER");
  keys.add("compulsory","KAPPA","Positive inverse-length assignment sharpness");
  keys.add("compulsory","CORE_RADIUS","0","Nonnegative radius with unit spatial weight");
  keys.add("compulsory","OUTER_RADIUS","0","Radius where spatial weight becomes zero");
  keys.add("optional","PAIR_SWITCH","End of the unmodified exponential pair kernel; required for compact evaluation");
  keys.add("optional","PAIR_CUTOFF","End of the smooth pair taper; required for compact evaluation");
  keys.add("compulsory","EVALUATION","SPARSE","SPARSE or DENSE compact kernel, or EXACT exponential reference");
  keys.add("compulsory","SKIN","0","Nonnegative padding for safe sparse candidate reuse");
  keys.add("compulsory","NL_STRIDE","1","Maximum interval between sparse candidate rebuilds");
  keys.addFlag("GLOBAL",false,"Monitor all centers without a spatial region; omit REGION_CENTER and radii");
  keys.add("compulsory","REGION_GEOMETRY","SPHERE","SPHERE, SLAB_X, SLAB_Y or SLAB_Z with a fixed Cartesian normal");
  keys.add("compulsory","REGION_SIDE","INSIDE","INSIDE or OUTSIDE spatial weight; OUTSIDE requires global assignment coverage");
  keys.add("compulsory","DEFECT","TOTAL","TOTAL squared activity, or its POSITIVE or NEGATIVE occupancy-defect part");
  keys.addFlag("SERIAL",false,"Evaluate all assignments on each rank without internal MPI partitioning");
  keys.setValueDescription("scalar","the spatially weighted squared coordination defects");
}

VoronoiLocalCoordination::VoronoiLocalCoordination(const ActionOptions& ao):
  PLUMED_COLVAR_INIT(ao), centersCount_(0), assignedCount_(0), anchorIndex_(0),
  listStride_(1), kappa_(0), coreRadius_(0), outerRadius_(0), pairSwitch_(0),
  pairCutoff_(0), skin_(0), pbc_(true), serial_(false), haveList_(false), global_(false), savedStep_(0) {
  std::vector<AtomNumber> centers, assigned, anchor;
  parseAtomList("CENTERS",centers);
  parseAtomList("ASSIGNED",assigned);
  parseAtomList("REGION_CENTER",anchor);
  parseFlag("GLOBAL",global_);
  parse("DEFECT",defect_);
  std::string geometry,side;
  parse("REGION_GEOMETRY",geometry);
  parse("REGION_SIDE",side);
  if(geometry=="SLAB_X") {
    regionAxis_=0;
  } else if(geometry=="SLAB_Y") {
    regionAxis_=1;
  } else if(geometry=="SLAB_Z") {
    regionAxis_=2;
  } else if(geometry!="SPHERE") {
    error("REGION_GEOMETRY must be SPHERE, SLAB_X, SLAB_Y or SLAB_Z");
  }
  if(side=="OUTSIDE") {
    outside_=true;
  } else if(side!="INSIDE") {
    error("REGION_SIDE must be INSIDE or OUTSIDE");
  }
  if(global_ && (outside_ || regionAxis_>=0)) {
    error("GLOBAL does not accept a regional geometry or OUTSIDE weight");
  }
  if(defect_!="TOTAL" && defect_!="POSITIVE" && defect_!="NEGATIVE") {
    error("DEFECT must be TOTAL, POSITIVE or NEGATIVE");
  }
  if(centers.empty() || assigned.empty() || (!global_ && anchor.size()!=1) || (global_ && !anchor.empty())) {
    error("CENTERS and ASSIGNED must be nonempty; provide one REGION_CENTER for a spatial region, or omit it with GLOBAL");
  }
  std::set<unsigned> seen;
  for(const auto& atom:centers) {
    if(!seen.insert(atom.serial()).second) {
      error("CENTERS must not contain duplicate atoms");
    }
  }
  for(const auto& atom:assigned) {
    if(!seen.insert(atom.serial()).second) {
      error("ASSIGNED must be unique and disjoint from CENTERS");
    }
  }
  centersCount_=centers.size();
  assignedCount_=assigned.size();
  anchorIndex_=centersCount_+assignedCount_;
  parseVector("REFERENCE",reference_);
  if(reference_.size()==1) {
    reference_.resize(centersCount_,reference_[0]);
  }
  if(reference_.size()!=centersCount_) {
    error("REFERENCE needs one number or one number per CENTER");
  }
  for(double x:reference_) {
    if(!std::isfinite(x)) {
      error("REFERENCE must be finite");
    }
  }
  parse("KAPPA",kappa_);
  parse("CORE_RADIUS",coreRadius_);
  parse("OUTER_RADIUS",outerRadius_);
  parse("PAIR_SWITCH",pairSwitch_);
  parse("PAIR_CUTOFF",pairCutoff_);
  parse("EVALUATION",evaluation_);
  parse("SKIN",skin_);
  int stride=1;
  parse("NL_STRIDE",stride);
  if(!std::isfinite(kappa_) || kappa_<=0) {
    error("KAPPA must be finite and positive");
  }
  if(global_ && (coreRadius_!=0 || outerRadius_!=0)) {
    error("GLOBAL must not be combined with CORE_RADIUS or OUTER_RADIUS");
  }
  if(!global_ && (!std::isfinite(coreRadius_) || !std::isfinite(outerRadius_) || coreRadius_<0 || outerRadius_<=coreRadius_)) {
    error("Require finite 0 <= CORE_RADIUS < OUTER_RADIUS");
  }
  if(evaluation_!="EXACT" && evaluation_!="DENSE" && evaluation_!="SPARSE") {
    error("EVALUATION must be SPARSE, DENSE or EXACT");
  }
  if(evaluation_!="EXACT" && (!std::isfinite(pairSwitch_) || !std::isfinite(pairCutoff_) || pairSwitch_<=0 || pairCutoff_<=pairSwitch_)) {
    error("Compact evaluation requires finite 0 < PAIR_SWITCH < PAIR_CUTOFF");
  }
  if(!std::isfinite(skin_) || skin_<0 || stride<1 || (stride>1 && skin_==0)) {
    error("SKIN must be finite and nonnegative; NL_STRIDE must be positive and needs SKIN > 0 when greater than one");
  }
  if(!std::isfinite(outerRadius_+skin_) || !std::isfinite(pairCutoff_+skin_)) {
    error("Radius plus SKIN must be finite");
  }
  listStride_=stride;
  bool nopbc=false;
  parseFlag("NOPBC",nopbc);
  pbc_=!nopbc;
  parseFlag("SERIAL",serial_);
  checkRead();
  centers.insert(centers.end(),assigned.begin(),assigned.end());
  if(!global_) {
    centers.push_back(anchor[0]);
  }
  addValueWithDerivatives();
  setNotPeriodic();
  requestAtoms(centers);
  log.printf("  local reactive Voronoi: %u centers, %u assigned, %s evaluation\n",centersCount_,assignedCount_,evaluation_.c_str());
  log.printf("  region core %g outer %g; pair switch %g cutoff %g; skin %g stride %u\n",coreRadius_,outerRadius_,pairSwitch_,pairCutoff_,skin_,listStride_);
  log.printf("  defect %s; global region %s\n",defect_.c_str(),global_ ? "yes" : "no");
}

// Use the factored endpoint polynomial to avoid cancellation near zero support.
double VoronoiLocalCoordination::taper(double r,double inner,double outer,double& derivative) {
  derivative=0;
  if(r<=inner) {
    return 1;
  }
  if(r>=outer) {
    return 0;
  }
  const double u=(outer-r)/(outer-inner);
  derivative=-30*u*u*(1-u)*(1-u)/(outer-inner);
  return u*u*u*(10+u*(-15+6*u));
}

Vector VoronoiLocalCoordination::displacement(unsigned a,unsigned b) const {
  return pbc_ ? pbcDistance(getPosition(a),getPosition(b)) : delta(getPosition(a),getPosition(b));
}

double VoronoiLocalCoordination::regionDistance(const Vector& d) const {
  return regionAxis_<0 ? d.modulo() : std::fabs(d[regionAxis_]);
}

Vector VoronoiLocalCoordination::regionDirection(Vector d) const {
  if(regionAxis_>=0) {
    for(unsigned k=0; k<3; ++k)
      if(static_cast<int>(k)!=regionAxis_) {
        d[k]=0;
      }
  }
  return d;
}

bool VoronoiLocalCoordination::needsList() {
  if(!haveList_ || skin_==0 || getStep()<=savedStep_ || getStep()-savedStep_>=listStride_) {
    return true;
  }
  if(pbc_) {
    const Tensor box=getPbc().getBox();
    for(unsigned a=0; a<3; ++a) {
      for(unsigned b=0; b<3; ++b) {
        if(box(a,b)!=savedBox_(a,b)) {
          return true;
        }
      }
    }
  }
  for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
    const Vector d=pbc_ ? pbcDistance(savedPositions_[i],getPosition(i)) : delta(savedPositions_[i],getPosition(i));
    if(d.modulo2()>skin_*skin_*0.25) {
      return true;
    }
  }
  return false;
}

void VoronoiLocalCoordination::buildList() {
  rows_.clear();
  if(evaluation_!="SPARSE") {
    for(unsigned j=0; j<assignedCount_; ++j) {
      Row row;
      row.assigned=j;
      row.centers.resize(centersCount_);
      std::iota(row.centers.begin(),row.centers.end(),0);
      rows_.push_back(std::move(row));
    }
  } else {
    LinkCells cells(comm);
    cells.setCutoff(pairCutoff_+skin_);
    const auto& positions=getPositions();
    const View<const Vector> physical{positions.data(),anchorIndex_};
    if(pbc_) {
      cells.setupCells(physical,getPbc());
    } else {
      cells.setupCells(physical);
    }
    std::vector<unsigned> centerIndices(centersCount_), hydrogenIndices(assignedCount_);
    std::iota(centerIndices.begin(),centerIndices.end(),0);
    std::iota(hydrogenIndices.begin(),hydrogenIndices.end(),centersCount_);
    const auto centerCells=cells.getCollection(View<const Vector> {positions.data(),centersCount_},make_const_view(centerIndices));
    const auto assignedCells=cells.getCollection(View<const Vector> {positions.data()+centersCount_,assignedCount_},make_const_view(hydrogenIndices));
    std::vector<unsigned> neighbors(27);
    std::vector<bool> retained(assignedCount_,global_ || outside_);
    const double cutoff2=(pairCutoff_+skin_)*(pairCutoff_+skin_);
    for(unsigned i=0; !global_ && !outside_ && i<centersCount_; ++i) {
      if(regionDistance(displacement(anchorIndex_,i))>outerRadius_+skin_) {
        continue;
      }
      unsigned count=0;
      cells.addRequiredCells(cells.findMyCell(getPosition(i)),count,neighbors,pbc_);
      for(unsigned c=0; c<count; ++c) {
        for(unsigned j:assignedCells.getCellIndexes(neighbors[c])) {
          if(displacement(i,j).modulo2()<=cutoff2) {
            retained[j-centersCount_]=true;
          }
        }
      }
    }
    for(unsigned j=0; j<assignedCount_; ++j) {
      if(!retained[j]) {
        continue;
      }
      Row row;
      row.assigned=j;
      unsigned count=0;
      cells.addRequiredCells(cells.findMyCell(getPosition(j+centersCount_)),count,neighbors,pbc_);
      for(unsigned c=0; c<count; ++c) {
        for(unsigned i:centerCells.getCellIndexes(neighbors[c])) {
          if(displacement(i,j+centersCount_).modulo2()<=cutoff2) {
            row.centers.push_back(i);
          }
        }
      }
      std::sort(row.centers.begin(),row.centers.end());
      row.centers.erase(std::unique(row.centers.begin(),row.centers.end()),row.centers.end());
      rows_.push_back(std::move(row));
    }
  }
  savedPositions_=getPositions();
  savedBox_=getPbc().getBox();
  savedStep_=getStep();
  haveList_=true;
}

void VoronoiLocalCoordination::calculate() {
  for(const Vector& position:getPositions()) {
    for(unsigned k=0; k<3; ++k) {
      if(!std::isfinite(position[k])) {
        error("All candidate and REGION_CENTER coordinates must be finite");
      }
    }
  }
  if(pbc_) {
    const Tensor box=getPbc().getBox();
    bool valid=std::isfinite(box.determinant()) && std::fabs(box.determinant())>1e-14;
    for(unsigned a=0; a<3; ++a) {
      for(unsigned b=0; b<3; ++b) {
        valid=valid && std::isfinite(box(a,b));
      }
    }
    if(!valid) {
      error("PBC requires a finite nondegenerate box; use NOPBC for an unwrapped nonperiodic system");
    }
  }
  if(!haveList_ || (evaluation_=="SPARSE" && needsList())) {
    buildList();
  }
  chiScratch_.assign(centersCount_,0);
  dchiScratch_.assign(centersCount_,0);
  radiiScratch_.assign(centersCount_,0);
  occupancyScratch_.assign(centersCount_,0);
  auto& chi=chiScratch_;
  auto& dchi=dchiScratch_;
  auto& radii=radiiScratch_;
  auto& occupancy=occupancyScratch_;
  for(unsigned i=0; i<centersCount_; ++i) {
    if(global_) {
      chi[i]=1;
    } else {
      radii[i]=regionDistance(displacement(anchorIndex_,i));
      chi[i]=taper(radii[i],coreRadius_,outerRadius_,dchi[i]);
      if(outside_) {
        chi[i]=1-chi[i];
        dchi[i]=-dchi[i];
      }
    }
  }
  pairScratch_.clear();
  rangeScratch_.clear();
  auto& pairs=pairScratch_;
  auto& ranges=rangeScratch_;
  const unsigned rank=serial_ ? 0 : comm.Get_rank();
  const unsigned ranks=serial_ ? 1 : comm.Get_size();
  int invalid=0;
  for(const Row& row:rows_) {
    if(row.assigned%ranks!=rank) {
      continue;
    }
    const unsigned begin=pairs.size();
    double maximum=-std::numeric_limits<double>::infinity();
    bool relevant=global_ || outside_, unswitched=false;
    for(unsigned i:row.centers) {
      const Vector d=displacement(i,row.assigned+centersCount_);
      const double r=d.modulo();
      if(evaluation_!="EXACT" && r>=pairCutoff_) {
        continue;
      }
      if(!std::isfinite(r) || r<=std::numeric_limits<double>::epsilon()) {
        invalid=1;
        continue;
      }
      double derivative=0;
      const double c=evaluation_=="EXACT" ? 1 : taper(r,pairSwitch_,pairCutoff_,derivative);
      const double score=-kappa_*r+std::log(c);
      const double dlog=-kappa_+derivative/c;
      if(!std::isfinite(score) || !std::isfinite(dlog)) {
        invalid=1;
        continue;
      }
      relevant=relevant || chi[i]>0;
      unswitched=unswitched || r<=pairSwitch_;
      maximum=std::max(maximum,score);
      pairs.push_back({i,row.assigned,d,r,score,dlog});
    }
    if(!relevant) {
      pairs.resize(begin);
      continue;
    }
    if(evaluation_!="EXACT" && !unswitched) {
      invalid=2;
    }
    double normalization=0;
    for(unsigned p=begin; p<pairs.size(); ++p) {
      pairs[p].weight=std::exp(pairs[p].weight-maximum);
      normalization+=pairs[p].weight;
    }
    if(!std::isfinite(normalization) || normalization<=0) {
      invalid=1;
      pairs.resize(begin);
      continue;
    }
    for(unsigned p=begin; p<pairs.size(); ++p) {
      pairs[p].weight/=normalization;
      occupancy[pairs[p].center]+=pairs[p].weight;
    }
    ranges.emplace_back(begin,pairs.size());
  }
  if(!serial_) {
    comm.Max(invalid);
  }
  if(invalid) {
    error(invalid==2 ? "A retained ASSIGNED atom has no CENTER inside PAIR_SWITCH; enlarge PAIR_SWITCH and PAIR_CUTOFF and validate the kernel" : "Invalid assignment distance or normalization; check overlapping atoms and KAPPA/range parameters");
  }
  if(!serial_) {
    comm.Sum(occupancy);
  }
  sensitivityScratch_.assign(centersCount_,0);
  derivativeScratch_.assign(getNumberOfAtoms(),Vector());
  auto& sensitivity=sensitivityScratch_;
  auto& derivatives=derivativeScratch_;
  Tensor virial;
  double value=0;
  for(unsigned i=0; i<centersCount_; ++i) {
    const double raw=occupancy[i]-reference_[i];
    const double q=defect_=="POSITIVE" ? std::max(raw,0.0) :
                   (defect_=="NEGATIVE" ? std::min(raw,0.0) : raw);
    value+=chi[i]*q*q;
    sensitivity[i]=2*chi[i]*q;
    if((serial_ || rank==0) && radii[i]>0 && dchi[i]!=0) {
      const Vector gradient=(q*q*dchi[i]/radii[i])*regionDirection(displacement(anchorIndex_,i));
      derivatives[i]+=gradient;
      derivatives[anchorIndex_]-=gradient;
      // PLUMED box derivatives use displacement as the first tensor index.
      virial-=Tensor(displacement(anchorIndex_,i),gradient);
    }
  }
  for(const auto& range:ranges) {
    double mean=0;
    for(unsigned p=range.first; p<range.second; ++p) {
      mean+=sensitivity[pairs[p].center]*pairs[p].weight;
    }
    for(unsigned p=range.first; p<range.second; ++p) {
      const Pair& pair=pairs[p];
      const Vector gradient=(pair.weight*(sensitivity[pair.center]-mean)*pair.logarithmicDerivative/pair.distance)*pair.delta;
      derivatives[pair.center]-=gradient;
      derivatives[pair.assigned+centersCount_]+=gradient;
      virial-=Tensor(gradient,pair.delta);
    }
  }
  if(!serial_) {
    comm.Sum(&derivatives[0][0],3*derivatives.size());
    comm.Sum(virial);
  }
  if(!std::isfinite(value)) {
    error("Local reactive Voronoi value is non-finite");
  }
  for(unsigned i=0; i<derivatives.size(); ++i) {
    for(unsigned k=0; k<3; ++k) {
      if(!std::isfinite(derivatives[i][k])) {
        error("Local reactive Voronoi derivative is non-finite");
      }
    }
    setAtomsDerivatives(i,derivatives[i]);
  }
  setBoxDerivatives(virial);
  setValue(value);
}

}
}
