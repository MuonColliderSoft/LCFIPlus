// TrackSelector.h

#ifndef TrackSelector_h
#define TrackSelector_h 1

#include "lcfiplus.h"
#include <vector>
#include <cmath>

using namespace std;
using namespace lcfiplus;

namespace lcfiplus {

class TrackSelectorConfig {
 public:
  // cuts which are combined using the AND scheme
  double minD0;
  double maxD0;
  double minD0Err;
  double maxD0Err;
  double minD0Sig;
  double maxD0Sig;
  double minZ0;
  double maxZ0;
  double minZ0Err;
  double maxZ0Err;
  double minZ0Sig;
  double maxZ0Sig;
  double minD0Z0Sig;
  double maxD0Z0Sig;
  double minPt;
  double maxInnermostHitRadius;
  // cuts which are combined using the OR scheme, then AND'd with the AND schemes above
  int minTpcHits;
  double minTpcHitsMinPt;
  int minFtdHits;
  int minVtxHits;
  int minVtxPlusFtdHits;

  double minR0;
  double maxR0;

  TrackSelectorConfig() {
    minD0 = 0.;
    maxD0 = 1e+300;
    minD0Err = 0.;
    maxD0Err = 1e+300;
    minD0Sig = 0.;
    maxD0Sig = 1e+300;
    minZ0 = 0.;
    maxZ0 = 1e+300;
    minZ0Err = 0.;
    maxZ0Err = 1e+300;
    minZ0Sig = 0.;
    maxZ0Sig = 1e+300;
    minD0Z0Sig = 0.;
    maxD0Z0Sig = 1e+300;
    minPt = 0.;
    maxInnermostHitRadius = 1e+300;

    //added
    minR0 = 0.;
    maxR0 = 1e+300;

    //put to zero bc tracker geometry is different
    minTpcHits = 0;
    minTpcHitsMinPt = 0;
    minFtdHits = 0;
    minVtxHits = 0;
    minVtxPlusFtdHits = 0;
  }
};

class TrackSelector {
 public:
  vector<const Track*> operator () (const vector<const Track*>& tracks, TrackSelectorConfig& config, const Vertex* ip = 0) {
    vector<const Track*> ret;

    for (unsigned int i=0; i<tracks.size(); i++) {
      if (passesCut(tracks[i], config, ip))
        ret.push_back(tracks[i]);
    }

    return ret;
  }

  bool passesCut(const Track* trk, const TrackSelectorConfig& cfg, const Vertex* ip = 0) {

    // AND cuts
    if (fabs(trk->getD0()) < cfg.minD0) {
      if(verboseDebug) std::cout << "MinD0: " << fabs(trk->getD0()) << " / " << cfg.minD0 << std::endl;
      return false;
    }
    if (fabs(trk->getD0()) > cfg.maxD0) {
      if(verboseDebug) std::cout << "MaxD0: " << fabs(trk->getD0()) << " / " << cfg.maxD0 << std::endl;
      return false;
    }
    if (sqrt(trk->getCovMatrix()[tpar::d0d0]) < cfg.minD0Err) {
      if(verboseDebug) std::cout << "MinD0Err: " << sqrt(trk->getCovMatrix()[tpar::d0d0])  << " / " << cfg.minD0Err << std::endl;
      return false;
    }
    if (sqrt(trk->getCovMatrix()[tpar::d0d0]) > cfg.maxD0Err)  {
      if(verboseDebug) std::cout << "MaxD0Err: " << sqrt(trk->getCovMatrix()[tpar::d0d0])  << " / " << cfg.maxD0Err << std::endl;
      return false;
    }
    double d0sig = fabs(trk->getD0()) / sqrt(trk->getCovMatrix()[tpar::d0d0]);
    if ( d0sig < cfg.minD0Sig) {
      if(verboseDebug) std::cout << "MinD0Sig: " << d0sig << " / " << cfg.minD0Sig << std::endl;
      return false;
    }
    if ( d0sig > cfg.maxD0Sig) {
      if(verboseDebug) std::cout << "MaxD0Sig: " << d0sig << " / " << cfg.maxD0Sig << std::endl;
      return false;
    }

    double z0 = (ip ? fabs(trk->getZ0() - ip->getZ()) : fabs(trk->getZ0()) );
    if (z0 < cfg.minZ0) {
      if(verboseDebug) std::cout << "MinZ0: " << z0 << " / " << cfg.minZ0 << std::endl;
      return false;
    }
    if (z0 > cfg.maxZ0) {
      if(verboseDebug) std::cout << "MaxZ0: " << z0 << " / " << cfg.maxZ0 << std::endl;
      return false;
    }
    if (sqrt(trk->getCovMatrix()[tpar::z0z0]) < cfg.minZ0Err) {
      if(verboseDebug) std::cout << "MinZ0Err: " << sqrt(trk->getCovMatrix()[tpar::z0z0]) << " / " << cfg.minZ0Err << std::endl;
      return false;
    }
    if (sqrt(trk->getCovMatrix()[tpar::z0z0]) > cfg.maxZ0Err) {
      if(verboseDebug) std::cout << "MaxZ0Err: " << sqrt(trk->getCovMatrix()[tpar::z0z0]) << " / " << cfg.maxZ0Err << std::endl;
      return false;
    }
    double z0sig = z0 / sqrt(trk->getCovMatrix()[tpar::z0z0]);
    if ( z0sig < cfg.minZ0Sig) {
      if(verboseDebug) std::cout << "MinZ0Sig: " << z0sig << " / " << cfg.minZ0Sig << std::endl;
      return false;
    }
    if ( z0sig > cfg.maxZ0Sig) {
      if(verboseDebug) std::cout << "MaxZ0Sig: " << z0sig << " / " << cfg.maxZ0Sig << std::endl;
      return false;
    }

    if (sqrt(d0sig * d0sig + z0sig * z0sig) < cfg.minD0Z0Sig) {
      if(verboseDebug) std::cout << "MinD0Z0Sig: " << sqrt(d0sig * d0sig + z0sig * z0sig) << " / " << cfg.minD0Z0Sig << std::endl;
      return false;
    }
    if (sqrt(d0sig * d0sig + z0sig * z0sig) > cfg.maxD0Z0Sig) {
      if(verboseDebug) std::cout << "MaxD0Z0Sig: " << sqrt(d0sig * d0sig + z0sig * z0sig) << " / " << cfg.maxD0Z0Sig << std::endl;
      return false;
    }

    if (trk->Pt() < cfg.minPt) {
      if(verboseDebug) std::cout << "MinPt: " << trk->Pt() << " / " << cfg.minPt << std::endl;
      return false;
    }
    if (trk->getRadiusOfInnermostHit() > cfg.maxInnermostHitRadius) {
      if(verboseDebug) std::cout << "mMxInnermostHitRadius: " << trk->getRadiusOfInnermostHit()  << " / " << cfg.maxInnermostHitRadius << std::endl;
      return false;
    }

    //added
    if( sqrt( pow(trk->getD0(), 2) + pow(trk->getZ0(), 2) ) < cfg.minR0 ) {
      if(verboseDebug) std::cout << "MinR0: " << sqrt( pow(trk->getD0(), 2) + pow(trk->getZ0(), 2) ) << " / " << cfg.minR0 << std::endl;
      return false;
    }
    if( sqrt( pow(trk->getD0(), 2) + pow(trk->getZ0(), 2) ) > cfg.maxR0 ) {
      if(verboseDebug) std::cout << "MaxR0: " << sqrt( pow(trk->getD0(), 2) + pow(trk->getZ0(), 2) ) << " / " << cfg.maxR0 << std::endl;
      return false;
    }

    // OR cuts
    if (trk->getFtdHits() >= cfg.minFtdHits) return true;
    if (trk->getVtxHits() >= cfg.minVtxHits) return true;
    if (trk->getVtxHits() + trk->getFtdHits() >= cfg.minVtxPlusFtdHits) return true;
    if (trk->getTpcHits() >= cfg.minTpcHits && trk->Pt() > cfg.minTpcHitsMinPt) return true;

    return false;
  }

  //c-tor / d-tor
  TrackSelector() {}
  ~TrackSelector() {}
};
}

#endif //TrackSelector_h
