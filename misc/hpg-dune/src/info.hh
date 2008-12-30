#ifndef INFO_HH
#define INFO_HH

#include<vector>

struct RegionInfo {
  double conductivity;
  double heat; // heat source/sink
};

struct BoundaryInfo {
  enum Type {DIRICHLET=1, NEUMANN=2, MIXED=3};
  Type type;
  // coefficients
  // DIRICHLET: none yet, the temperature is assumed 0.0
  // NEUMANN: q - heat flux
  // MIXED: T - external temperature,
  //        h - heat transfer coefficient (due to external fluid/gas convection)
  double parameters[2];
};

extern std::vector<BoundaryInfo> boundaryInfo;
extern std::vector<RegionInfo> regionsInfo;

inline RegionInfo& getRegionInfo(int region) {
  static RegionInfo def = { 1.0, 0.0 };
  
  if (region<=0 || region>regionsInfo.size())
  return def;
  else
    return regionsInfo[region-1];
}

inline BoundaryInfo& getBoundaryInfo(int boundaryId) {
  static BoundaryInfo def = {BoundaryInfo::DIRICHLET, 0.0, 0.0};

  if (boundaryId>boundaryInfo.size()) // assume dirichlet if missing
    return def;
  else
    return boundaryInfo[boundaryId-1];
}

#endif
