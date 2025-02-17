/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

// include
// Associated header
#include "RBDyn/ZMP.h"

// RBDyn
#include "RBDyn/MultiBodyConfig.h"

namespace rbd
{

Eigen::Vector3d computeCentroidalZMP(MultiBodyConfig & mbc,
                                     Eigen::Vector3d & com,
                                     Eigen::Vector3d & comA,
                                     double altitude)
{
  double g;
  Eigen::Vector3d zmp;
  g = mbc.gravity(2);
  zmp(0) = com(0) - ((com(2) - altitude) * comA(0)) / (comA(2) + g);
  zmp(1) = com(1) - ((com(2) - altitude) * comA(1)) / (comA(2) + g);
  zmp(2) = altitude;
  return zmp;
}

Eigen::Vector3d computeCentroidalZMPNoGravity(Eigen::Vector3d & com, Eigen::Vector3d & comA, double altitude)
{
  Eigen::Vector3d zmp;
  zmp(0) = com(0) - ((com(2) - altitude) * comA(0)) / (comA(2));
  zmp(1) = com(1) - ((com(2) - altitude) * comA(1)) / (comA(2));
  zmp(2) = altitude;
  return zmp;
}

Eigen::Vector3d computeCentroidalZMPComplete(MultiBodyConfig & mbc,
                                             Eigen::Vector3d & com,
                                             Eigen::Vector3d & comA,
                                             double altitude,
                                             sva::ForceVecd wr_external,
                                             double mass)
{
  Eigen::Vector3d zmp;
  double g;
  g = mbc.gravity(2);
  double denom;
  denom = mass * g - wr_external.force()(2);
  zmp(0) = com(0) - ((com(2) - altitude) * comA(0) * (mass * g / denom)) / (comA(2)) + wr_external.couple()(1) / denom
           + com(2) * wr_external.force()(0) / denom;
  zmp(1) = com(1) - ((com(2) - altitude) * comA(1) * (mass * g / denom)) / (comA(2))-wr_external.couple()(0) / denom
           + com(2) * wr_external.force()(1) / denom;
  zmp(2) = altitude;
  return zmp;
}

} // namespace rbd
