/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

// includes
// std
#include <vector>

// RBDyn
#include <rbdyn/config.hh>

#include "Joint.h"

namespace rbd
{
class MultiBody;
struct MultiBodyConfig;

/**
 * Integrate a constant rotation acceleration.
 * @param qi Initial orientation.
 * @param wi Initial rotation speed.
 * @param wD Constant acceleration.
 * @param step Integration step.
 */
RBDYN_DLLAPI std::pair < Eigen::Quaterniond, bool> SO3Integration(
                                                      const Eigen::Quaterniond & qi,
                                                      const Eigen::Vector3d & wi,
                                                      const Eigen::Vector3d & wD,
                                                      double step,
                                                      double relEps = 1e-12,
                                                      double absEps = std::numeric_limits<double>::epsilon(),
                                                      bool breakOnWarning = false);

/**
 * Integrate joint configuration.
 * @param type Joint type.
 * @param alpha Joint velocity vector.
 * @param alphaD Joint acceleration vector.
 * @param step Integration step.
 * @param q Joint configuration vector.
 */
RBDYN_DLLAPI void eulerJointIntegration(Joint::Type type,
                                        const std::vector<double> & alpha,
                                        const std::vector<double> & alphaD,
                                        double step,
                                        std::vector<double> & q);

/**
 * Use the euler method to integrate.
 * @param mb MultiBody used has model.
 * @param mbc Use alphaD, alpha and q. Fill alpha and q.
 * @param step Integration step.
 */
RBDYN_DLLAPI void eulerIntegration(const MultiBody & mb, MultiBodyConfig & mbc, double step);

/// safe version of @see eulerIntegration.
RBDYN_DLLAPI void sEulerIntegration(const MultiBody & mb, MultiBodyConfig & mbc, double step);

} // namespace rbd
