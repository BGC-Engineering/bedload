
#' Converting dimensionless flux to flux
#'
#' this function converts dimensionless bedload transport flux values to fluxes
#' in cubic meters per second per meter
#' @param q_star dimensionless bedload flux
#' @param d median bed sediment grain diameter (m)
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param Gs the submerged specific gravity of sediment. Defaults to 1.65
#' @export
qb_conversion <- function(q_star, d, g = 9.81, Gs = 1.65){
  qb <- q_star * sqrt(Gs * g * d^3)
  return(qb)
}



#' Calculating the cross section average shear stress
#'
#' this function calculates the average shear stress acting on the channel
#' boundary for a given cross section.
#' @param R hydraulic radius for the flow (m)
#' @param S gradient of the water surface slope (m/m)
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @export
tau_1D <- function(R, S, g = 9.81, rho = 1000){
  tau <- g * rho * R * S
  return(tau)
}



#' Estimating the mean flow velocity
#'
#' This function uses the Ferguson (2007) flow resistance law to estimate
#' the mean velocity of the flow using the bed surface D84
#' @param R hydraulic radius for the flow (m)
#' @param S gradient of the water surface slope (m/m)
#' @param d84 84th percentile diameter of bed sediment (m)
#' @param a1 constant fit by Ferguson (2007). Defaults to 6.5
#' @param a2 constant fit by Ferguson (2007). Defaults to 2.5
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @export
u_ferg <- function(R, S, d84, a1 = 6.5, a2 = 2.5, g = 9.81){
  Res = a1 * a2 * (R / d84) / (a1^2 + a2^2 * (R / d84)^(5/3))^(1/2)
  U = Res * sqrt(g * R * S)
  return(U)
}



#' Estimating the dimensionless shear stress
#'
#' This function calculates the dimensionless shear stress using the bed surface
#' median grain size.
#' @param tau the cross section-averaged shear stress (Pa/m2).
#' @param d median bed sediment grain diameter (m)
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @param rho_s sediment density. Defaults to 2650 kg/m3
#' @export
t_star <- function(tau, d, g = 9.81, rho = 1000, rho_s = 2650){
  Gs <- (rho_s - rho) / rho
  t_star <- tau / (g * rho * Gs  * d)
  return(t_star)
}



#' Estimating critical dimensionless shear stress
#'
#' Using the equation published by Soulsby (1997), this function estimates the
#' critical dimensionless shear stress as a function of grain size.
#' @param d median bed sediment grain diameter (m)
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @param rho_s sediment density. Defaults to 2650 kg/m3
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param kv kinematic viscosity of water. Defaults to 1e-6
#' @export
t_crit_soul <- function(d, rho = 1000, rho_s = 2650, g = 9.81, kv = 1e-6) {
  # Estimates the critical Shields parameter for initiation of sediment motion
  # using the hiding-function relationship of Soulsby & Whitehouse (1997).
  #
  # Parameters
  #
  # d     : numeric, grain diameter (m)
  # rho   : numeric, fluid density (kg/m^3). Default 1000.
  # rho_s : numeric, sediment density (kg/m^3). Default 2650.
  # g     : numeric, gravitational acceleration (m/s^2). Default 9.81.
  # kv    : numeric, kinematic viscosity (m^2/s). Default 1e-6.
  #
  # Returns
  #
  # numeric, critical Shields parameter (dimensionless).
  #   Accounts for grain size through the dimensionless grain diameter (d_star),
  #   transitioning from hiding-dominated (fine) to exposure-dominated (coarse)
  #   behaviour.

  Gs     <- (rho_s - rho) / rho
  d_s    <- d_star(d, Gs, g, kv)
  0.3 / (1 + 1.2 * d_s) + 0.055 * (1 - exp(-0.020 * d_s))

  # NOTE: code cleaned up using CLAUDE AI, but performance validated against previous TR'd version of the equation
}



#' Estimating critical dimensionless shear stress
#'
#' Using the equation published by Van Rijn (1984), this function estimates the
#' critical dimensionless shear stress as a function of grain size.
#' @param d median bed sediment grain diameter (m)
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @param rho_s sediment density. Defaults to 2650 kg/m3
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param kv kinematic viscosity of water. Defaults to 1e-6
#' @export
t_crit_vr <- function(d, rho = 1000, rho_s = 2650, g = 9.81, kv = 1e-6) {
  # Estimates the critical Shields parameter for initiation of sediment motion
  # using the five-regime piecewise relationship of Van Rijn (1984).
  #
  # Parameters
  #
  # d     : numeric, grain diameter (m)
  # rho   : numeric, fluid density (kg/m^3). Default 1000.
  # rho_s : numeric, sediment density (kg/m^3). Default 2650.
  # g     : numeric, gravitational acceleration (m/s^2). Default 9.81.
  # kv    : numeric, kinematic viscosity (m^2/s). Default 1e-6.
  #
  # Returns
  #
  # numeric, critical Shields parameter (dimensionless).
  #   Piecewise function of dimensionless grain size (d_star), covering five
  #   hydraulic regimes from viscous (d_star <= 4) to fully rough (d_star > 150).

  Gs  <- (rho_s - rho) / rho
  d_s <- d_star(d, Gs, g, kv)

  if      (d_s <= 4)          0.24 * d_s^-1
  else if (d_s <= 10)         0.14 * d_s^-0.64
  else if (d_s <= 20)         0.04 * d_s^-0.1
  else if (d_s <= 150)        0.013 * d_s^0.29
  else                        0.056

  # NOTE: code cleaned up using CLAUDE AI, but performance validated against previous TR'd version of the equation
}



#' Estimating  dimensionless particle size
#'
#' this function estimates the dimensionless particle size
#' @param d median bed sediment grain diameter (m)
#' @param Gs the submerged specific gravity of sediment. Defaults to 1.65
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param kv kinematic viscosity of water. Defaults to 1e-6
#' @export
d_star <- function(d, Gs = 1.65, g = 9.81, kv = 1e-6){
  d_star <- d * (Gs * g / kv^2)^(1/3)
  return(d_star)
}



#' Estimating bedload transport using Meyer Peter and Muller
#'
#' This function uses a dimensionless form of the traditional Meyer Peter and
#' Muller equation to estimate the dimensionless bedload flux for a given shear
#' stress and characteristic grain size.
#' @param R hydraulic radius for the flow (m)
#' @param S gradient of the water surface slope (m/m)
#' @param d median bed sediment grain diameter (m)
#' @param t_crit critical dimensionless shear stress. Defaults values from Soulsby
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @param rho_s sediment density. Defaults to 2650 kg/m3
#' @export
mpm <- function(R, S, d, t_crit = t_crit_soul(d, rho, rho_s, g), g = 9.81, rho = 1000, rho_s = 2650) {
  # Computes dimensionless bedload transport rate using the
  # Meyer-Peter & Müller (1948) equation.
  #
  # Parameters
  #
  # R      : numeric, hydraulic radius (m)
  # S      : numeric, channel slope (dimensionless)
  # d      : numeric, grain diameter (m), typically D50
  # t_crit : numeric, critical Shields parameter (dimensionless). Default 0.047.
  # g      : numeric, gravitational acceleration (m/s^2). Default 9.81.
  # rho    : numeric, fluid density (kg/m^3). Default 1000.
  # rho_s  : numeric, sediment density (kg/m^3). Default 2650.
  #
  # Returns
  #
  # numeric, dimensionless bedload transport rate (q_star).
  #   Returns 0 when shear stress is below the critical threshold.

  tau <- tau_1D(R, S, g, rho)
  t_s <- t_star(tau, d, g, rho, rho_s)

  if (t_s <= t_crit) return(0)

  8 * (t_s - t_crit)^1.5

  # NOTE: code cleaned up using CLAUDE AI, but performance validated against previous TR'd version of the equation
}



#' Estimating bedload transport using Wong and Parker (2006)
#'
#' This function uses a dimensionless form of a modified version of Meyer Peter
#' and Muller equation to estimate the dimensionless bedload flux for a given
#' shear stress and characteristic grain size. Using the original Meyer Peter and Muller equation will produce higher
#' transport estimates and is more conservative than using this modified version.
#' @param R hydraulic radius for the flow (m)
#' @param S gradient of the water surface slope (m/m)
#' @param d median bed sediment grain diameter (m)
#' @param t_crit critical dimensionless shear stress. Defaults to values from Soulsby
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @param rho_s sediment density. Defaults to 2650 kg/m3
#' @export
wp <- function(R, S, d, t_crit = t_crit_soul(d, rho, rho_s, g), g = 9.81, rho = 1000, rho_s = 2650) {
  # Computes dimensionless bedload transport rate using the
  # Wong & Parker (2006) equation.
  #
  # Parameters
  #
  # R      : numeric, hydraulic radius (m)
  # S      : numeric, channel slope (dimensionless)
  # d      : numeric, grain diameter (m), typically D50
  # t_crit : numeric, critical Shields parameter (dimensionless). Default 0.047.
  # g      : numeric, gravitational acceleration (m/s^2). Default 9.81.
  # rho    : numeric, fluid density (kg/m^3). Default 1000.
  # rho_s  : numeric, sediment density (kg/m^3). Default 2650.
  #
  # Returns
  #
  # numeric, dimensionless bedload transport rate (q_star).
  #   Returns 0 when shear stress is below the critical threshold.

  tau <- tau_1D(R, S, g, rho)
  t_s <- t_star(tau, d, g, rho, rho_s)

  if (t_s <= t_crit) return(0)

  4 * (t_s - t_crit)^1.5

  # NOTE: code cleaned up using CLAUDE AI, but performance validated against previous TR'd version of the equation
}



#' Estimating bedload transport using Van Rijn (1984)
#'
#' This function uses a dimensionless form of Van Rijn's (1984) bedload transport
#' function for use with sand bed rivers.
#' @param R hydraulic radius for the flow (m)
#' @param S gradient of the water surface slope (m/m)
#' @param d median bed sediment grain diameter (m)
#' @param t_crit critical dimensionless shear stress
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @param rho_s sediment density. Defaults to 2650 kg/m3
#' @param kv kinematic viscosity of water. Defaults to 1e-6
#' @export
vr <- function(R, S, d, t_crit = t_crit_vr(d, rho, rho_s, g, kv),
               g = 9.81, rho = 1000, rho_s = 2650, kv = 1e-6) {
  # Computes dimensionless bedload transport rate using the Van Rijn (1984) approach.
  #
  # Parameters
  #
  # R      : numeric, hydraulic radius (m)
  # S      : numeric, channel slope (dimensionless)
  # d      : numeric, grain diameter (m), typically D50
  # t_crit : numeric, critical Shields parameter (dimensionless).
  #            Defaults to the value returned by t_crit_vr() using the
  #            physical constants below.
  # g      : numeric, gravitational acceleration (m/s^2). Default 9.81.
  # rho    : numeric, fluid density (kg/m^3). Default 1000.
  # rho_s  : numeric, sediment density (kg/m^3). Default 2650.
  # kv     : numeric, kinematic viscosity (m^2/s). Default 1e-6.
  #
  # Returns
  #
  # numeric, dimensionless bedload transport rate (q_star).
  #   Returns 0 when shear stress is below the critical threshold.

  tau   <- tau_1D(R, S, g, rho)
  t_s   <- t_star(tau, d, g, rho, rho_s)
  Gs    <- (rho_s - rho) / rho
  d_s   <- d_star(d, Gs, g, kv)

  # return zero transport below the critical threshold rather than a negative value
  if (t_s <= t_crit) return(0)

  (0.053 / d_s^0.3) * (t_s / t_crit - 1)^2.1

  # NOTE: code cleaned up using CLAUDE AI, but performance validated against previous TR'd version of the equation
}



#' Estimating bedload transport using Yalin (1963)
#'
#' This function uses a dimensionless form of Yalin's (1963) bedload transport
#' function for use with sand bed rivers and gravel bed rivers.
#' @param R hydraulic radius for the flow (m)
#' @param S gradient of the water surface slope (m/m)
#' @param d median bed sediment grain diameter (m)
#' @param t_crit critical dimensionless shear stress, calculated using Soulsby (1997)
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @param rho_s sediment density. Defaults to 2650 kg/m3
#' @param kv kinematic viscosity of water. Defaults to 1e-6
#' @export
yln <- function(R, S, d, t_crit = t_crit_soul(d), g = 9.81, rho = 1000, rho_s = 2650, kv = 1e-6){
  # Computes dimensionless bedload transport rate using the Yalin (1963) equation.
  #
  # Parameters
  #
  # R      : numeric, hydraulic radius (m)
  # S      : numeric, channel slope (dimensionless)
  # d      : numeric, grain diameter (m), typically D50
  # t_crit : numeric, critical Shields parameter (dimensionless).
  #            Defaults to the value returned by t_crit_soul().
  # g      : numeric, gravitational acceleration (m/s^2). Default 9.81.
  # rho    : numeric, fluid density (kg/m^3). Default 1000.
  # rho_s  : numeric, sediment density (kg/m^3). Default 2650.
  # kv     : numeric, kinematic viscosity (m^2/s). Default 1e-6. Unused directly
  #            but retained for consistency with vr() signature.
  #
  # Returns
  #
  # numeric, dimensionless bedload transport rate (q_star).
  #   Returns 0 when shear stress is below the critical threshold.

  tau   <- tau_1D(R, S, g, rho)
  t_s   <- t_star(tau, d, g, rho, rho_s)
  Gs    <- (rho_s - rho) / rho

  if (t_s <= t_crit) return(0)

  r     <- t_s / t_crit - 1
  sigma <- 2.45 * sqrt(t_crit) / (Gs + 1)^0.4

  0.635 * r * sqrt(t_s) * (1 - (1 / (sigma * r)) * log(1 + sigma * r))
  # NOTE: code cleaned up using CLAUDE AI, but performance validated against previous TR'd version of the equation
}



#' Estimating bedload transport using Einstein and brown (1950)
#'
#' This function uses a dimensionless form of Einstein and brown (1950) bedload
#' transport function for use with sand bed rivers and gravel bed rivers.
#' @param R hydraulic radius for the flow (m)
#' @param S gradient of the water surface slope (m/m)
#' @param d median bed sediment grain diameter (m)
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @param rho_s sediment density. Defaults to 2650 kg/m3
#' @param kv kinematic viscosity of water. Defaults to 1e-6
#' @export
eb <- function(R, S, d, g = 9.81, rho = 1000, rho_s = 2650, kv = 1e-6) {
  # Computes dimensionless bedload transport rate using the
  # Einstein-Brown equation (Brown, 1950).
  #
  # Parameters
  #
  # R      : numeric, hydraulic radius (m)
  # S      : numeric, channel slope (dimensionless)
  # d      : numeric, grain diameter (m), typically D50
  # g      : numeric, gravitational acceleration (m/s^2). Default 9.81.
  # rho    : numeric, fluid density (kg/m^3). Default 1000.
  # rho_s  : numeric, sediment density (kg/m^3). Default 2650.
  # kv     : numeric, kinematic viscosity (m^2/s). Default 1e-6.
  #
  # Returns
  #
  # numeric, dimensionless bedload transport rate (q_star).
  #   Two regimes: exponential form below t_s = 0.19, power law above.
  #   No critical threshold — transport is non-zero for any positive shear stress.

  tau <- tau_1D(R, S, g, rho)
  t_s <- t_star(tau, d, g, rho, rho_s)
  Gs  <- (rho_s - rho) / rho
  d_s <- d_star(d, Gs, g, kv)
  K   <- sqrt(2/3 + 36 / d_s^3) - sqrt(36 / d_s^3)

  if (t_s >= 0.19) {
    40 * K * t_s^3
  } else {
    (K * exp(-0.391 / t_s)) / 0.465
  }

  # NOTE: code cleaned up using CLAUDE AI, but performance validated against previous TR'd version of the equation
}

#' Estimating bed material transport using Eaton and Church (2011)
#'
#' This function uses a dimensionless stream power equation from Eaton and Church
#' (2011) to estimate bed material load transport in sand bed rivers and gravel
#' bed rivers.
#' @param R hydraulic radius for the flow (m)
#' @param S gradient of the water surface slope (m/m)
#' @param d median bed sediment grain diameter (m)
#' @param d84 84th percentile diameter of bed sediment (m). Defaults to 2 * d
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @param rho_s sediment density. Defaults to 2650 kg/m3
#' @param U mean flow velocity. Defaults to flows calculated using Ferguson (2007)
#' @param a1 constant fit by Ferguson (2007). Defaults to 6.5
#' @param a2 constant fit by Ferguson (2007). Defaults to 2.5
#' @param kv kinematic viscosity of water. Defaults to 1e-6
#' @param t_crit critical dimensionless shear stress, calculated using Soulsby (1997)
#' @export
ec <- function(R, S, d,
               d84    = 2 * d,
               g      = 9.81,
               rho    = 1000,
               rho_s  = 2650,
               kv     = 1e-6,
               a1     = 6.5,
               a2     = 2.5,
               U      = u_ferg(R, S, d84, a1, a2, g),
               t_crit = t_crit_soul(d, rho, rho_s, g, kv)) {
  # Computes dimensionless bedload transport rate using the Eaton & Church (2011)
  # stream power efficiency approach.
  #
  # Parameters
  #
  # R      : numeric, hydraulic radius (m)
  # S      : numeric, channel slope (dimensionless)
  # d      : numeric, median grain diameter (m), D50
  # d84    : numeric, 84th percentile grain diameter (m). Default 2 * d.
  # g      : numeric, gravitational acceleration (m/s^2). Default 9.81.
  # rho    : numeric, fluid density (kg/m^3). Default 1000.
  # rho_s  : numeric, sediment density (kg/m^3). Default 2650.
  # kv     : numeric, kinematic viscosity (m^2/s). Default 1e-6.
  # a1     : numeric, Ferguson resistance coefficient, smooth end. Default 6.5.
  # a2     : numeric, Ferguson resistance coefficient, rough end. Default 2.5.
  # U      : numeric, mean flow velocity (m/s). Defaults to u_ferg(R, S, d84, a1, a2, g).
  # t_crit : numeric, critical Shields parameter (dimensionless).
  #            Defaults to t_crit_soul(d, rho, rho_s, g, kv).
  #
  # Returns
  #
  # numeric, dimensionless bedload transport rate (q_star).
  #   Returns 0 when dimensionless stream power is zero or negative.

  tau     <- tau_1D(R, S, g, rho)
  Gs      <- (rho_s - rho) / rho

  # dimensionless stream power
  omega_s <- tau * U / (rho * (g * Gs * d)^1.5)

  if (omega_s <= 0) return(0)

  # critical depth and reference resistance at transport initiation
  d_crit  <- t_crit * Gs * d / S
  Res     <- u_ferg(d_crit, S, d84, a1, a2, g) / sqrt(g * d_crit * S)

  # reference dimensionless stream power at threshold
  om_crit <- Res * t_crit^1.5

  # transport efficiency (Eaton & Church 2011)
  E_star  <- (0.92 - 0.25 * sqrt(om_crit / omega_s))^9
  q_star  <- E_star * omega_s

  pmax(q_star, 0)

  # NOTE: code cleaned up using CLAUDE AI, but performance validated against previous TR'd version of the equation
}

#' Estimating bed material transport using Recking (2013)
#'
#' This function uses a dimensionless shear stress equation from Recking
#' based on the D84 of the bed surface, and is the version for use with
#' field data.
#' @param R hydraulic radius for the flow (m)
#' @param S gradient of the water surface slope (m/m)
#' @param d median bed sediment grain diameter (m)
#' @param d84 84th percentile diameter of bed sediment (m). Defaults to 2 * d
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @param rho_s sediment density. Defaults to 2650 kg/m3
#' @param rp option to select equation suited to riffle pool morphology, default is FALSE
#' @export
rk <- function(R, S, d, d84 = 2 * d, g = 9.81, rho = 1000, rho_s = 2650, rp = FALSE) {
  # Computes dimensionless bedload transport rate using the
  # Recking (2013) equation.
  #
  # Parameters
  #
  # R      : numeric, hydraulic radius (m)
  # S      : numeric, channel slope (dimensionless)
  # d      : numeric, grain diameter (m), typically D50
  # d84    : numeric, 84th percentile grain diameter (m). Default 2 * d.
  # g      : numeric, gravitational acceleration (m/s^2). Default 9.81.
  # rho    : numeric, fluid density (kg/m^3). Default 1000.
  # rho_s  : numeric, sediment density (kg/m^3). Default 2650.
  # rp     : logical, if TRUE uses the reach-averaged (riffle-pool) form of the
  #            critical Shields parameter; if FALSE uses the plane-bed form.
  #            Default FALSE.
  #
  # Returns
  #
  # numeric, dimensionless bedload transport rate (q_star).
  #   No critical threshold — transport is non-zero for any positive shear stress.
  #   Note: Shields stress is computed using d84, not d.

  tau  <- tau_1D(R, S, g, rho)
  t_s  <- t_star(tau, d84, g, rho, rho_s)

  t_sm <- if (rp) {
    (5 * S + 0.06) * (d84 / d)^(4.4 * sqrt(S) - 1.5)
  } else {
    1.5 * S^0.75
  }

  14 * t_s^2.5 / (1 + (t_sm / t_s)^4)

  # NOTE: code cleaned up using CLAUDE AI, but performance validated against previous TR'd version of the equation
}




#' Compare transport functions
#'
#' This function compares the transport functions for a specified grain size
#' and channel gradient over a reasonable range of shear stress values using the
#' default values for all functions
#' @param d median bed sediment grain diameter (m)
#' @param S gradient of the water surface slope (m/m)
#' @param d84 84th percentile diameter of bed sediment (m). Defaults to 2 * d
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @param rho_s sediment density. Defaults to 2650 kg/m3
#' @export
compare_fun <- function(d, S, d84 = 2 * d, g = 9.81, rho = 1000, rho_s = 2650) {
  # Plots dimensionless bedload transport as a function of bed shear stress for
  # seven transport equations, over a range of tau from 0.25 to 5 times the
  # reference critical shear stress (tau_c = 0.05 * g * (rho_s - rho) * d).
  #
  # Parameters
  #
  # d     : numeric, median grain diameter (m), D50
  # S     : numeric, channel slope (dimensionless)
  # d84   : numeric, 84th percentile grain diameter (m). Default 2 * d.
  #           Used by ec() and rk().
  # g     : numeric, gravitational acceleration (m/s^2). Default 9.81.
  # rho   : numeric, fluid density (kg/m^3). Default 1000.
  # rho_s : numeric, sediment density (kg/m^3). Default 2650.
  #
  # Returns
  #
  # NULL, invisibly. Called for its side effect of producing a plot.

  # --- shear stress and hydraulic radius range ---
  tau_critical <- 0.05 * g * (rho_s - rho) * d
  tau_range    <- seq(0.25 * tau_critical, 5 * tau_critical, length.out = 200)
  R_range      <- tau_range / (g * rho * S)

  # --- evaluate all transport equations over R_range ---
  compute <- function(fun, grain = d, ...) {
    as.numeric(lapply(R_range, FUN = fun, S = S, d = grain, g = g,
                      rho = rho, rho_s = rho_s, ...))
  }

  qs_mpm <- compute(mpm)
  qs_wp  <- compute(wp)
  qs_smt  <- compute(smt)
  qs_yln <- compute(yln)
  qs_ec  <- compute(ec,  d84 = d84)
  qs_eb  <- compute(eb)
  qs_vr  <- compute(vr)
  qs_rk  <- compute(rk,  d84 = d84, rp = TRUE)

  # --- plot: filter to positive mpm values for axis range ---
  filt <- qs_mpm > 0
  plot(tau_range[filt], qb_conversion(qs_mpm[filt], d),
       type  = "l",
       log   = "xy",
       main  = paste("grain size =", d * 1000, "mm, gradient =", S, "m/m"),
       col   = "red",
       lty   = 3,
       xlab  = expression(tau ~ (Pa)),
       ylab  = expression("bedload flux" ~ (m^3/s/m)),
       cex.main = 0.7)

  lines(tau_range, qb_conversion(qs_wp,  d),    col = "darkorange", lty = 3)
  lines(tau_range, qb_conversion(qs_smt,  d),    col = "darkorange", lty = 1)
  lines(tau_range, qb_conversion(qs_yln, d),    col = "blue")
  lines(tau_range, qb_conversion(qs_ec,  d),    col = "purple",     lty = 2)
  lines(tau_range, qb_conversion(qs_eb,  d),    col = "darkred")
  lines(tau_range, qb_conversion(qs_vr,  d),    col = "black",      lty = 2)
  lines(tau_range, qb_conversion(qs_rk,  d84),  col = "darkgreen",  lty = 3)

  legend("bottomright",
         inset  = 0.02,
         legend = c("Eaton & Church", "Einstein-Brown", "Meyer-Peter & Müller",
                    "Recking", "Van Rijn", "Wong & Parker", "Smart", "Yalin"),
         col    = c("purple", "darkred", "red", "darkgreen", "black", "darkorange","darkorange", "blue"),
         lty    = c(2, 1, 3, 3, 2, 3, 1, 1),
         pch    = NA,
         cex    = 0.7)

  invisible(NULL)
  # NOTE: code cleaned up using CLAUDE AI, but performance validated against previous TR'd version of the equation
}


#' Load cross section data
#'
#'This function loads profile data that is saved in a text file (XYZ format) and
#'calculates the distance along the profile from the first data point in the
#'file. The function assumes that the first row is header information, and that
#'the data in column 1 correspond to x coordinates (or eastings), column 2
#'corresponds to y coordinates (northings), and column 3 corresponds to elevation.
#' @param xs_nm path to a csv file with 3 columns containing the profile data
#' @param dx optional dist. increment to resample section at regular intervals (m)
#' @param rev option to reverse the direction in which distance is measured
#' @param trim option to identify left & right banks and trim section
#' @param mkplt option to make a plot of the data (original & resampled/trimed)
#' @export
load_section <- function(xs_nm, dx = NA, rev = FALSE, trim = FALSE, mkplt = FALSE){
  #load data and calculate distance along the profile
  xs_data <- utils::read.csv(file = xs_nm, skip = 1)
  colnames(xs_data) <- c("E", "N", "Z")
  xs_data$dist <- c(0, cumsum(sqrt(diff(xs_data$N)^2 + diff(xs_data$E)^2)))
  if(rev){
    xs_data$dist <- max(xs_data$dist, na.rm = T) - xs_data$dist
  }
  if(mkplt){
    graphics::plot(xs_data$dist, xs_data$Z,
                   type = "l",
                   col = "darkgreen",
                   asp = 2,
                   lty = 1,
                   lwd = 2,
                   xlab = "distance from left bank (m)",
                   ylab = "elevation above sea level (m)")
    graphics::legend("bottomright",
                     inset = c(0.01, 0.02),
                     legend = c("Original", "Trimmed", "Resampled"),
                     col = c("darkgreen", "darkorange", "darkred"),
                     lty = c(1, 1, NA),
                     pch = c(NA, NA, 19),
                     cex = 0.7)
  }
  if(trim){
    #find the thalweg of the channel, where long profile intersects the cross section
    center <- which(xs_data$Z == min(xs_data$Z, na.rm = T))[1]

    #find left bank & delete ground points left of bank
    filt <- xs_data$dist < xs_data$dist[center]
    leftbank <- which(xs_data$Z == max(xs_data$Z[filt], na.rm = T))[1]
    filt <- xs_data$dist < xs_data$dist[leftbank]
    xs_data$Z[filt] <- NA  #delete topography outside the main channel

    #find the right bank
    filt <- xs_data$dist > xs_data$dist[center]
    rightbank <- which(xs_data$Z == max(xs_data$Z[filt], na.rm = T))[1]
    filt <- xs_data$dist > xs_data$dist[rightbank]
    xs_data$Z[filt] <- NA #delete topography outside the main channel
    if(mkplt){
      graphics::lines(xs_data$dist, xs_data$Z, lty = 1, lwd = 2, col = "darkorange")
    }
  }
  #drop the data outside the banks
  xs_data <-xs_data[!is.na(xs_data$Z),]

  if(is.na(dx)){
    return(data.frame(E = xs_data$E,
                      N = xs_data$N,
                      Z = xs_data$Z,
                      dist = xs_data$dist))
  } else{
    sample_d <- seq(from = min(xs_data$dist),
                    to = max(xs_data$dist),
                    by = dx)
    sample_z <- stats::approx(x = xs_data$dist,
                              y = xs_data$Z,
                              xout = sample_d)[[2]]
    sample_x <- stats::approx(x = xs_data$dist,
                              y = xs_data$E,
                              xout = sample_d)[[2]]
    sample_y <- stats::approx(x = xs_data$dist,
                              y = xs_data$N,
                              xout = sample_d)[[2]]
    if(mkplt){
      graphics::points(sample_d, sample_z, pch = 19, cex = 0.6, col = "darkred")
    }
    return(data.frame(E = sample_x,
                      N = sample_y,
                      Z = sample_z,
                      dist = sample_d))
  }
}



#' Calculate flow depths
#'
#'this function uses a water surface elevation (WSE) and the location index for
#'the channel thalweg to calculate the depths in the main channel.  disconnected
#'areas outside the main channel that are below the WSE are located and assigned
#'a depth of NA. This is particularly important on a fan
#' @param xs_data a data frame containing $dist (distances) and $z (elevations)
#' @param wse elevation of the water surface of interest (m)
#' @param mkplt logical option to plot the data, default is FALSE
#' @export
calc_depths <- function(xs_data, wse, mkplt = F){
  thalweg <- which(xs_data$Z == min(xs_data$Z))[1]
  depth <- wse - xs_data$Z
  depth[depth < 0] <- NA
  for (j in seq_along(xs_data$dist)){
    if(j < thalweg){
      if(sum(is.na(depth[j:thalweg])) > 0){
        depth[j] <- NA  #set all depths separated from main channel by a dry patch to zero
      }
    }
    if(j > thalweg){
      if(sum(is.na(depth[thalweg:j])) > 0){
        depth[j] <- NA  #set all depths separated from main channel by a dry patch to zero
      }
    }
  }
  if(mkplt){
    graphics::plot(xs_data$dist,
                   xs_data$Z,
                   type = "l",
                   col = "darkgreen",
                   asp = 2,
                   main = "locations where water depths have been calculated",
                   xlab = "distance from left bank (m)",
                   ylab = "elevation above sea level (m)")
    graphics::abline(h = wse, lty = 2, col = "blue")
    # Filter out NA values for plotting points
    valid_points <- !is.na(depth)
    graphics::points(xs_data$dist[valid_points], xs_data$Z[valid_points],
                     pch = 19, col = "blue")
  }
  return(depth)
}



#' Estimate flow hydraulics for range of flood stages
#'
#'this function estimates hydraulic properties at regular WSE elevation intervals
#'these data can be used to interpolate values for target flow values.
#'
#' @param xs_data a data frame containing $dist (distances) and $z (elevations)
#' @param S gradient of the water surface slope (m/m)
#' @param d median bed sediment grain diameter (m)
#' @param d84 84th percentile diameter of bed sediment (m). Defaults to 2 * d
#' @param ymax the maximum flow depth to model, measured at the thalweg, defaults to 1.5 m
#' @param dx the increment by which to change stage. defaults to 0.05 m
#' @param t_crit critical dimensionless shear stress, calculated using Soulsby (1997)
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @param rho_s sediment density. Defaults to 2650 kg/m3
#' @param a1 constant fit by Ferguson (2007). Defaults to 6.5
#' @param a2 constant fit by Ferguson (2007). Defaults to 2.5
#' @param kv kinematic viscosity of water. Defaults to 1e-6
#' @param fun choice of bedload function: 'eb','ec', 'rk', 'mpm', 'vr','wp', 'yln'. Defaults to 'ec' ('mpm', 'wp', 'eb', 'yln', 'rk')
#' @param rp option for 'rk' select equation suited to riffle pool morphology, default is FALSE
#' @export
estimate_hydraulics <- function(xs_data,
                                S,
                                d,
                                d84 = 2 * d,
                                ymax = 3.5,
                                dx = 0.05,
                                t_crit = t_crit_soul(d),
                                g = 9.81,
                                rho = 1000,
                                rho_s = 2650,
                                a1 = 6.5,
                                a2 = 2.5,
                                kv = 1e-6,
                                fun = "ec",
                                rp = FALSE){
  Gs <- (rho_s - rho) / rho
  #specify the range of target WSE, and calculate dx and dp
  y_range <- seq(0, ymax, dx)  #make calculations for a range of depths for each section
  dx <- (c(0,diff(xs_data$dist)) + c(diff(xs_data$dist),0) ) /2
  px <- sqrt(diff(xs_data$dist)^2 + diff(xs_data$Z)^2)
  dp <- (c(0, px) + c(px, 0)) / 2
  invert <- min(xs_data$Z, na.rm = T)
  thalweg <- which(xs_data$Z == min(xs_data$Z, na.rm = T))[1]
  xs_hyd <- data.frame(flow = array(NA, length(y_range)),
                       area = array(NA, length(y_range)),
                       width = array(NA, length(y_range)),
                       perimeter = array(NA, length(y_range)),
                       velocity = array(NA, length(y_range)),
                       transport = array(NA, length(y_range)),
                       wse = array(NA, length(y_range)),
                       t_star = array(NA, length(y_range)))

  for(j in seq_along(y_range)){
    wse <- invert + y_range[j]
    depths <- calc_depths(xs_data, wse)
    #plot(depths)
    A <- sum(depths * dx, na.rm = T)
    W <- sum(dx[!is.na(depths)], na.rm = T)
    P <- sum(dp[!is.na(depths)], na.rm = T)
    R <- A / P
    U <- u_ferg(R, S, d84, a1, a2, g)
    tau <- tau_1D(R, S, g, rho)
    t_star <- t_star(tau, d, g, rho, rho_s)
    if(fun == "eb"){
      q_star <- eb(R, S, d, g, rho, rho_s, kv)
    }
    if(fun == "ec"){
      q_star <- ec(R, S, d, d84, g, rho, rho_s, U, t_crit)
    }
    if(fun == "rk"){
      d84_transp <- d * 2  #recalculate based on d
      q_star <- rk(R, S, d, d84_transp, g, rho, rho_s, rp)
    }
    if(fun == "mpm"){
      q_star <- mpm(R, S, d, t_crit, g, rho, rho_s)
    }
    if(fun == "wp"){
      q_star <- wp(R, S, d, t_crit, g, rho, rho_s)
    }
    if(fun == "vr"){
      q_star <- vr(R, S, d, t_crit, g, rho, rho_s, kv)
    }
    if(fun == "yln"){
      q_star <- yln(R, S, d, t_crit, g, rho, rho_s, kv)
    }
    if(fun == "rk"){
      Qb <- qb_conversion(q_star, d84_transp, g, Gs) * W
    } else {
      Qb <- qb_conversion(q_star, d, g, Gs) * W
    }
    Q <- A * U
    xs_hyd[j,] <- c(Q, A, W, P, U, Qb, wse, t_star)
  }
  return(xs_hyd)
}



#' Estimate transport capacity for a hydrograph
#'
#'this function uses a data frame created by the estimate_hydraulics() function
#'to approximate the volume of bed material that could be transported during a
#'hydrograph. It also estimates the geomorphically active volume of transport
#'based on a threshold dimensionless shear stress.
#' @param xs_hyd a data frame created by estimate_hydraulics()
#' @param Q a vector containing the flows that define the hydrograph
#' @param t a vector that contains the corresponding times for each Q value (hours)
#' @param tc the threshold dimensionless shear stress for channel change
#' @param mkplt option to make a plot of the data, defaults to FALSE
#' @export
estimate_capacity <- function(xs_hyd, Q, t, tc = 0.08, mkplt = FALSE){
  #this function uses the output from estimate_hydraulics to approximate the
  #transport rate, peak width, area, velocity, and equivalent manning n value
  #and integrates transport capacity over a  hydrograph produced by HEC-HMS

  min2sec <- 60
  Qmax <- max(Q)  #extract and save the peak flow
  Qt <- stats::approx(x = t, y = Q, xout = seq(from = min(t), to = max(t), by = 1/60)) #resample to 1 hour
  Qbt <- stats::approx(x = xs_hyd$flow,  #estimate transport every minute
                       y = xs_hyd$transport,
                       xout = Qt$y)
  t_star <- stats::approx(x = xs_hyd$flow,
                          y = xs_hyd$t_star,
                          xout = Qt$y)
  filt <- t_star$y >= tc  #only those values where d50 is fully mobile
  vol_transp <- sum(Qbt$y * min2sec, na.rm = T)  #add up all transport values
  vol_active <- sum(Qbt$y[filt] * min2sec, na.rm = T)#sum transport capacity once bed fully mobile
  width <- stats::approx(x = xs_hyd$flow,  #find the channel width at peak flow
                         y = xs_hyd$width,
                         xout = Qmax)[[2]]
  vel <- stats::approx(x = xs_hyd$flow,  #find the peak velocity
                       y = xs_hyd$velocity,
                       xout = Qmax)[[2]]
  area <- stats::approx(x = xs_hyd$flow,  #find the peak area
                        y = xs_hyd$area,
                        xout = Qmax)[[2]]
  wse <- stats::approx(x = xs_hyd$flow,  #find the peak area
                       y = xs_hyd$wse,
                       xout = Qmax)[[2]]
  d <- area / width
  output <- data.frame(Qmax, vol_transp, vol_active, width, vel, area, wse)

  if(mkplt){
    graphics::par(mar = c(5,5,1,1))
    graphics::plot(Qt$x, Qbt$y,
                   type = "l",
                   col = "blue",
                   xlab = "time (hrs)",
                   ylab = expression(transport~rate~(m^3/s))
    )
    graphics::lines(Qt$x[filt], Qbt$y[filt],
                     pch = 19,
                     cex = 0.5,
                     lwd = 3,
                     col = grDevices::rgb(1,0,0,0.5)
    )
    graphics::legend("topright",
                     inset = c(0.01, 0.02),
                     legend = c("total", "active"),
                     col = c(rgb(0,0,1), rgb(1,0,0,0.5)),
                     pch = c(NA, 19),
                     lty = c(1,NA)
    )
  }

  return(output)
}

####NEW FUNCTIONS ADDED MAY 2025 ############

#' Estimating bedload transport using Smart (1984)
#'
#' This function uses the Smart (1984) equation to estimate the dimensionless
#' bedload flux for steep channels (slopes 3-25 percent). The equation extends
#' Meyer-Peter & Müller to steep slopes by incorporating an explicit slope
#' term and a grain size non-uniformity factor. It is appropriate for coarse
#' sediment (d > 0.4 mm) in steep gravel- or cobble-bed channels.
#'
#' @param R hydraulic radius for the flow (m)
#' @param S gradient of the water surface slope (m/m)
#' @param d median bed sediment grain diameter (m)
#' @param d90_30 ratio of the 90th and the 30th percentile of the gsd. Defaults to 3.4 (from a log normal distribution)
#' @param d84 84th percentile diameter of bed sediment (m). Defaults to 2 * d
#' @param t_crit critical dimensionless shear stress. Defaults to 0.047
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @param rho_s sediment density. Defaults to 2650 kg/m3
#' @param a1 constant fit by Ferguson (2007). Defaults to 6.5
#' @param a2 constant fit by Ferguson (2007). Defaults to 2.5
#' @param U mean flow velocity. Defaults to flows calculated using Ferguson (2007)

#' @export
smt <- function(R, S, d,
                d90_30 = 3.4,
                d84    = 2*d,
                t_crit = 0.047,
                g      = 9.81,
                rho    = 1000,
                rho_s  = 2650,
                a1     = 6.5,
                a2     = 2.5,
                U      = u_ferg(R, S, d84, a1, a2, g)) {
  # Computes dimensionless bedload transport rate using the
  # Meyer-Peter & Müller (1948) equation.
  #
  # Parameters
  #
  # R      : numeric, hydraulic radius (m)
  # S      : numeric, channel slope (dimensionless)
  # d      : numeric, grain diameter (m), typically D50
  # t_crit : numeric, critical Shields parameter (dimensionless). Default 0.047.
  # g      : numeric, gravitational acceleration (m/s^2). Default 9.81.
  # rho    : numeric, fluid density (kg/m^3). Default 1000.
  # rho_s  : numeric, sediment density (kg/m^3). Default 2650.
  #
  # Returns
  #
  # numeric, dimensionless bedload transport rate (q_star).
  #   Returns 0 when shear stress is below the critical threshold.

  tau <- tau_1D(R, S, g, rho)
  t_s <- t_star(tau, d, g, rho, rho_s)

  if (t_s <= t_crit) return(0)
  C <- U / (sqrt(g * R * S))  #calculate the flow resistance term

  #calculate Smart's (1984) dimensionless transport rate using Eq. 8 from
  #the paper
  4 * (d90_30)^0.2 * S^0.6 * C * sqrt(t_s) * (t_s - t_crit)
  #Not yet formally reviewed, but performs very similarly to Wong and Parker at 0.02 m/m & 50 mm d50
}

#' Integrating  Wilcock and Crowe (2003) over a hydrograph
#'
#' This function integrates the Wilcock and Crowe transport function over a
#' hydrograph for a given cross section to estimate the total fractional
#' transport volumes for the flood event. Output includes the inputted GSD
#' information for each size class (Di, cpf, Dm, Fi), as well as the volumetric
#' transport rate (Qbi), and the fraction of the computed bedload in each size
#' class (Li)
#' @param hydrograph a data frame created by scs_hydrograph() (Q in m3/s and t in hours)
#' @param xs_data a data frame containing topography expressed as $dist (distances) and $z (elevations)
#' @param S gradient of the water surface slope (m/m)
#' @param gsd a data frame describing the grain size distribution created with function wolman_gsd()
#' @param Fsand proportion of the surface that is covered by sand sized sediment (affects entrainment, typically varies from 0 to 0.2)
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @param rho_s sediment density. Defaults to 2650 kg/m3
#' @export
fractional_capacity <-function(hydrograph, xs_data, S, gsd, Fsand, g = 9.81, rho = 1000, rho_s = 2650){
  #define variables
  hours2sec <- 60 *60
  dt <- mean(diff(hydrograph$x))
  # extract some reference grain sizes and perform 1D hydraulics
  d50 <- approx(x = gsd$cpf, y = gsd$Di, 0.5)[[2]]
  d84 <- approx(x = gsd$cpf, y = gsd$Di, 0.84)[[2]]
  # estimate the flow depth for the peak discharge for estimate_hydraulics func
  dmax <- 0.3 * max(hydrograph$y) ^ 0.33
  xs_hyd <- estimate_hydraulics(xs_data,
                                S,
                                d = (d50 / 1000),
                                d84 = (d84 / 1000),
                                ymax = dmax * 3,
                                fun = "wp")

  # estimate the hydraulic radius and channel width for all flows in the hydrograph
  xs_hyd$R <- xs_hyd$area / xs_hyd$perimeter
  R <- approx(x = xs_hyd$flow,
              y = xs_hyd$R,
              xout = hydrograph$y,
              rule = 2)[[2]]
  W <- approx(x = xs_hyd$flow,
              y = xs_hyd$width,
              xout = hydrograph$y,
              rule = 2)[[2]]

  # set up the storage matrix
  r <- length(gsd$dm)
  c <- length(R) + 1
  out_array <- array(data = NA, dim = c(r,c))
  out_array[,1] <- gsd$dm

  # loop over all flows in the hydrograph
  for(i in seq_along(R)){
    tmp<- wc(R[i], S, gsd, Fsand, g, rho, rho_s)
    out_array[, (i+1)] <- tmp$qbi * W[i] * dt * hours2sec
  }
  # sum over all timesteps
  gsd$Qbi <- c(rowSums(out_array[,-1]) )
  gsd$Li <- gsd$Qbi / sum(gsd$Qbi, na.rm = T)
  return(gsd)
}



#' Estimating fractional bedload transport using Wilcock and Crowe (2003) in m3/second
#'
#' This function uses the Wilcock and Crowe transport function on standard data
#' from a surface Wolman count to estimate the fraction-specific transport rate
#' for a given shear stress. Output is a data frame containing inputted GSD
#' information for each size class (Di, cpf, Dm, Fi), as well as the dimensionless
#' transport rate (Wi), the volumetric transport rate per unit width (qbi), and
#' the fraction of the computed bedload in each size class (Li)
#' @param R hydraulic radius for the flow (m)
#' @param S gradient of the water surface slope (m/m)
#' @param gsd a data frame describing the grain size distribution created with function wolman_gsd()
#' @param Fsand proportion of the surface that is covered by sand sized sediment (affects entrainment, typically varies from 0 to 0.2)
#' @param counts number of stones in each size class size class bounds (!!!must start with a 0!!!)
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @param rho_s sediment density. Defaults to 2650 kg/m3
#' @export
wc <- function(R, S, gsd, Fsand, g = 9.81, rho = 1000, rho_s = 2650) {
  # Computes fractional bedload transport rates using the
  # Wilcock & Crowe (2003) surface-based mixture transport equation.
  #
  # Parameters
  #
  # R      : numeric, hydraulic radius (m)
  # S      : numeric, channel slope (dimensionless)
  # gsd    : data frame describing the surface grain size distribution, with columns:
  #            Di  - grain size (mm)
  #            Fi  - fraction of surface area in each size class
  #            cpf - cumulative proportion finer (0-1)
  #            dm  - geometric mean grain size for each fraction (mm)
  # Fsand  : numeric, proportion of sand on the bed surface (0-1).
  #            Used to compute the reference Shields stress via t_ref_star().
  #            Can be estimated from bank erosion probability or proximity to
  #            sandy sediment sources; replaces the commented-out approx() call.
  # g      : numeric, gravitational acceleration (m/s^2). Default 9.81.
  # rho    : numeric, fluid density (kg/m^3). Default 1000.
  # rho_s  : numeric, sediment density (kg/m^3). Default 2650.
  #
  # Returns
  #
  # A copy of gsd with three additional columns:
  #   Wi  - dimensionless transport parameter W* for each size fraction
  #   qbi - volumetric transport rate per unit width (m^2/s) for each fraction
  #   Li  - fractional contribution of each size class to total transport (0-1)

  tau    <- tau_1D(R, S, g, rho)
  u_star <- sqrt(tau / rho)

  # --- grain size statistics ---
  gsd <- gsd[!duplicated(gsd$cpf), ]
  gsd <- gsd[order(gsd$cpf), ]
  d50 <- approx(x = gsd$cpf, y = gsd$Di, xout = 0.5, ties = "ordered")$y  # surface median (mm)

  # --- size-specific critical shear stress ---
  shields_50 <- t_ref_star(Fsand)
  b          <- 0.67 / (1 + exp(1.5 - gsd$dm / d50))
  tr_i       <- g * shields_50 * (rho_s - rho) * (d50 / 1000) * (gsd$dm / d50)^b

  # --- dimensionless transport parameter W* (two-regime, Wilcock & Crowe 2003) ---
  phi  <- tau / tr_i
  W_i  <- ifelse(phi < 1.35,
                 0.002 * phi^7.5,
                 14 * (1 - 0.895 / sqrt(phi))^4.5)

  # --- volumetric transport rate per unit width for each fraction ---
  qb_i <- W_i * gsd$Fi * u_star^3 / (g * (rho_s / rho - 1))

  # --- fractional contribution to total transport ---
  L_i  <- qb_i / sum(qb_i)

  # --- append results to gsd and return ---
  df      <- gsd
  df$Wi   <- W_i
  df$qbi  <- qb_i
  df$Li   <- L_i
  df
}


#' Estimating reference dimensionless shear stress for Wilcock & Crowe (2003)
#'
#' Using the equation published by Wilcock & Crow, this function estimates the
#' reference dimensionless shear stress based on the fraction of sand on the bed.
#' @param Fsand Proportion of sand on the bed surface
#' @export
t_ref_star <- function(Fsand){
  ts_ref <- 0.021 + 0.015 * exp(-20 * Fsand)
  return(ts_ref)
}



#' Processing Wolman count data into grain size distribution
#'
#' This function uses standard data from a surface Wolman count to estimate
#' cumulative percent finer for each size, the fraction of the distribution in
#' a given size class, and the median size of each size class. Output is a data
#' frame containing the upper limit of each size class (Di), the cumulative
#' proportion finer than the size listed in Di (cpf), the median size for each
#' size class (Dm), and the fraction of the distribution in each class (Fi)
#' @param Di a list of grain sizes representing the upper bounds of the size classes
#' @param counts the number of stones sampled in the szie class defined by the upper bound values, Di
#' @export
wolman_gsd <- function(Di, counts){
  #check that the input data are listed in ascending order
  if(Di[1] == max(Di)) {
    Di <- rev(Di)
    counts <- rev(counts)
  }
  cpf <- cumsum(counts) / sum(counts)  #cumulative distribution
  #define the size class data
  x <- log2(Di)  #temp var in wentworth sizes
  dm <- 2^(x[-1] - diff(x)/2)  #median size of each size class
  F_i <- diff(cpf)  #fraction in each size class
  df <- data.frame(Di, cpf, dm = c(min(Di),dm), Fi = c(0, F_i))
  df <- df[order(-df$Di),]
  return(df)
}



#' Processing Wolman count data into grain size distribution
#'
#' This function uses standard data from a surface Wolman count to estimate
#' cumulative percent finer for each size, the fraction of the distribution in
#' a given size class, and the median size of each size class. Output is a data
#' frame containing the upper limit of each size class (Di), the cumulative
#' proportion finer than the size listed in Di (cpf), the median size for each
#' size class (Dm), and the fraction of the distribution in each class (Fi)
#' @param D50 an estimate of the median bed surface grain size (mm)
#' @param sp (optional) the standard deviation of the log-normal distribution (in phi units), default of 1 gives a D84 twice the D50
#' @param Map (optional) LOGICAL, when TRUE, the distribution is graphed
#' @export
sim_gsd <- function(D50, sp = 1, plot_gsd = FALSE) {
  # Simulates a surface grain size distribution by sampling from a log-normal
  # distribution in phi units, then binning into half-phi size classes.
  #
  # Parameters
  #
  # D50      : numeric, median grain diameter (mm).
  # sp       : numeric, spread of the distribution in phi units (analogous to
  #              sorting). Default 1.
  # plot_gsd : logical, if TRUE plots the cumulative grain size distribution.
  #              Default FALSE. Renamed from 'Map' for clarity.
  #
  # Returns
  #
  # A data frame with one row per size class, sorted from coarse to fine:
  #   Di  - upper bound of size class (mm)
  #   cpf - cumulative proportion finer (0-1)
  #   dm  - geometric mean grain size of size class (mm)
  #   Fi  - proportion of distribution in size class (0-1)

  # --- phi-space mean and truncation limits (+/- 2.5 phi) ---
  phi_mean  <- log2(D50)
  lim_lower <- 2^(phi_mean - 2.5 * sp)
  lim_upper <- 2^(phi_mean + 2.5 * sp)

  # --- simulate and truncate grain size sample ---
  grain_sizes <- 2^rnorm(1e5, mean = phi_mean, sd = sp)
  grain_sizes <- grain_sizes[grain_sizes >= lim_lower & grain_sizes <= lim_upper]

  # --- bin into half-phi size classes ---
  phi_breaks  <- 2^seq(phi_mean - 2.5 * sp, phi_mean + 2.5 * sp, by = 0.25 * sp)
  grain_dist  <- hist(grain_sizes, breaks = phi_breaks, plot = FALSE)

  # --- compute proportions and cumulative distribution ---
  p   <- grain_dist$counts / sum(grain_dist$counts)
  cdf <- cumsum(p)

  # --- assemble output, sorted coarse to fine ---
  results <- data.frame(
    Di  = grain_dist$breaks[-1],
    cpf = cdf,
    dm  = grain_dist$mids,
    Fi  = p
  )
  results <- results[order(-results$Di), ]

  # --- optional plot ---
  if (plot_gsd) {
    plot(results$Di, results$cpf,
         type = "o", log  = "x",
         xlim = c(lim_lower, lim_upper),
         xlab = "Grain size (mm)",
         ylab = "Proportion finer")
  }

  results
}



#' Creating a hydrograph based on the Standard SCS Dimensionless Unit Hydrograph
#'
#' This function uses the Soil Conservation Service Standard Unit Hydrograph to
#' simulate a simple hydrograph associated with a rainfall event. The hydrograph
#' depends on a parameter controlling the peak flow that will be reached and a
#' parameter controlling how quickly the hydrograph will rise and fall. The
#' shape of the hydrograph remains the same.
#'
#' The dimensionless unit hydrograph used by the SCS was developed by Victor
#' Mockus and was derived based on a large number of unit hydrographs from basins
#' that varied in characteristics such as size and geographic location. The unit
#' hydrographs were averaged and the final product was made dimensionless by
#' considering the ratios of Q/Qpeak (flow/peak flow) on the ordinate axis and
#' t/tpeak (time/time to peak) on the abscissa. This final, dimensionless unit
#' hydrograph, which is the result of averaging a large number of individual
#' dimensionless unit hydrographs, has a time-to-peak located at approximately
#' 20 percent of its time base and an inflection point at 1.7 times the time-to-peak.
#'
#' @param Qpeak the peak instantaneous discharge to be reached during the hydrograph (cumecs)
#' @param tpeak the lag time between the peak rainfall input and the peak runoff (hours)
#' @param dt a parameter that controls the time increment for data output (in hours)
#' @export
scs_hydrograph <- function(Qpeak, tpeak, dt){
  # values from https://www.nohrsc.noaa.gov/technology/gis/uhg_manual.html
  t_over_tp <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,
                 1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.5,5)

  Q_over_Qp <- c(0,0.03,0.1,0.19,0.31,0.47,0.66,0.82,0.93,0.99,1,0.99,0.93,0.86,0.78,
                 0.68,0.56,0.46,0.39,0.33,0.28,0.207,0.147,0.107,0.077,0.055,0.04,
                 0.029,0.021,0.015,0.011,0.005,0)

  #scale the values
  t <- tpeak * t_over_tp
  Q <- Qpeak * Q_over_Qp

  hydrograph <- approx(x = t, y = Q, xout = seq(0, max(t), dt), method = "linear")
  return(hydrograph)
}



#' Create a cross section based on measured trapezoid dimensions
#'
#'This function creates profile data including the distance along the profile
#'from the left bank.
#' @param w_top measured or estimated top width of the trapezoid (m)
#' @param w_bot measured or estimated bottom width of the trapezoid (m)
#' @param H_bank measures or estimated vertical height of the trapezoid (m)
#' @param z_fp elevation of the top of the trapezoid (i.e. the floodplain or terrace elevation, m)
#' @param dx a parameter controlling the spacing between exported data points on the section
#' @export
make_trapezoid <- function(w_top, w_bot, H_bank, Z_fp, dx = 0.1){
  # CREATE A CHANNEL CROSS SECTION with a trapexzoidal geometry
  dx <- 0.1
  dw <- (w_top - w_bot) / 2
  dist <- c(0, dw, dw + w_bot, w_top)
  stad <- c(0, H_bank, H_bank, 0)
  z <- Z_fp - stad

  xs_data <- approx(x = dist,
                    y = z,
                    xout = seq(min(dist), max(dist), dx))
  xs_data <- as.data.frame(xs_data)
  colnames(xs_data) <- c("dist", "Z")
  return(xs_data)
}



#' the GBEM function: gravel-bed river bank erosion model
#'
#' This function implements the gravel bed bank erosion model (GBEM) that is the
#' core function used by STOCHASIM and the BGC stochastic bank erosion model.  It
#' limits the amount of bank erosion that can occur in a given timestep based on
#' the bedload transport capacity of the stream proximate to the channel banks.
#' The output from this function is a vector containing: (1) the predicted widening
#' based only on a stability criterion; (2) the predicted widening based on
#' the stability criterion, constrained by the bedload transport capacity; (4)
#' the volume of bedload that could be transported during the time step t in
#' the near-bank region; and (4) the critical depth for bank erosion
#'
#' @param Q discharge carried by the stream (m3/s)
#' @param t time for which Q acts on the stream channel (hrs)
#' @param n Manning's n value for the main channel
#' @param D84 84th percentile of the surface grain size distribution (m)
#' @param D50 50th percentile of the grain size distribution (m)
#' @param W water surface width at the beginning of time interval t
#' @param S energy gradient of the stream channel (m/m)
#' @param H effective rooting depth for vegetation; grassy banks, no trees / shrubs H = 0.35; 1 to 5 percent tree / shrub cover H = 0.50; 5 to 50 percent tree / shrub cover H =  0.90; more than 50 percent tree / shrub cover H = 1.1
#' @param fun choice of transport function. Default is Eaton and Church ("ec")
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @param rho_s sediment density. Defaults to 2650 kg/m3
#' @export
gbem <- function(Q, t, n, D84, D50, W, S, H = 0, fun = "ec",
                 g = 9.81, rho_s = 2650, rho = 1000) {
  # Estimates lateral channel widening and updated flow depth for a single
  # flood event using the gravel-bed erosion model (GBEM) framework.
  # Widening is limited by the lesser of the geotechnically predicted widening
  # and the volume of bedload transport in the bank zone.
  #
  # Parameters
  #
  # Q     : numeric, peak discharge (m^3/s)
  # t     : numeric, flood duration (hours)
  # n     : numeric, Manning's roughness coefficient
  # D84   : numeric, 84th percentile grain diameter (m). Used for bank
  #           stability threshold.
  # D50   : numeric, median grain diameter (m). Used for bedload transport.
  # W     : numeric, current channel width (m)
  # S     : numeric, channel slope (dimensionless)
  # H     : numeric, bank height (m). Default 0.
  # fun   : character, bedload transport equation to use. Default "ec"
  #           (Eaton & Church 2011). Passed to find_q_b().
  # g     : numeric, gravitational acceleration (m/s^2). Default 9.81.
  # rho_s : numeric, sediment density (kg/m^3). Default 2650.
  # rho   : numeric, fluid density (kg/m^3). Default 1000.
  #
  # Returns
  #
  # A named numeric vector with four elements:
  #   dw      - predicted channel widening (m); 0 if channel is stable
  #   depth   - updated flow depth after widening (m)
  #   v_b     - volumetric bedload transport during event (m^3/m width)
  #   d_crit  - critical flow depth for bank erosion threshold (m)

  # --- constants ---
  shields_c84    <- 0.02
  sec_per_hour   <- 3600
  travel_angle   <- 30   # fahrboschung angle for small failures in sand/gravel (degrees)

  # --- step 1: critical threshold for channel widening ---
  t_c84   <- shields_c84 * g * (rho_s - rho) * D84
  d_crit  <- find_d_crit(H, t_c84, S)
  v_crit  <- d_crit^(2/3) * S^0.5 / n

  # --- step 2: flow depth at current width and bedload transport volume ---
  d    <- ((n * Q) / (W * S^0.5))^(3/5)
  q_b  <- find_q_b(d, S, D50, D84, fun, g, rho, rho_s)
  v_b  <- q_b * t * sec_per_hour

  # --- step 3: widening, limited by transport volume in bank zone ---
  # if stable, no widening; otherwise widen to stable width, capped by
  # available transport volume: dw = v_b / tan(travel_angle), derived from
  # the bank zone width being proportional to bank height via the fahrboschung
  # angle (see Hassan et al. 1992)
  dw <- if (d < d_crit) {
    0
  } else {
    W_stable <- Q / (d_crit * v_crit)
    min(W_stable - W, v_b / tan(travel_angle * pi / 180))
  }

  # --- step 4: updated width and flow depth ---
  W_new <- W + dw
  depth <- ((n * Q) / (W_new * S^0.5))^(3/5)

  c(dw = dw, depth = depth, v_b = v_b, d_crit = d_crit)
}


#' Estimate the cumulative bank erosion during a flood hydrograph in an incized channel
#'
#' This function implements the GBEM model for estimating bank erosion using
#' a simple adjustment parameter to adapt the standard model for streams where
#' the channel is incised, and confined by a terrace. The function outputs a
#' data frame containing variables for time since the start of the hydrograph
#' (time), discharge values for each time (Q), the predicted bank erosion
#' distance for the time step (erosion), the volume of sediment transport that
#' the stream can convey during the timestep (vb), the critical depth for bank
#' erosion (dc), channel width at the end of each timestep (wi)
#'
#' @param Q a list of discharge values at regular time intervals defining the flood hydrograph (m3/s)
#' @param t time interval between discharge measurements (hrs)
#' @param n Manning's n value for the main channel
#' @param D84 84th percentile of the surface grain size distribution (m)
#' @param D50 50th percentile of the grain size distribution (m)
#' @param W water surface width at the beginning of time interval t
#' @param S energy gradient of the stream channel (m/m)
#' @param H effective rooting depth for vegetation grassy banks, no trees / shrubs H = 0.35;  1 to 5percent tree / shrub cover H = 0.50; 5 to 50percent tree / shrub cover H =  0.90; more than 50percent tree / shrub cover H = 1.1
#' @param Z height of the terrace above the channel bottom (m)
#' @param fun choice of transport function. Default is Eaton and Church ("ec")
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @param rho_s sediment density. Defaults to 2650 kg/m3
#' @export
bank_erode_trc <- function(Q, t, n, D84, D50, W, S, H = 0, Z, fun = "ec", g = 9.81, rho_s = 2650, rho = 1000){

  event <- cbind(seq_along(Q), Q)
  extra <- matrix(data = NA,
                  nrow = nrow(event),
                  ncol = 5)
  event <- cbind(event, extra)
  wi <- W
  for( i in seq(1, nrow(event), 1) ){
    event[i, 3:6] <- gbem(event[i,2],t , n, D84, D50, wi, S, H, fun)
    event[i,3] <- event[i, 3] * event[i,4] / Z #correct based on ratio of flow depth and terrace height
    event[i, 5] <- event[i,5] * wi #scale up to total volumes (m3)
    wi <- wi +  event[i,3] #widen the channel
    event[i, 7] <- wi
  }
  event <- as.data.frame(event[,c(1,2,4,3,5,6,7)])
  colnames(event) <- c('time', 'Q', 'depth', 'erosion', 'Vb', 'dc', 'wi')
  return(event)
}



#' Estimate the cumulative bank erosion  during a flood hydrograph for a floodplain channel
#'
#' This function implements the GBEM model for estimating bank erosion using
#' the standard model for streams where the floodplain is assumed to be equal
#' to the critical depth based on the bed texture and slope. The function outputs a
#' data frame containing variables for time since the start of the hydrograph
#' (time), discharge values for each time (Q), the predicted bank erosion
#' distance for the time step (erosion), the volume of sediment transport that
#' the stream can convey during the timestep (vb), the critical depth for bank
#' erosion (dc), channel width at the end of each timestep (wi)
#'
#' @param Q a list of discharge values at regular time intervals defining the flood hydrograph (m3/s)
#' @param t time interval between discharge measurements (hrs)
#' @param n Manning's n value for the main channel
#' @param D84 84th percentile of the surface grain size distribution (m)
#' @param D50 50th percentile of the grain size distribution (m)
#' @param W water surface width at the beginning of time interval t
#' @param S energy gradient of the stream channel (m/m)
#' @param H effective rooting depth for vegetation
#' @param Z height of the terrace above the channel bottom (m)
#' @param fun choice of transport function. Default is Eaton and Church ("ec")
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @param rho_s sediment density. Defaults to 2650 kg/m3
#' @export
bank_erode_fp <- function(Q, t, n, D84, D50, W, S, H = 0, fun = "ec"){

  event <- cbind(seq_along(Q), Q)
  extra <- matrix(data = NA,
                  nrow = nrow(event),
                  ncol = 5)
  event <- cbind(event, extra)
  wi <- W
  for( i in seq(1, nrow(event), 1) ){
    event[i, 3:6] <- gbem(event[i,2],t , n, D84, D50, wi, S, H, fun)
    event[i, 5] <- event[i,5] * wi #scale up to total volumes (m3)
    wi <- wi +  event[i,3] #widen the channel
    event[i, 7] <- wi
  }
  event <- as.data.frame(event[,c(1,2,4,3,5,6,7)])
  colnames(event) <- c('time', 'Q', 'depth', 'erosion', 'Vb', 'dc', 'wi')
  return(event)
}



#' find the critical depth at which bank erosion will begin
#'
#' this is a utility function used by GBEM to determine the typical floodplain
#' thickness above the channel bed.
#'
#' @param H effective rooting depth for vegetation
#' @param t_c84 stress required to entrain the 84th percentile of the bed surface
#' @param S energy gradient of the stream channel (m/m)
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @export
find_d_crit <- function(H, t_c84, S, g = 9.81, rho = 1000, tol = 0.0001){

  d_crit <- t_c84 / (g * rho * S)  #max depth that can be maintained (no veg)

  if(H > 0){   #if mu > 0, adjust threshold and depth
    #bounds <- c(1, est_mu(H, d_crit)) * d_crit
    bounds <- c(1, 4) * d_crit
    d_test <- mean(bounds)
    d_target <- t_c84 * est_mu(H, d_test) / (g * rho * S)
    converg <- (d_test - d_target) / d_target
    while(abs(converg) > tol){
      if (converg > 0) {
        bounds[2] <- d_test #do this if depth.test > depth.target
      }else if (converg < 0) {
        bounds[1] <- d_test #do this if depth.target > depth.test
      }
      d_test <- mean(bounds)
      d_target <- t_c84 * est_mu(H, d_test) / (g * rho * S)
      converg <- (d_test - d_target) / d_target
    }
    d_crit <- d_test
  }
  return(d_crit)
}



#' Calculating volumetric unit transport rate for a specified function
#'
#' This function takes the hydraulic data required to calculate the dimensionless
#' shear stress and uses a specified bedload transport function to calculate
#' the volumetric transport rate per unit width
#'
#' @param R hydraulic radius for the flow (m)
#' @param S gradient of the water surface slope (m/m)
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param D84 84th percentile of the surface grain size distribution (m)
#' @param D50 50th percentile of the grain size distribution (m)
#' @param fun choice of transport function. Default is Eaton and Church ("ec")
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @param rho_s sediment density. Defaults to 2650 kg/m3
#' @export
find_q_b <- function(R, S, D50, D84, fun = "ec", g = 9.81, rho = 1000, rho_s = 2650){
  #fun can be "ec", "eb", "wp", "yln", "vr"
  Gs <- (rho_s - rho) / rho
  if(fun == "vr"){
    qb_star <- vr(R, S, D50, rho = rho, rho_s = rho_s, g = g)
    qb <- qb_conversion(qb_star, D50, g = g, Gs = Gs)
  }
  if(fun == "yln"){
    qb_star <- yln(R, S, D50, rho = rho, rho_s = rho_s, g = g)
    qb <- qb_conversion(qb_star, D50, g = g, Gs = Gs)
  }
  if(fun == "wp"){
    qb_star <- wp(R, S, D50, rho = rho, rho_s = rho_s, g = g)
    qb <-  qb_conversion(qb_star, D50, g = g, Gs = Gs)
  }
  if(fun == "mpm"){
    qb_star <- mpm(R, S, D50, rho = rho, rho_s = rho_s, g = g)
    qb <-  qb_conversion(qb_star, D50, g = g, Gs = Gs)
  }
  if(fun == "eb"){
    qb_star <- eb(R, S, D50, rho = rho, rho_s = rho_s, g = g)
    qb <- qb_conversion(qb_star, D50, g = g, Gs = Gs)
  }
  if(fun == "ec"){
    if(R > 0 & S >0){
      qb_star <- ec(R, S, D50, D84, rho = rho, rho_s = rho_s, g = g)
    } else {
      qb_star = 0
    }
    qb <- qb_conversion(qb_star, D50, g = g, Gs = Gs)
  }
  if(fun == "rk"){
    qb_star <- rk(R, S, D50, D84, rho = rho, rho_s = rho_s, g = g)
    qb <- qb_conversion(qb_star, D84, g = g, Gs = Gs)
  }
  return(qb)
}



#' estimate relative bank strength
#'
#' this is a utility function used by GBEM to determine the relative bank strength
#' given the channel depth and rooting depth. It returns a variable that is the
#' bank strength relative to the bed strengths.
#'
#' @param H effective rooting depth for vegetation
#' @param d channel depths
#' @export
est_mu <- function(H, d){
  if(H / d > 0.94){
    mu <- 4  #set an upper threshold of 4
  } else {
    a = 0.85
    b = 0.87
    mu = 1 / (1 - a*(H/d))^b  #from Eaton and Millar 2017
  }
  return(mu)
}




#' Estimating critical stream power for bed entrainment
#'
#' Using the equation published by Ferguson (2005) or Petit et al. (2005), this
#' function estimates the critical stream power as a function of grain size,
#' and flow depth. These empirical equations are valid for D50 ranging from
#' 0.020 to 0.150 m
#' @param D50 median bed sediment grain diameter (m)
#' @param d flow depth (m)
#' @export
crit_power <- function(D50, d, use_ferguson = FALSE) {
  # Estimates the critical unit stream power required to initiate sediment
  # transport, using one of two formulations.
  #
  # Parameters
  #
  # D50          : numeric, median grain diameter (m)
  # d            : numeric, flow depth (m). Only used when use_ferguson = TRUE.
  # use_ferguson : logical, if TRUE uses the Ferguson (2005) depth-dependent
  #                  formulation; if FALSE uses the empirical regression of
  #                  Hassan et al. (1992). Default FALSE.
  #
  # Returns
  #
  # numeric, critical unit stream power (W/m^2).

  if (use_ferguson) {
    2860 * D50^1.5 * log10(12 * d / D50)
  } else {
    0.130 * (D50 * 1000)^1.438
  }
}


#' Calculating the cross section average unit stream power
#'
#' this function calculates the unit stream power acting on the channel
#' boundary for a given cross section.
#' @param q specific discharge (m3/s/m)
#' @param S gradient of the water surface slope (m/m)
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @export
unit_power <- function(q, S, g = 9.81, rho = 1000){
  unit_p <- q * rho * g * S
  return(unit_p)
}



#' Calculating the step length for sediment transport of coarse bed material
#'
#' this function calculates the average distance for bedload transport during
#' a flow of a given duration
#' @param D50 median bed sediment grain diameter (m)
#' @param q a list of specific discharge values for a hydrograph (m3/s/m)
#' @param d flow depth (m)
#' @param S gradient of the water surface slope (m/m)
#' @param dt time interval between discharge measurements (hrs)
#' @param flashy a logical variable (true if the hydrograph is short and flashy)
#' @param use_ferguson logical variable for selecting the ferguson equation for calculating the critical stream power
#' @param g the acceleration of gravity. defaults to 9.81 m/s2
#' @param rho fluid density. Defaults to 1000 kg/m3
#' @export
step_length <- function(D50, q, d, S, dt, flashy = FALSE, use_ferguson = TRUE,
                        g = 9.81, rho = 1000) {
  # Estimates the virtual travel distance of a representative grain over a
  # flood event, using the relationships of Hassan et al. (1992).
  #
  # Hassan, M.A., M. Church, and P.J. Ashworth. 1992. Virtual rate and mean
  # distance of travel of individual clasts in gravel-bed channels. Earth
  # Surface Processes and Landforms 17(6): 617-627.
  #
  # Parameters
  #
  # D50          : numeric, median grain diameter (m)
  # q            : numeric vector, unit discharge time series (m^2/s)
  # d            : numeric, flow depth (m), used in crit_power()
  # S            : numeric, channel slope (dimensionless)
  # dt           : numeric vector, time step durations (s), same length as q
  # flashy       : logical, if TRUE uses the peak unit stream power form
  #                  (suitable for short, flashy events); if FALSE integrates
  #                  virtual velocity over the full hydrograph. Default FALSE.
  # use_ferguson : logical, passed to crit_power() to select the resistance
  #                  equation used for critical stream power. Default TRUE.
  # g            : numeric, gravitational acceleration (m/s^2). Default 9.81.
  # rho          : numeric, fluid density (kg/m^3). Default 1000.
  #
  # Returns
  #
  # numeric, estimated virtual travel distance (m) over the flood event.

  omega_crit <- crit_power(D50, d, use_ferguson)

  if (flashy) {
    # peak unit stream power form (Hassan et al. 1992, eq. for flashy events)
    omega_peak <- unit_power(max(q, na.rm = TRUE), S, g, rho)
    0.0283 * (omega_peak - omega_crit)^1.44
  } else {
    # integrate virtual velocity over the hydrograph
    omega    <- unit_power(q, S, g, rho)
    virt_v   <- 0.00188 * (omega - omega_crit)^1.62
    sum(virt_v * dt, na.rm = TRUE)
  }
}



#' Estimate flood quantiles for annual daily maximum flow
#'
#' This function loads data from WSC or USGS stations and fits a generalized
#' extreme value (GEV) distribution to the data using maximum likelihood
#' estimation.
#' @param station a string containing the station name (e.g. '08GA022' or '14048000')
#' @param rp a vector containing the return periods for which to estimate flood quantiles
#' @param qual a value from 0 to 1 specifying how complete an annual dataset should be in order to estimate an annual maximum discharge value for that year
#' @export
quick_ffa_ann <- function(station, rp = c(2, 5, 10, 20, 50, 100, 200), qual = 0.6){
  #station <- '14048000'
  if(nchar(station) == 7) {
    #download the data from the hydat database
    q_df <- RivRetrieve::canada(station,
                   variable = 'discharge')
  }else {
    #download the data from USGS
    q_df <- RivRetrieve::usa(station,
                variable = 'discharge')
  }
  q_df$year <- as.numeric(format(q_df$Date, format = '%Y'))
  q_df$day <- as.numeric(format(q_df$Date, format = '%j'))
  #make a data frame for calculating the average/max/min values for a year
  ann_df <- data.frame(years = unique(q_df$year))
  ann_df$Q_max <- NA #max annual discharge
  ann_df$qual <- NA

  #collect the mean/ max/min annual discharges
  for (j in seq_along(ann_df$years)) {
    filt <- q_df$year == ann_df$years[j]
    yearly_q <- q_df[filt,]
    ann_df$Q_max[j] <- max(yearly_q$Q, na.rm = TRUE)
    ann_df$qual[j] <- sum(filt) / 366
  }

  #get rid of bad data
  ann_df <- ann_df[which(ann_df$qual > qual),]

  #fit GEVs to data
  q_max <- sort(ann_df$Q_max, decreasing = TRUE)
  ann_n <- length(q_max)
  ann_rp <- (1 + ann_n) / seq(1, ann_n)
  ann_yt <- -log(-log( 1 -  1/ ann_rp))
  ann_gev <- extRemes::fevd(q_max,
                  type = 'GEV')

  # use fits to calculate flows for select return periods
  #rp <- c(2, 5, 10, 20, 25, 50, 100, 200, 500)
  mod_yt <-  -log(-log( 1 -  1/ rp))
  ann_mod_gev <- extRemes::qevd(p = 1 - 1 / rp,
                      loc = ann_gev$results$par[1],
                      scale = ann_gev$results$par[2],
                      shape = ann_gev$results$par[3])

  # plot the data
  plot(ann_yt,
       q_max,
       type = "p",
       pch = 20,
       col = rgb(0,0,1),
       ylim = c(0, 1.25*max(q_max)),
       xlim = c(-2, 6),
       main = paste('Flood frequency analysis for Station',station, sep = " "),
       ylab = expression(Discharge~(m^3/s)), xlab = "Return Period (yrs)",
       xaxt = 'n')
  # add the GEV fit
  lines(mod_yt, ann_mod_gev, col = 'blue')

  # add the return period labels to the axis
  axis(side = 1, at = mod_yt, labels = rp)

  fit_data <- data.frame(return_period = rp,
                         estimate = signif(ann_mod_gev, digits = 3)
  )
  return(fit_data)
}



#' Estimate flood quantiles for dual mechanism daily maximum flows
#'
#' This function loads data from WSC or USGS stations and fits a generalized
#' extreme value (GEV) distribution to the data using maximum likelihood
#' estimation.
#' @param station a string containing the station name (e.g. '08GA022' or '14048000')
#' @param start a string specifying the beginning of the melt window in MM-DD format
#' @param rp a vector containing the return periods for which to estimate flood quantiles
#' @param qual a value from 0 to 1 specifying how complete an annual dataset should be in order to estimate an annual maximum discharge value for that year
#' @export
quick_ffa_dual <- function(station, start = '04-01', end = '07-01', rp = c(2, 5, 10, 20, 50, 100, 200), qual = 0.6){
  #station <- '14048000'
  if(nchar(station) == 7) {
    #download the data from the hydat database
    q_df <- RivRetrieve::canada(station,
                                variable = 'discharge')
  }else {
    #download the data from USGS
    q_df <- RivRetrieve::usa(station,
                             variable = 'discharge')
  }
  q_df$year <- as.numeric(format(q_df$Date, format = '%Y'))
  q_df$day <- as.numeric(format(q_df$Date, format = '%j'))

  #make a data frame for calculating the average/max/min values for a year
  ann_df <- data.frame(years = unique(q_df$year))
  ann_df$Q_max <- NA #max annual discharge
  ann_df$qual <- NA

  #collect the mean/ max/min annual discharges
  for (j in seq_along(ann_df$years)) {
    filt <- q_df$year == ann_df$years[j]
    yearly_q <- q_df[filt,]
    ann_df$Q_max[j] <- max(yearly_q$Q, na.rm = TRUE)
    ann_df$qual[j] <- sum(filt) / 366
  }

  #get rid of bad data
  ann_df <- ann_df[which(ann_df$qual > qual),]
  ann_n <- length(ann_df$Q_max)
  ann_rp <- (1 + ann_n) / seq(1, ann_n)
  ann_yt <- -log(-log( 1 -  1/ ann_rp))

  begin_melt <- as.Date(start, format = '%m-%d')
  end_melt <- as.Date(end, format = '%m-%d')
  melt_window <- c(as.numeric(format(begin_melt, format = '%j')),
                     as.numeric(format(end_melt, format = '%j')))

    #collect the snowmelt and rain annual discharges
    for (j in seq_along(ann_df$years)) {
      filt <- q_df$year == ann_df$years[j]
      yearly_q <- q_df[filt,]
      season <- yearly_q$day > melt_window[1] & yearly_q$day < melt_window[2]
      ann_df$Qmelt[j] <- max(yearly_q$Q[season], na.rm = TRUE)
      ann_df$Qrain[j] <- max(yearly_q$Q[!season], na.rm = TRUE)
    }

    #fit the rain distribution
    q_rain = sort(ann_df$Qrain, decreasing = TRUE)
    rain_gev <- extRemes::fevd(q_rain,
                     type = 'GEV')
    rain_mod_gev <- extRemes::qevd(p = 1 - 1 / rp,
                         loc = rain_gev$results$par[1],
                         scale = rain_gev$results$par[2],
                         shape = rain_gev$results$par[3])

    #fit the melt distribution
    q_melt = sort(ann_df$Qmelt, decreasing = TRUE)
    melt_gev <- extRemes::fevd(q_melt,
                     type = 'GEV')
    melt_mod_gev <- extRemes::qevd(p = 1 - 1 / rp,
                         loc = melt_gev$results$par[1],
                         scale = melt_gev$results$par[2],
                         shape = melt_gev$results$par[3])

    #create an empirical distribution using these two GEVs
    n <- 90000 #number of values needed to simulate empirical dist.
    seed <- 17
    empirical_sim <- pmax(extRemes::revd(n,
                               loc = melt_gev$results$par[1],
                               scale = melt_gev$results$par[2],
                               shape = melt_gev$results$par[3],
                               type = 'GEV'),
                          extRemes::revd(n,
                               loc = rain_gev$results$par[1],
                               scale = rain_gev$results$par[2],
                               shape = rain_gev$results$par[3],
                               type = 'GEV'))
    empirical_sim <- sort(empirical_sim, decreasing = TRUE)
    emp_rp <- (1 + n) / seq(1, n)
    emp_yt <- -log(-log( 1 -  1/ emp_rp))

    mod_yt <-  -log(-log( 1 -  1/ rp))
    #plot the data
    plot(ann_yt,
         q_rain,
         type = "p",
         pch = 19,
         col = rgb(0,0,1),
         ylim = c(0, 1.5*max(ann_df$Q_max)),
         xlim = c(-2, 6),
         main = paste('Flood frequency analysis for station',station, sep = " "),
         ylab = expression(Discharge~(m^3/s)), xlab = "Return Period (yrs)",
         xaxt = 'n')
    points(ann_yt,
           q_melt,
           pch = 19,
           col = 'red'
    )
    # add the GEV fit based on snowmelt data
    lines(mod_yt, melt_mod_gev, col = 'red', lty = 2)

    # add the GEV fit based on rain domin data
    lines(mod_yt, rain_mod_gev, col = 'blue', lty = 2)

    # add the empirical fit based on rain and snow
    lines(emp_yt, empirical_sim, col = "darkorange")

    # add the return period labels to the axis
    axis(side = 1, at = mod_yt, labels = rp)

    legend(
      "topleft",
      inset = c(0.01, 0.015),
      legend = c('snowmelt', 'rain-dominated',  'dual-mech. FFA'),
      pch = c(19, 19, NA, NA),
      lty = c(2, 2, 1, 1),
      ncol = 1,
      cex = 1,
      col = c("red", 'blue', 'darkorange')
    )

    dual_Q <- approx(x = emp_rp,
                     y = empirical_sim,
                     xout = rp)
    output <- data.frame(return_period = rp, estimate = signif(dual_Q$y, digits = 3))
    return(output)
}


#' Display timing for annual daily maximum and minimum flows
#'
#' This function loads data from WSC or USGS stations and extracts the annual
#' daily maximum flow and minimum flow, and then plots it against day of the year
#' so that peak flow and low flow timing can be quickly assessed
#'
#' @param station a string containing the station name (e.g. '08GA022' or '14048000')
#' @param qual a value from 0 to 1 specifying how complete an annual dataset should be in order to estimate an annual maximum discharge value for that year
#' @export
max_min_timing <- function(station, qual = 0.6){
  #station <- '14048000'
  if(nchar(station) == 7) {
    #download the data from the hydat database
    q_df <- RivRetrieve::canada(station,
                                variable = 'discharge')
  }else {
    #download the data from USGS
    q_df <- RivRetrieve::usa(station,
                             variable = 'discharge')
  }
  q_df$year <- as.numeric(format(q_df$Date, format = '%Y'))
  q_df$day <- as.numeric(format(q_df$Date, format = '%j'))
  #make a data frame for calculating the average/max/min values for a year
  ann_df <- data.frame(years = unique(q_df$year))
  ann_df$Q_max <- NA #max annual discharge
  ann_df$qual <- NA

  #collect the mean/ max/min annual discharges
  for (j in seq_along(ann_df$years)) {
    filt <- q_df$year == ann_df$years[j]
    yearly_q <- q_df[filt,]
    ann_df$Q_max[j] <- max(yearly_q$Q, na.rm = TRUE)
    ann_df$Q_min[j] <- min(yearly_q$Q, na.rm = TRUE)
    ann_df$qual[j] <- sum(filt) / 366
  }

  #get rid of bad data
  ann_df <- ann_df[which(ann_df$qual > qual),]

  #find the dates for peaks
  maxdate <- array(data = NA, dim = nrow(ann_df))
  mindate <- array(data = NA, dim = nrow(ann_df))
  for(j in seq_along(ann_df$years)){
    maxfilt <- which( q_df$year == ann_df$years[j] & q_df$Q == ann_df$Q_max[j])[1]
    maxdate[j] <- as.character(q_df$Date[maxfilt])
    minfilt <- which( q_df$year == ann_df$years[j] & q_df$Q == ann_df$Q_min[j])[1]
    mindate[j] <- as.character(q_df$Date[minfilt])
  }
  ann_df$maxdate <- as.Date(maxdate)
  ann_df$mindate <- as.Date(mindate)

  months <- as.Date(paste('1900-', seq(1,12),'-01', sep = ""))
  months <-c(months, "1900-12-31")
  month_labels <-  c('Jan-01', 'Feb-01', 'Mar-01', 'Apr-01', 'May-01', 'Jun-01', 'Jul-01', 'Aug-01', 'Sep-01', 'Oct-01', 'Nov-01', 'Dec-01', 'Dec-31')

  par(mar = c(5,5,1,1))
  plot(as.numeric(format(ann_df$maxdate, format = '%j')),
       ann_df$Q_max,
       type = 'p',
       pch = 19,
       cex = 1,
       col = 'blue',
       xlab = "Date",
       ylab = expression(Discharge~(m^3/s)),
       xlim = c(0,366),
       ylim = c(0, max(ann_df$Q_max)),
       xaxt = 'n')
  points(as.numeric(format(ann_df$mindate, format = '%j')),
         ann_df$Q_min,
         type = 'p',
         pch = 19,
         cex = 1,
         col = 'red',)
  axis(side = 1,
       cex.axis = 1,
       at = as.numeric(format(months, format = "%j")),
       labels = month_labels)
  abline(v = as.numeric(format(months, format = "%j")),
         lwd = 0.5,
         col = 'grey')
  legend(
    "topright",
    inset = c(0.01, 0.02),
    legend = c('low flows', 'high flows'),
    pch = 19,
    ncol = 2,
    cex = 1,
    col = c("red", 'blue')
  )
  return(ann_df)
}



#' Display flow regime
#'
#' This function loads data from WSC or USGS stations and extracts the mean
#' daily flow for every day of the year, then generates a plot showing the
#' flow regime.
#'
#' @param station a string containing the station name (e.g. '08GA022' or '14048000')
#' @param window (optional) the size of the window to use if you wish to smooth the flow regime by calculating a rolling average
#' @export
flow_regime <- function(station, window = 1){
  #station <- '14048000'
  if(nchar(station) == 7) {
    #download the data from the hydat database
    q_df <- RivRetrieve::canada(station,
                                variable = 'discharge')
  }else {
    #download the data from USGS
    q_df <- RivRetrieve::usa(station,
                             variable = 'discharge')
  }
  q_df$year <- as.numeric(format(q_df$Date, format = '%Y'))
  q_df$day <- as.numeric(format(q_df$Date, format = '%j'))

  #calculate the rolling mean from the daily data for the entire period of record
  q_df$rolling <- data.table::frollmean(q_df$Q,
                            window,
                            align = 'center'
  )
  q_df$rolling[is.infinite(q_df$rolling)] <- NA

  med_Q <- aggregate(q_df$rolling,
                     by = list(q_df$day),
                     na.rm = TRUE,
                     FUN = quantile,
                     probs = c(0.5))

  range <- aggregate(q_df$rolling,
                     by = list(q_df$day),
                     na.rm = TRUE,
                     FUN = quantile,
                     probs = c(0.25, 0.75))

  par(mar = c(5,5,3,1))
  #the median daily flow based on all years of record
  plot(med_Q,
       col = "darkorange",
       pch = 19,
       main = paste('median daily flows for all years of record for Station',
                    station),
       xlab = "Date",
       ylab = expression(Discharge~(m^3/s)),
       ylim = c(0, max(range$x, na.rm = T)),
       cex = 0.75,
       xaxt = "n")
  lines(range$x[,1])
  lines(range$x[,2])
  months <- as.Date(paste('1900-', seq(1,12),'-01', sep = ""))
  months <-c(months, "1900-12-31")
  month_labels <-  c('Jan-01', 'Feb-01', 'Mar-01', 'Apr-01', 'May-01', 'Jun-01', 'Jul-01', 'Aug-01', 'Sep-01', 'Oct-01', 'Nov-01', 'Dec-01', 'Dec-31')
  axis(side = 1,
       cex.axis = 1,
       at = as.numeric(format(months, format = "%j")),
       labels = month_labels)
  abline(v = as.numeric(format(months, format = "%j")),
         lwd = 0.5,
         col = 'grey')
  legend('topright',
         inset = c(0.01, 0.02),
         legend = c(paste(window,'-day median flow', sep = ''), 'inter-quartile range'),
         lty = c(0,1),
         pch = c(19,NA),
         cex = 1.0,
         col = c("darkorange", 'black'))
}


#' Bedload transport for steep channels using Smart & Jaeggi (1983) / Smart (1984)
#'
#' Computes the dimensionless bedload transport rate for steep channels using
#' the Smart (1984) uniform-sediment equation, calibrated for slopes between
#' 0.5% and 20% and grain sizes greater than 0.4 mm. The critical Shields
#' parameter is estimated using the Soulsby & Whitehouse (1997) relationship,
#' making it suitable for sand as well as gravel. The Chezy flow resistance
#' factor is computed internally using the Ferguson (2007) variable power
#' equation, consistent with \code{u_ferg()}, and non-dimensionalised for
#' compatibility with \code{qb_conversion()}.
#'
#' @param R hydraulic radius for the flow (m)
#' @param S channel slope (m/m)
#' @param d median bed sediment grain diameter, D50 (m)
#' @param d84 84th percentile grain diameter of bed sediment (m). Used in the
#'   Ferguson resistance equation. Defaults to 2 * d.
#' @param t_crit critical Shields parameter (dimensionless). Defaults to the
#'   value returned by \code{t_crit_soul(d, rho, rho_s, g, kv)}, which accounts
#'   for grain-size-dependent entrainment thresholds.
#' @param a1 Ferguson resistance coefficient for the smooth/deep flow limit.
#'   Defaults to 6.5.
#' @param a2 Ferguson resistance coefficient for the rough/shallow flow limit.
#'   Defaults to 2.5.
#' @param g acceleration due to gravity (m/s^2). Defaults to 9.81.
#' @param rho fluid density (kg/m^3). Defaults to 1000.
#' @param rho_s sediment density (kg/m^3). Defaults to 2650.
#' @param kv kinematic viscosity of the fluid (m^2/s). Defaults to 1e-6.
#'   Passed to \code{t_crit_soul()} when using the default critical threshold.
#' @return A numeric value giving the dimensionless bedload transport rate
#'   (q_star) in the Einstein normalisation, compatible with
#'   \code{qb_conversion()}. Returns 0 when the Shields stress is at or below
#'   the critical threshold.
#' @references
#'   Smart, G.M. (1984). Sediment transport formula for steep channels.
#'   \emph{Journal of Hydraulic Engineering}, 110(3), 267--276.
#'
#'   Smart, G.M. and Jaeggi, M.N.R. (1983). Sediment transport on steep slopes.
#'   \emph{Mitteilungen der Versuchsanstalt fur Wasserbau, Hydrologie und
#'   Glaziologie}, ETH Zurich, No. 64.
#'
#'   Soulsby, R.L. and Whitehouse, R.J.S. (1997). Threshold of sediment motion
#'   in coastal environments. \emph{Proceedings of the Pacific Coasts and Ports
#'   Conference}, Christchurch, New Zealand, 149--154.
#' @seealso \code{\link{u_ferg}}, \code{\link{t_crit_soul}},
#'   \code{\link{tau_1D}}, \code{\link{t_star}}, \code{\link{qb_conversion}}
#' @export
sj <- function(R, S, d, d84 = 2 * d,
               t_crit = t_crit_soul(d, rho, rho_s, g, kv),
               a1 = 6.5, a2 = 2.5,
               g = 9.81, rho = 1000, rho_s = 2650, kv = 1e-6) {

  tau <- tau_1D(R, S, g, rho)
  t_s <- t_star(tau, d, g, rho, rho_s)

  if (t_s <= t_crit) return(0)

  Gs  <- (rho_s - rho) / rho

  # dimensional Chezy coefficient (m^0.5/s), then non-dimensionalised
  # by sqrt(Gs * g * d) for compatibility with the Einstein q_star convention
  U   <- u_ferg(R, S, d84, a1, a2, g)
  C   <- (U / sqrt(R * S)) / sqrt(Gs * g * d)

  4.2 * S^0.6 * C * (t_s - t_crit) * sqrt(t_s)
}

#' Total sediment transport using Engelund & Hansen (1967)
#'
#' Computes the dimensionless total sediment transport rate (bedload plus
#' suspended load) using the Engelund & Hansen (1967) equation. Can be used
#' as a drop-in replacement for bedload-only functions in this package,
#' but note that it will overestimate bedload transport as it includes the
#' suspended load component. Most appropriate for sand-bed channels with
#' substantial suspended load (d50 between 0.19 and 0.93 mm).
#'
#' @param R hydraulic radius for the flow (m)
#' @param S channel slope (m/m)
#' @param d median bed sediment grain diameter, D50 (m)
#' @param d84 84th percentile grain diameter of bed sediment (m). Used in the
#'   Ferguson resistance equation to compute the friction factor. Defaults to
#'   2 * d.
#' @param a1 Ferguson resistance coefficient for the smooth/deep flow limit.
#'   Defaults to 6.5.
#' @param a2 Ferguson resistance coefficient for the rough/shallow flow limit.
#'   Defaults to 2.5.
#' @param g acceleration due to gravity (m/s^2). Defaults to 9.81.
#' @param rho fluid density (kg/m^3). Defaults to 1000.
#' @param rho_s sediment density (kg/m^3). Defaults to 2650.
#' @return A numeric value giving the dimensionless total sediment transport
#'   rate (q_star). No critical threshold — returns a non-zero value for any
#'   positive shear stress. Note that this includes suspended load and will
#'   overestimate bedload transport when used in a bedload-only context.
#' @references
#'   Engelund, F. and Hansen, E. (1967). \emph{A Monograph on Sediment
#'   Transport in Alluvial Streams}. Teknisk Forlag, Copenhagen, Denmark.
#' @seealso \code{\link{u_ferg}}, \code{\link{tau_1D}}, \code{\link{t_star}}
#' @export
eh <- function(R, S, d, d84 = 2 * d, a1 = 6.5, a2 = 2.5,
               g = 9.81, rho = 1000, rho_s = 2650) {

  tau <- tau_1D(R, S, g, rho)
  t_s <- t_star(tau, d, g, rho, rho_s)

  # Darcy-Weisbach friction factor from Ferguson velocity equation
  # f = 8 * g * R * S / U^2
  U   <- u_ferg(R, S, d84, a1, a2, g)
  f   <- 8 * g * R * S / U^2

  0.05 / f * t_s^2.5
}
