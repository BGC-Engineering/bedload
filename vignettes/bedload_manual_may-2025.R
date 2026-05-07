## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(here)
#source(here('bedload', 'R','bedload_functions.R'))
library(bedload)

## -----------------------------------------------------------------------------
compare_fun(d = 0.05, 
            S = 0.02)

## -----------------------------------------------------------------------------
compare_fun(d = 0.05, 
            S = 0.0002)

## -----------------------------------------------------------------------------

example_section <- make_trapezoid(w_top = 10,
                                  w_bot = 5,
                                  H_bank = 1.5,
                                  Z_fp = 100)
plot(example_section$dist, example_section$Z,
     col = "darkgreen", 
     xlab = "distance (m)",
     ylab = "elevation (m)",
     asp = 1)
example_slope <- 0.05
example_d <- 0.04
example_result <- estimate_hydraulics(xs_data = example_section,
                                      S = example_slope,
                                      d = example_d)
head(example_result)


## -----------------------------------------------------------------------------
# enter data that you have from the field on the site of interest
test_d <- 0.05  # median surface texture, ideally from a Wolman sample
test_d84 <- 0.12  #the 84th percentile of the bed surface
test_slope <- 0.015  # reach average bed slope from lidar or map
top_width <- 12  # bank to bank width at the site
bot_width <- 3.5  # width of the channel bed 
bank_ht <- 1.2  # vertical high of the floodplain above the channel bed

# make a cross section for analysis based on field measurements
#(you can also is load_section() to load a more detailt topo profile)
test_section <- make_trapezoid(w_top = top_width,
                               w_bot = bot_width,
                               H_bank = bank_ht,
                               Z_fp = 100)

# use the scs_hydrograph() function to produce a hydrograph
test_Q_peak <- 12 # in cubic meters per second
test_time_lag <- 6 #in hours
test_dt <- 15/60 #produce discharge output every 15 minutes
test_hydrograph <- scs_hydrograph(Qpeak = test_Q_peak,
                                  tpeak = test_time_lag, 
                                  dt = test_dt)
plot(test_hydrograph)

# now we are ready to run estimate_hydraulics...
test_df <- estimate_hydraulics(xs_data = test_section,
                               S = test_slope,
                               d = test_d,
                               d84 = test_d84,
                               ymax = bank_ht * 2,
                               fun = "mpm"
                               )
# NOTE: the optional parameter ymax was assigned a value twice the channel depth...this is not strictly necessary, but ymax controls the range of water depths (and discharge values) analysed by the function, and the next step requires a data frame with discharge values at least as high as the peak flow being assessed.

# now we can run the integration over the hydrograph
test_output <- estimate_capacity(xs_hyd = test_df, 
                                 Q = test_hydrograph$y,
                                 t = test_hydrograph$x)

print(test_output)

# and we can show the peak water surface elevation on the cross section to check that things make sense

plot(test_section, asp = 2)
abline(h = test_output$wse, col = "blue", lty = 2)

## -----------------------------------------------------------------------------
#example gsd data
field_sizes <- c(16, 23, 32, 45, 64, 91, 128, 181, 256, 360)
field_counts <-c(0, 4, 7, 15, 18, 20, 17, 11, 9, 5)

test_gsd <- wolman_gsd(Di = field_sizes, 
                  counts = field_counts)

#assume a percent sand for the bed surface
test_fsand <- 8 / 100  

test_output <- fractional_capacity(hydrograph = test_hydrograph,
                                    xs_data = test_section,
                                    S = test_slope,
                                    gsd = test_gsd,
                                    Fsand = test_fsand)

print(test_output)

## -----------------------------------------------------------------------------
total_transport_volume <- sum(test_output$Qbi)
print(round(total_transport_volume))


## -----------------------------------------------------------------------------
plot(test_output$dm, rev(cumsum(rev(test_output$Fi))),
     type = "b",
     lty = 2,
     pch = 20,
     log = "x",
     col = "darkgreen",
     ylim = c(0, 1),
     xlab = "Grain size (mm)",
     ylab = "Cumulative proportion finer")
points(test_output$dm, rev(cumsum(rev(test_output$Li))),
       pch = 20,
      col = "darkorange")
lines(test_output$dm, rev(cumsum(rev(test_output$Li))),
      lty = 2,
      col = "darkorange")

## -----------------------------------------------------------------------------
#pick a D50 and log-normal standard deviation that gives the same D50 and D84 as the Wolman sample above
D50 <- approx(x = test_output$cpf,
              y = test_output$Di,
              xout = 0.5)[[2]]
sp <- 1.08 

test_gsd <- sim_gsd(D50, sp)

test_output <- fractional_capacity(hydrograph = test_hydrograph,
                                    xs_data = test_section,
                                    S = test_slope,
                                    gsd = test_gsd,
                                    Fsand = test_fsand)
total_transport_volume <- sum(test_output$Qbi)
print(round(total_transport_volume))

plot(test_output$dm, rev(cumsum(rev(test_output$Fi))),
     type = "b",
     lty = 2,
     pch = 20,
     log = "x",
     col = "darkgreen",
     ylim = c(0, 1),
     xlab = "Grain size (mm)",
     ylab = "Cumulative proportion finer")
points(test_output$dm, rev(cumsum(rev(test_output$Li))),
       pch = 20,
      col = "darkorange")
lines(test_output$dm, rev(cumsum(rev(test_output$Li))),
      lty = 2,
      col = "darkorange")

## -----------------------------------------------------------------------------
#specify flood details
flood_mag <- 40
flood_dur <- 6

#specify river reach details
river_n <- 0.04
river_d50 <- 0.04
river_d84 <- 2.1 * river_d50
river_w0 <- 25
river_slope <- 0.007
river_H <- 0.5

#run the model
out_put <- gbem(flood_mag, flood_dur, river_n, river_d84, river_d50, river_w0, river_slope, river_H)

#update the river width based on the modelled bank erosion
river_w1 <- river_w0 + out_put[2]
print(river_w1)

## -----------------------------------------------------------------------------
#change the bank strength
river_H <- river_H / 2

#run the model using the stable width from the last model run
out_put <- gbem(flood_mag, flood_dur, river_n, river_d84, river_d50, river_w1, river_slope, H = river_H)
#update the river width based on the modelled bank erosion
river_w2 <- river_w1 + out_put[1]
print(river_w2)

## -----------------------------------------------------------------------------
#specify the flood hydrograph
flood_peak <- 60
time2peak <- 6
dt <- 0.25
test_hydrograph <- scs_hydrograph(flood_peak, time2peak, dt)[[2]]

#specify river reach details
river_n <- 0.04
river_d50 <- 0.04
river_d84 <- 2.1 * river_d50
river_w0 <- 25
river_slope <- 0.007
river_H <- 0.5

out_put <- bank_erode_fp(test_hydrograph, dt, river_n, river_d84, river_d50, river_w0, river_slope, river_H)
head(out_put)

## -----------------------------------------------------------------------------
total_width_increase <- sum(out_put$erosion)
print(total_width_increase)

total_transport_volume <- sum(out_put$Vb)
print(total_transport_volume)

plot(out_put$Q, out_put$Vb,
     col = "darkorange",
     pch = 21,
     type = "b",
     xlab = expression(Discharge~(m^3/s)),
     ylab = expression(Transport~volume~(m^3/s)))
                       

## -----------------------------------------------------------------------------
#repeat analysis with higher bank strength
river_H <- 5
out_put <- bank_erode_fp(test_hydrograph, dt, river_n, river_d84, river_d50, river_w0, river_slope, river_H)

total_transport_volume <- sum(out_put$Vb)
print(total_transport_volume)

plot(out_put$Q, out_put$Vb,
     col = "darkorange",
     pch = 21,
     type = "b",
     xlab = expression(Discharge~(m^3/s)),
     ylab = expression(Transport~volume~(m^3/s)))
                       


## -----------------------------------------------------------------------------
#use the flood hydrograph created above for bank_erode_fp


#use the same river reach details as for bank_erode_fp
river_H <- 0.5
river_Z <- 2.5

out_put <- bank_erode_trc(test_hydrograph, dt, river_n, river_d84, river_d50, river_w0, river_slope, river_H, river_Z)
head(out_put)

## -----------------------------------------------------------------------------
total_width_increase <- sum(out_put$erosion)
print(total_width_increase)

total_transport_volume <- sum(out_put$Vb)
print(total_transport_volume)

plot(out_put$Q, out_put$Vb,
     col = "darkorange",
     pch = 21,
     type = "b",
     xlab = expression(Discharge~(m^3/s)),
     ylab = expression(Transport~volume~(m^3/s)))
                       

## -----------------------------------------------------------------------------
#simulate the hydrograph
Q <- scs_hydrograph(Qpeak = 7, 
                    tpeak = 8, 
                    dt = 1)[[2]]
q <- Q / 6  #channel is approximately 6 m wide
d <- 0.6  #estimated bankfull depth
D50 <- 0.05
Slope <- 0.02

#use same values for stream evaluated in bank_erode_pf() and bank_erode_trc()
step_length(D50, q, d, Slope, 1, flashy = FALSE, use_ferguson = TRUE)

## -----------------------------------------------------------------------------
# if running for the first time, also run...
# tidyhydat::download_hydat()

wsc_test <- quick_ffa_ann('08GA022')

print(wsc_test)

usgs_test <-quick_ffa_ann('12150800')

print(usgs_test)

## -----------------------------------------------------------------------------
# if running for the first time, also run...
# tidyhydat::download_hydat()

wsc_test <- quick_ffa_dual('08MF065',
                           start = '04-01',
                           end = '07-01')

print(wsc_test)

usgs_test <-quick_ffa_dual('12451000',
                           start = '04-01',
                           end = '07-01')

print(usgs_test)

