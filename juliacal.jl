using CSV
using Interpolations
using SpecialFunctions

"""
p95, p68, calprob, medage = juliacal(c14age, c14err, calcurve, yeartype)

juliacal v 0.1
2019-10-30
bryan.lougheed@geo.uu.se

Radiocarbon (14C) calibration in Julia. Probabilistic age calibration with credible
intervals using Bayesian highest posterior density (HPD). Based on the Matlab/Octave 
function 'matcal' (Lougheed and Obrochta 2016), which was ported to Julia by
Bryan Lougheed in 2019.

Will produce credible intervals, median age and an array of age probabilities.
Calibration plot not yet available, still a work in progress!

Required input
==============
c14age:    Lab 14C determination in 14C yr BP.

c14err:    Lab 14C determination uncertainty (1 sigma) in 14C yr.

calcurve:  String specifying calibration curve to use, select from
            the following (not case sensitive):
            "IntCal13", "Marine13", "SHCal13", "IntCal09", "Marine09"
            "IntCal04", "Marine04", "SHCal04", "IntCal98", "Marine98"

yeartype:  String specifying how to report calibrated age.
           Choices are "CalBP" or "BCE/CE". (Not case sensitive)

Optional input parameters
=========================

resage:    Optional (parameter name and value). Specify reservoir
           age in 14C yr. R(t) in the case of atmospheric calibration
           curve, delta-R in the case of marine curve. (default = 0)
           e.g. resage=320 for a reservoir age of 320

reserr:    Optional (parameter name and value). Specify a 1 sigma
           uncertainty for your chosen resage (default = 0)
           e.g. reserr=50 for an uncertainty of 50

"""
function juliacal(c14age, c14err, calcurve, yeartype; resage=0, reserr=0)

# Process reservoir age
if isnan(resage) == true
	resage = 0
end
if isnan(reserr) == true
	reserr = 0
end
c14ageorig = c14age
c14errorig = c14err
c14age = c14age - resage
c14err = sqrt(c14err^2 + reserr^2)

# convert to f14c space
f14age = exp(c14age/-8033)
f14err = f14age*c14err/8033

# import calibration curve
d = CSV.read("IntCal13.14c", delim=",", datarow=12, header=["calbp","c14age","c14err","d14c","d14cerr"])
curvecal = reverse(d.calbp)
curve14c = reverse(d.c14age)
curve14cerr = reverse(d.c14err)
curvef14 = exp.(curve14c./-8033) # convert to f14c space
curvef14err = curvef14.*curve14cerr./8033 # convert to f14c space

# linear interp
interpres = 1
hicurvecal = collect(minimum(curvecal):interpres:maximum(curvecal))
hicurvef14 = interpolate((curvecal,), curvef14, Gridded(Linear()))(hicurvecal)
hicurvef14err = interpolate((curvecal,), curvef14err, Gridded(Linear()))(hicurvecal)

# calculate probability for every hicurvecal year in F14 space
calprob = fill(NaN,length(hicurvecal),2)
calprob[:,1] = hicurvecal
# equation from e.g. p.261 in Bronk Ramsey, 2008. doi:10.1111/j.1475-4754.2008.00394.x
calprob[:,2] = exp.(-((f14age .- hicurvef14).^2)./(2 .* (f14err^2 .+ hicurvef14err.^2))) ./ ((f14err^2 .+ hicurvef14err.^2).^0.5)
calprob[:,2] = calprob[:,2] ./ sum(calprob[:,2]) # 

# do HPD to find 1 sig and 2 sig credible intervals
hpd = collect(calprob)
hpd = hpd[sortperm(hpd[:,2]),:] # sort rows by second column
hpd = [hpd fill(NaN,length(hicurvecal),1)]
hpd[:,3] = cumsum(hpd[:,2])

# 1 sig
hpd68 = hpd[hpd[:,3] .>= 1-erf(1/sqrt(2)), :]
hpd68 = hpd68[sortperm(hpd68[:,1]),:] # sort by first column
ind1 = findall(diff(hpd68[:,1]) .> 1)
if isempty(ind1) == true
	p68 = fill(NaN,1,3)
	p68[1,1] = hpd68[end,1]
	p68[1,2] = hpd68[1,1]
	p68[1,3] = sum(hpd68[1:end,dims=2])
else
	indy1 = fill(NaN,length(ind1)*2,1)
	for i = 1:length(ind1)
		indy1[i*2-1,1] = ind1[i]
		indy1[i*2,1] = ind1[i]+1
	end
	indy1 = trunc.(Int,[ 1  indy1' length(hpd68[:,1]) ])
	p68 = fill(NaN,length(2:2:length(indy1)),3)
	for i = 2:2:length(indy1)
		p68[trunc(Int,i/2),1] = hpd68[indy1[i],1]
		p68[trunc(Int,i/2),2] = hpd68[indy1[i-1],1]
		p68[trunc(Int,i/2),3] = sum(hpd68[indy1[i-1]:indy1[i],2])
	end
	p68 = reverse(p68, dims=1)
end

# 2 sig
hpd95 = hpd[hpd[:,3] .>= 1-erf(2/sqrt(2)), :]
hpd95 = hpd95[sortperm(hpd95[:,1]),:] # sort by first column
ind1 = findall(diff(hpd95[:,1]) .> 1)
if isempty(ind1) == true
	p95 = fill(NaN,1,3)
	p95[1,1] = hpd95[end,1]
	p95[1,2] = hpd95[1,1]
	p95[1,3] = sum(hpd95[1:end,dims=2])
else
	indy1 = fill(NaN,length(ind1)*2,1)
	for i = 1:length(ind1)
		indy1[i*2-1,1] = ind1[i]
		indy1[i*2,1] = ind1[i]+1
	end
	indy1 = trunc.(Int,[ 1  indy1' length(hpd95[:,1]) ])
	p95 = fill(NaN,length(2:2:length(indy1)),3)
	for i = 2:2:length(indy1)
		p95[trunc(Int,i/2),1] = hpd95[indy1[i],1]
		p95[trunc(Int,i/2),2] = hpd95[indy1[i-1],1]
		p95[trunc(Int,i/2),3] = sum(hpd95[indy1[i-1]:indy1[i],2])
	end
	p95 = reverse(p95, dims=1)
end

# get median age
_, indmed = findmin(abs.(cumsum(calprob[:,2]).-0.5))
medage = trunc(Int,(calprob[indmed[1],1]))

# convert to BCE/CE if asked for
if (lowercase(yeartype) == "bce/ce") == true
	medage = (medage-1950) * -1
	calprob[:,1] = (calprob[:,1].-1950) .* -1
	p95_4[:,1:2] = (p95_4[:,1:2].-1950) .* -1
	p68_2[:,1:2] = (p68_2[:,1:2].-1950) .* -1
end

return p95, p68, calprob, medage

end # end function


	
