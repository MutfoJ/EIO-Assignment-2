*==============================================================================
* EIO Assignment 2: Market Power and Differentiated Goods
* Authors: Dragos Florin Vasile & Wong Hei Wong
* Date: March 2026
* Dataset: EU Car Market (cars_ps.dta)
*
* One-level nested logit demand (Berry 1994) with class-varying nesting
* parameters. Multi-product Bertrand-Nash markups and merger simulation
* using the mergersim package (Bjornerstedt & Verboven, 2014).
*==============================================================================

clear all
set more off
cap log close

global datadir "c:\Users\vdrag\OneDrive\Escritorio\Estudios\TSE\Semestre 2\EIO\EIO Assignment 2"
global outdir "$datadir\Final"
cap mkdir "$outdir"

log using "$outdir\assignment2.log", replace

*==============================================================================
* SECTION 1: DATA PREPARATION
*==============================================================================

use "$datadir\cars_ps.dta", clear

* Market = country x year
egen mkt = group(country year)

* Market size: households = pop / 3
gen M = pop / 3

* Market shares
bysort mkt: egen total_q = total(qu)
gen sj = qu / M
gen s0 = (M - total_q) / M
assert s0 > 0

* Dependent: ln(s_j / s_0)
gen lnsj_s0 = ln(sj) - ln(s0)

* Within-group shares (class: 1=small, 2=medium, 3=luxury)
bysort mkt class: egen group_q = total(qu)
gen sjg = qu / group_q

* Class-specific ln(s_{j|g}) for class-varying sigma
gen ln_sjg_small  = ln(sjg) * (class == 1)
gen ln_sjg_medium = ln(sjg) * (class == 2)
gen ln_sjg_luxury = ln(sjg) * (class == 3)

*==============================================================================
* SECTION 2: BLP INSTRUMENTS
*==============================================================================

* Sums of product characteristics (Berry, Levinsohn, Pakes 1995)
local chars "horsepower fuel width height"
foreach var of local chars {
    bysort mkt: egen sum_`var'_mkt = total(`var')
    bysort mkt firm: egen sum_`var'_firm = total(`var')
    bysort mkt class: egen sum_`var'_grp = total(`var')
    bysort mkt firm class: egen sum_`var'_fgrp = total(`var')

    gen iv_`var'_rival = sum_`var'_mkt - sum_`var'_firm
    gen iv_`var'_own   = sum_`var'_firm - `var'
    gen iv_`var'_grp   = sum_`var'_grp - `var'
    gen iv_`var'_fgrp  = sum_`var'_fgrp - `var'
}

bysort mkt: gen n_mkt = _N
bysort mkt firm: gen n_firm_mkt = _N
bysort mkt class: gen n_grp = _N
bysort mkt firm class: gen n_fgrp = _N

gen iv_n_rival = n_mkt - n_firm_mkt
gen iv_n_own   = n_firm_mkt - 1
gen iv_n_grp   = n_grp - 1
gen iv_n_fgrp  = n_fgrp - 1

global ivlist "iv_horsepower_rival iv_fuel_rival iv_width_rival iv_height_rival iv_horsepower_own iv_fuel_own iv_width_own iv_height_own iv_horsepower_grp iv_fuel_grp iv_width_grp iv_height_grp iv_horsepower_fgrp iv_fuel_fgrp iv_width_fgrp iv_height_fgrp iv_n_rival iv_n_own iv_n_grp iv_n_fgrp"

*==============================================================================
* EXERCISE 3: ESTIMATION
*==============================================================================

* ---- 3a. OLS with model fixed effects ----
di _n "============================================="
di "EXERCISE 3a: OLS ESTIMATION"
di "============================================="

areg lnsj_s0 princ ln_sjg_small ln_sjg_medium ln_sjg_luxury, ///
    absorb(co) cluster(mkt)
estimates store ols

* ---- 3b. IV (2SLS with absorbed model FE) ----
di _n "============================================="
di "EXERCISE 3b: IV ESTIMATION (2SLS with model FE)"
di "============================================="

ivreghdfe lnsj_s0 (princ ln_sjg_small ln_sjg_medium ln_sjg_luxury = $ivlist), ///
    absorb(co) cluster(mkt)
estimates store iv

scalar alpha_iv    = _b[princ]
scalar sigma_s_iv  = _b[ln_sjg_small]
scalar sigma_m_iv  = _b[ln_sjg_medium]
scalar sigma_l_iv  = _b[ln_sjg_luxury]

di _n "IV estimates:"
di "alpha       = " alpha_iv
di "sigma_small = " sigma_s_iv
di "sigma_med   = " sigma_m_iv
di "sigma_lux   = " sigma_l_iv

* ---- 3c. First stage regressions ----
di _n "============================================="
di "FIRST STAGE REGRESSIONS"
di "============================================="

foreach depvar in princ ln_sjg_small ln_sjg_medium ln_sjg_luxury {
    di _n "--- First stage for `depvar' ---"
    areg `depvar' $ivlist, absorb(co) cluster(mkt)
    testparm $ivlist
    di "Joint F-stat = " r(F) "  p-value = " r(p)
}

* ---- OLS vs IV comparison ----
di _n "============================================="
di "OLS vs IV COMPARISON"
di "============================================="
estimates table ols iv, se stats(N r2) b(%9.4f) se(%9.4f)

*==============================================================================
* EXERCISE 4: ELASTICITIES (using IV estimates)
*==============================================================================

di _n "============================================="
di "EXERCISE 4: OWN AND CROSS-PRICE ELASTICITIES"
di "============================================="

* Class-specific sigma (bound luxury at 0 if negative)
gen sigma_g = sigma_s_iv * (class == 1) ///
            + sigma_m_iv * (class == 2) ///
            + max(0, sigma_l_iv) * (class == 3)

* Own-price elasticity (Berry 1994, nested logit):
* eta_jj = -[alpha/(1-sigma_g)] * [1 - sigma_g*s_{j|g} - (1-sigma_g)*s_j] * p_j
gen own_elas = (alpha_iv / (1 - sigma_g)) * ///
              (1 - sigma_g * sjg - (1 - sigma_g) * sj) * princ

* Cross-price: same group
* eta_kj = -alpha * p_j * [sigma_g/(1-sigma_g)*s_{j|g} + s_j]
gen cross_same = (-alpha_iv) * princ * ///
                (sigma_g / (1 - sigma_g) * sjg + sj)

* Cross-price: different group
* eta_kj = -alpha * p_j * s_j
gen cross_diff = (-alpha_iv) * princ * sj

di _n "--- BMW (firm==2) ---"
di "Own-price elasticities:"
sum own_elas if firm == 2
di "Cross-price (same group):"
sum cross_same if firm == 2
di "Cross-price (diff group):"
sum cross_diff if firm == 2

di _n "--- Mercedes (firm==12) ---"
di "Own-price elasticities:"
sum own_elas if firm == 12
di "Cross-price (same group):"
sum cross_same if firm == 12
di "Cross-price (diff group):"
sum cross_diff if firm == 12

*==============================================================================
* EXERCISE 5: MARKUPS - BMW & GERMANY 1999
*==============================================================================

di _n "============================================="
di "EXERCISE 5: MARKUPS - BMW & GERMANY 1999"
di "============================================="

* Multi-product Bertrand-Nash markups via mergersim.
* mergersim recovers marginal costs by inverting the ownership-weighted
* demand Jacobian: c = p + (Theta odot Delta(p))^{-1} * q(p)
* This properly accounts for firms owning multiple products.

xtset co mkt

* mergersim accepts a single sigma; use sales-weighted average of bounded sigmas
sum sigma_g [w=qu] if country == 3 & year == 1999
scalar sigma_avg = r(mean)
di "Calibrated alpha = " alpha_iv
di "Sales-weighted sigma = " sigma_avg

mergersim init, nests(class) price(princ) quantity(qu) marketsize(M) ///
    alpha(`=alpha_iv') sigmas(`=sigma_avg')

mergersim market if country == 3 & year == 1999, firm(firm)

* Multi-product markup = price - recovered marginal cost
gen markup = princ - M_costs if country == 3 & year == 1999
gen lerner = markup / princ if country == 3 & year == 1999

di _n "BMW products in Germany 1999:"
list type class princ M_costs markup lerner ///
    if firm == 2 & country == 3 & year == 1999, noobs

di _n "Summary of BMW markups in Germany 1999:"
sum markup lerner if firm == 2 & country == 3 & year == 1999, detail

di _n "Mercedes products in Germany 1999:"
list type class princ M_costs markup lerner ///
    if firm == 12 & country == 3 & year == 1999, noobs

di _n "All firms in Germany 1999:"
tabstat markup lerner if country == 3 & year == 1999, ///
    by(firm) stat(mean sd min max) format(%9.4f) nototal

*==============================================================================
* EXERCISE 6: MERGER SIMULATION (BMW & MERCEDES, GERMANY 1999)
*==============================================================================

di _n "============================================="
di "EXERCISE 6: MERGER SIMULATION"
di "============================================="

* ---- Pre- and post-merger HHI ----
preserve
keep if country == 3 & year == 1999
bysort firm: egen firm_sales = total(qu)
egen total_sales = total(qu)
gen firm_share = firm_sales / total_sales * 100
collapse (first) firm_share, by(firm)
gen firm_share_sq = firm_share^2
qui sum firm_share_sq
scalar hhi_pre = r(sum)
qui sum firm_share if firm == 2
scalar share_bmw = r(mean)
qui sum firm_share if firm == 12
scalar share_merc = r(mean)
scalar hhi_post = hhi_pre - share_bmw^2 - share_merc^2 ///
                + (share_bmw + share_merc)^2
di "Pre-merger HHI  = " hhi_pre
di "Post-merger HHI = " hhi_post
di "Delta HHI       = " hhi_post - hhi_pre
restore

* ---- Scenario 1: No cost efficiencies ----
di _n "============================================="
di "SCENARIO 1: NO COST EFFICIENCIES"
di "============================================="

cap drop M_* markup lerner

mergersim init, nests(class) price(princ) quantity(qu) marketsize(M) ///
    alpha(`=alpha_iv') sigmas(`=sigma_avg')

mergersim market if country == 3 & year == 1999, firm(firm)

di _n "Pre-merger marginal costs for BMW and Mercedes:"
list type firm class princ M_costs ///
    if country == 3 & year == 1999 & (firm == 2 | firm == 12), noobs

mergersim simulate if country == 3 & year == 1999, ///
    firm(firm) buyer(2) seller(12) detail keepvars

gen price_ch_pct = (M_price_ch / princ) * 100 ///
    if country == 3 & year == 1999

di _n "Post-merger price changes (no efficiencies):"
list type firm class princ M_price_ch price_ch_pct ///
    if country == 3 & year == 1999 & (firm == 2 | firm == 12), noobs

di "Mean price increase for merging firms:"
sum price_ch_pct if country == 3 & year == 1999 & (firm == 2 | firm == 12)
di "Mean price increase for all firms:"
sum price_ch_pct if country == 3 & year == 1999

* ---- Scenario 2: 10% cost efficiencies for both firms ----
di _n "============================================="
di "SCENARIO 2: 10% COST EFFICIENCIES"
di "============================================="

cap drop M_* price_ch_pct

mergersim init, nests(class) price(princ) quantity(qu) marketsize(M) ///
    alpha(`=alpha_iv') sigmas(`=sigma_avg')

mergersim market if country == 3 & year == 1999, firm(firm)

mergersim simulate if country == 3 & year == 1999, ///
    firm(firm) buyer(2) seller(12) ///
    buyereff(0.10) sellereff(0.10) detail keepvars

gen price_ch_pct = (M_price_ch / princ) * 100 ///
    if country == 3 & year == 1999

di _n "Post-merger price changes (10% efficiencies):"
list type firm class princ M_price_ch price_ch_pct ///
    if country == 3 & year == 1999 & (firm == 2 | firm == 12), noobs

di "Mean price change for merging firms:"
sum price_ch_pct if country == 3 & year == 1999 & (firm == 2 | firm == 12)
di "Mean price change for all firms:"
sum price_ch_pct if country == 3 & year == 1999

*==============================================================================
di _n "Assignment 2 completed."
log close
