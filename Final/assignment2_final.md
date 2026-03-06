---
title: "Market Power and Differentiated Goods"
subtitle: "Empirical Industrial Organization -- Assignment 2"
author: "Dragos Florin Vasile & Wong Hei Wong"
date: "March 2026"
geometry: margin=0.9in
fontsize: 11pt
linestretch: 1.10
urlcolor: blue
header-includes: |
  \usepackage{booktabs}
  \usepackage{float}
  \floatplacement{table}{H}
---

# Model Specification

We estimate a one-level nested logit demand model following Berry (1994) on a panel of the EU car market covering five countries and multiple years. Each consumer chooses at most one car to maximize utility, with the outside option being not to purchase; the total number of potential consumers in each market is $M = \text{population}/3$. Cars are grouped into three nests defined by the `class` variable: small, medium, and luxury.

Applying the Berry inversion, the estimating equation is:
$$\ln(s_j) - \ln(s_0) = \text{FE}_{co} - \alpha\, p_j + \sum_g \sigma_g \ln(s_{j|g})\cdot\mathbf{1}[j\in g] + \xi_j$$
where $\text{FE}_{co}$ are model (brand$\times$model) fixed effects that absorb all time-invariant product characteristics like horsepower, fuel efficiency, and size. The key parameters are the price sensitivity $\alpha > 0$ and three class-specific nesting parameters $\sigma_g \in [0,1]$. A larger $\sigma_g$ means consumers substitute more within that class relative to switching across classes. Standard errors are clustered at the market (country$\times$year) level.

# Identification and Instruments

Once we absorb model fixed effects, the parameters $\alpha$ and $\sigma_g$ are identified from variation in prices and within-group shares of the *same* car model across different countries and years. This relies on the idea that, say, the VW Golf faces different competitive environments in France versus Germany, leading to different prices and group shares even though the car itself is essentially the same.

Estimating this by OLS is problematic for two reasons. First, prices are endogenous: a positive demand shock $\xi_j$ raises sales, which leads firms to charge higher prices, biasing the price coefficient upward and making demand appear less elastic. Second, $\ln s_{j|g}$ is endogenous for similar reasons: a positive shock to one product increases its group share, biasing $\sigma_g$ toward one.

To address this, we construct BLP-style instruments (Berry, Levinsohn, and Pakes, 1995). For each product, we compute sums of observed characteristics (horsepower, fuel, width, height) across: *(i)* rival firms' products, *(ii)* own-firm other products, *(iii)* same-group products from other firms, and *(iv)* own-firm same-group products, plus product counts, totaling 20 excluded instruments. These sums capture the competitive pressure a product faces without directly entering utility. Instruments *(iii)* and *(iv)* are especially important for identifying $\sigma_g$, since they capture within-group competitive intensity.

# Estimation Results

**Table 1: OLS vs. IV Estimates**

| Parameter | OLS | IV (2SLS) |
|:--------------------------|---------------------:|---------------------:|
| $\alpha$ (price)          | $-0.951\;(0.053)$   | $-1.231\;(0.136)$   |
| $\sigma_{\text{small}}$   | $0.980\;(0.003)$    | $0.948\;(0.028)$    |
| $\sigma_{\text{medium}}$  | $0.963\;(0.005)$    | $0.946\;(0.074)$    |
| $\sigma_{\text{luxury}}$  | $0.747\;(0.025)$    | $-0.286\;(0.332)$   |
| $N$                       | 11,483               | 11,468               |
| $R^2$                     | 0.958                | 0.867                |

*Standard errors clustered by market in parentheses.*

Comparing the two columns confirms the endogeneity problems discussed in Section 2. The IV price coefficient is substantially larger in absolute value than OLS ($-1.231$ vs. $-0.951$): once we instrument for price, demand is revealed to be more elastic, exactly as we would expect if OLS was biased by the positive correlation between prices and unobserved quality. Similarly, the nesting parameters for small and medium cars decrease under IV ($\sigma_{\text{small}}$: $0.980 \to 0.948$; $\sigma_{\text{medium}}$: $0.963 \to 0.946$), consistent with upward endogeneity bias in the within-group share.

The luxury $\sigma$ becomes negative and insignificant under IV ($-0.286$, s.e. $0.332$). This is not uncommon for small groups; Brenkers and Verboven (2006) obtain similar results, and it means the luxury segment effectively behaves as standard logit. We bound it at zero for elasticities and markups.

For the first stage, the F-statistics on the excluded instruments are 8.57 (price), 21.15 (small), 14.41 (medium), and 6.01 (luxury). The Kleibergen--Paap Wald $F$-statistic is 5.86, indicating some weak instrument concern, particularly for price and luxury.

# Elasticities

We derive the elasticities from the IV estimates, bounding $\sigma_{\text{luxury}}$ at zero. Following Berry (1994), the own-price elasticity for product $j$ in group $g$ is:
$$\eta_{jj} = \frac{\alpha\, p_j}{1-\sigma_g}\bigl[1 - \sigma_g s_{j|g} - (1-\sigma_g) s_j\bigr]$$
For products $k$ and $j$ in the *same* group, the cross-price elasticity is $\eta_{kj} = \frac{\alpha\, p_j}{1-\sigma_g}[\sigma_g s_{j|g} + (1-\sigma_g) s_j]$, while for products in *different* groups it simplifies to $\eta_{kj} = \alpha\, p_j\, s_j$. This is the key feature of the nested logit: products within the same class are much closer substitutes than products across classes.

**Table 3: Mean Elasticities (IV Estimates)**

| Firm | Own-price | Cross (same grp) | Cross (diff. grp) |
|:---------|----------:|------------------:|-------------------:|
| BMW      |  $-15.00$ |           $0.418$ |           $0.001$  |
| Mercedes |   $-2.63$ |           $0.015$ |           $0.003$  |

BMW's much higher own-price elasticity (in absolute value) largely reflects its cheaper Rover-branded products, which have high within-group shares and thus are very sensitive to price changes. Mercedes products are more differentiated in the luxury segment and face less elastic demand. The within-group cross-elasticities are orders of magnitude larger than the across-group ones, confirming the nested logit's ability to capture realistic substitution patterns that the simple logit would miss.

# Markups: BMW and Germany, 1999

To compute markups, we use the `mergersim` package (Bjornerstedt and Verboven, 2014), which recovers marginal costs from the multi-product Bertrand--Nash first-order conditions:
$$c = p + (\Theta \odot \Delta(p))^{-1}\, q(p)$$
where $\Theta$ is the ownership matrix and $\Delta(p)$ the matrix of demand derivatives. This accounts for multi-product firms like BMW internalizing cannibalization between their own products when setting prices. Since `mergersim` requires a single nesting parameter, we use a sales-weighted average $\hat\sigma = 0.854$. The table below shows the distribution of markups for all BMW products sold in Germany in 1999:

**Table 4: BMW Product Markups (Germany, 1999)**

| Product | Class | Price | Markup | Lerner |
|:---------------------|:---------|------:|-------:|-------:|
| Rover Mini           | Small    | 0.445 |  0.120 |  0.270 |
| Rover 200            | Small    | 0.511 |  0.120 |  0.235 |
| Rover 400            | Small    | 0.582 |  0.120 |  0.206 |
| BMW 3                | Medium   | 0.778 |  0.138 |  0.177 |
| Rover RH(620,623)    | Medium   | 0.820 |  0.138 |  0.168 |
| Rover 75             | Medium   | 0.957 |  0.138 |  0.144 |
| BMW 5                | Luxury   | 1.285 |  0.148 |  0.115 |

BMW's mean Lerner index is 18.8%. A clear pattern emerges: cheaper Rover products have *higher* relative markups (21--27%) because they face less elastic demand in their segments, while the expensive BMW 5 has a lower Lerner of 11.5%. For comparison, Mercedes products in Germany in 1999 average a Lerner of 27.0%, with the C Klasse at 29.1% and the E Klasse at 22.3%. Mercedes' higher markups reflect its stronger market position in the luxury segment.

# Merger Simulation: BMW--Mercedes in Germany, 1999

We simulate a hypothetical merger between BMW and Mercedes in Germany 1999 using `mergersim`. The program: (1) recovers pre-merger marginal costs from the Bertrand--Nash FOCs, (2) modifies the ownership matrix so BMW and Mercedes products belong to one firm, changing which cross-price effects are internalized, and (3) solves for the new post-merger equilibrium prices using Newton's method, holding marginal costs constant (or adjusting them for cost efficiencies). We consider two scenarios: no efficiencies, and a 10% marginal cost reduction for all products of the merged firm.

The pre-merger HHI in Germany 1999 is 1,636; post-merger it rises to 1,810, giving a $\Delta$HHI of 174.

**Table 5: Merger Effects, No Efficiencies vs. 10% Cost Reduction**

| Metric | No Efficiencies | 10% Cost Reduction |
|:------------------------------------|-----------------:|-------------------:|
| Mean $\Delta p$ (merging firms)     |         $+3.8\%$ |           $-5.3\%$ |
| Mean $\Delta p$ (all firms)         |         $+0.4\%$ |           $-0.8\%$ |
| BMW 5 $\Delta p$                    |        $+14.9\%$ |          $+11.2\%$ |
| Mercedes C Klasse $\Delta p$        |         $+9.8\%$ |           $+7.8\%$ |
| $\Delta$ Consumer surplus           |      $-31{,}525$ |         $-1{,}905$ |
| $\Delta$ Producer surplus           |      $+17{,}679$ |        $+37{,}900$ |

Without efficiencies, the merged entity internalizes the competition between BMW and Mercedes, leading to price increases that are concentrated in the luxury segment where both firms directly compete: the BMW 5 rises by 14.9% and the Mercedes C Klasse by 9.8%. Products in the small and medium classes see more modest increases (1--4%), since merger-related market power effects are weaker there. Consumer surplus falls by 31,525 units while producer surplus rises by 17,679, yielding a net total welfare loss of roughly 13,846.

With 10% cost reductions, the picture changes considerably: average prices for the merging firms actually *fall* by 5.3%, and most products become cheaper. However, the luxury segment remains problematic: BMW 5 and Mercedes C Klasse still see price increases of 11.2% and 7.8% respectively, because the market power effect in that segment outweighs the cost savings. Consumer surplus still drops, though much less ($-1{,}905$), while producer surplus rises substantially ($+37{,}900$).

**Should the authority clear the merger?** Without efficiencies, the answer is clearly no: both consumers and total welfare are worse off. Even with a generous 10% cost reduction, consumer surplus still falls and luxury-segment prices rise significantly. The $\Delta$HHI of 174 exceeds standard thresholds that trigger scrutiny. A competition authority should block this merger unless the firms can credibly demonstrate cost efficiencies, particularly in the luxury segment where competitive overlap is greatest. At minimum, behavioral remedies or divestitures of overlapping luxury models would be needed to offset the anticompetitive effects.
