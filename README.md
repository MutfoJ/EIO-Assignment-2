# EIO Assignment 2: Market Power and Differentiated Goods

**Authors:** Dragos Florin Vasile & Wong Hei Wong  
**Course:** Empirical Industrial Organization, TSE (M1)  
**Date:** March 2026

## Overview

We estimate a nested logit demand model (Berry, 1994) for the EU car market, recover multi-product markups, and simulate a BMW-Mercedes merger using the `mergersim` package.

## Structure

```
PS2.pdf                       # Problem set
cars_ps.dta                   # EU car market dataset
Final/
  assignment2_final.do        # Stata code
  assignment2_final.log       # Stata output log
  assignment2_final.md        # Report (Markdown source)
  assignment2_final.pdf       # Report (PDF)
```

## Requirements

- Stata 18 with `ivreghdfe`, `reghdfe`, `ftools`, `ivreg2`, `ranktest`, and `mergersim` (st0349)
- Pandoc + XeLaTeX for PDF generation from Markdown
