# COMBAT_eQTL

## Citation

If you use these data, please cite our work:

**Authors:** Jennifer Astley, Andrew Kwok, Benjamin Hollis, COvid-19 Multi-omics Blood ATlas (COMBAT) Consortium, Calliope A. Dendrou, Stephen Sansom, Julian C. Knight, Alexander J Mentzer, Yang Luo, Luke Jostins-Dean  
**Title:** "Single-cell eQTL mapping identifies disease-specific genetic control of COVID-19 infection"  
**Status:** In preparation  

*(DOI will be added once available.)*

## Abstract

The impact of genotype on gene expression can depend on both cellular and organismal context. Here, we leverage an extensive blood atlas of genotyped patients with varying severity of infection produced by the COVID-19 Multi-omics Blood ATlas Consortium to study the role of genetic regulation on gene expression in a context-specific manner. We analyze single-cell transcriptomic and genome-wide genetic data from >600,000 cells and 100 donors. Across 15 cell types, we identify 3,562 independent cis expression quantitative trait loci (eQTLs) in high linkage disequilibrium (r²>0.8) with 538 infectious and inflammatory disease-associated risk variants, including rheumatoid arthritis (RA) and inflammatory bowel disease (IBD). Notably, we find eQTLs absent from a general population dataset, such as *REL*, *IRF5* and *TRAF*, all of which are differentially regulated by infection and whose variants are associated with RA or IBD. We also identify infection-modified eQTLs, including *RPS26* and *ADAM10*, implying that the regulatory sequence context of these genes may play a role in specific immune cell subsets in infection. Our work demonstrates that the overriding effect of genetics on gene expression in blood immune cells is independent of infection status or severity. However, small numbers of eQTLs are modified by infection, and these differences can illustrate potentially important immune biology.

## Repository Structure

- `data/` – Data for plotting main and supplementary figures.  
- `figures/` – Generated plots for main and supplementary figures.  
- `make_manuscript_figures.R` – Script to recreate figures from the manuscript.
