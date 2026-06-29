# NMF Rank Stability Revision Analysis

Generated: 2026-06-28T17:36:59+00:00

## Purpose

This sensitivity analysis tests whether the CAF-core spatial ecotype result depends on arbitrarily fixing NMF rank to 4.

## Design

- Input: 143 spatial samples by 16 CAF-core enrichment features from `spatial_ecotype_sample_summary.csv`.
- Ranks tested: 2-8.
- Repeats per rank: 50 randomized NNDSVDar initializations plus one NNDSVDa reference fit matching the original Stage 22 implementation.
- Metrics: reconstruction explained fraction, adjusted Rand index across runs, consensus PAC, consensus cophenetic correlation and component cosine reproducibility.

## Rank Summary

- Rank 2: explained 0.960; delta from previous nan; ARI 0.991; PAC 0.014; component cosine 1.000; converged 50/50; best-run smallest component n=23.
- Rank 3: explained 0.971; delta from previous 0.011; ARI 0.998; PAC 0.003; component cosine 1.000; converged 50/50; best-run smallest component n=14.
- Rank 4: explained 0.980; delta from previous 0.009; ARI 0.993; PAC 0.012; component cosine 1.000; converged 50/50; best-run smallest component n=4.
- Rank 5: explained 0.985; delta from previous 0.005; ARI 0.997; PAC 0.003; component cosine 1.000; converged 50/50; best-run smallest component n=3.
- Rank 6: explained 0.989; delta from previous 0.004; ARI 0.993; PAC 0.013; component cosine 1.000; converged 50/50; best-run smallest component n=1.
- Rank 7: explained 0.991; delta from previous 0.003; ARI 0.981; PAC 0.020; component cosine 1.000; converged 50/50; best-run smallest component n=0.
- Rank 8: explained 0.993; delta from previous 0.002; ARI 0.968; PAC 0.024; component cosine 1.000; converged 50/50; best-run smallest component n=0.

## NNDSVDa Reference Sweep

- Rank 2: explained 0.960; delta from previous nan; converged=True; smallest component n=18.
- Rank 3: explained 0.971; delta from previous 0.011; converged=True; smallest component n=8.
- Rank 4: explained 0.980; delta from previous 0.009; converged=True; smallest component n=4.
- Rank 5: explained 0.985; delta from previous 0.005; converged=True; smallest component n=5.
- Rank 6: explained 0.989; delta from previous 0.004; converged=True; smallest component n=2.
- Rank 7: explained 0.991; delta from previous 0.003; converged=True; smallest component n=0.
- Rank 8: explained 0.993; delta from previous 0.002; converged=True; smallest component n=0.

## Rank-4 Component Labels

- rank4_NMF1: DC/APC / immune core / IFN/MHC.
- rank4_NMF2: EMT/invasion / hypoxia / myCAF.
- rank4_NMF3: plasma cell / immune core / T cell.
- rank4_NMF4: basal-like / tumor aggressive / hypoxia.

## Interpretation

Rank 4 explains 0.980 of the non-negative CAF-core feature matrix on average across randomized starts and 0.980 in the NNDSVDa reference fit, with ARI 0.993, PAC 0.012 and component cosine reproducibility 1.000.
The most stable assignment rank by ARI/cosine is rank 3, and the lowest-PAC consensus rank is rank 3; therefore rank 4 should be presented as a biologically interpretable working resolution rather than as the mathematically unique optimum.
The manuscript-safe claim is that the main CAF-core axes are robust across a rank sweep, while the exact number of ecotype labels is a resolution choice supported by reconstruction, stability and interpretability.

## Manuscript Use

- Add as an Extended Data or Supplementary Figure supporting rank choice.
- In Methods, state that NMF ranks 2-8 were rerun with 50 random initializations per rank.
- In Results, avoid saying four ecotypes are intrinsically discrete; say rank 4 resolved four recurrent axes used for downstream interpretation.
