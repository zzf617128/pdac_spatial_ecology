# Results reorganization audit

Date: 2026-06-29

## Purpose

The manuscript-modification task list requested that the Results section be tightened into four coherent modules rather than reading as a stack of many separate analyses.

## Completed Reorganization

The active revised manuscript now uses four Results modules:

1. `CAF-myeloid cores define a reproducible spatial architecture in PDAC`
2. `Metastatic-site remodeling suggests lymph-node immune decoupling`
3. `Cell-state and interface analyses nominate CAF-core biological axes`
4. `Independent validation and pathology context bound the model`

## Structural Change

- Previous Results structure: six Results subsections, including separate sections for deep spatial architecture and H&E morphology.
- Revised Results structure: four Results subsections aligned to the task-list story modules.
- Approximate Results length after reorganization: 1,888 word-like tokens across 19 paragraphs.

## What Changed Scientifically

- Spatial specificity, threshold sensitivity, contiguous-null, alternative-anchor and marker-overlap analyses were consolidated into the first module.
- GSE272362 primary/liver/LN results, all-LN maps and leave-one-out sensitivity were consolidated into the metastatic-site module.
- NMF ecotypes, candidate axes, targeted genes, marker/reference/NNLS support, mechanism prioritization, core-to-interface transitions and geometry were consolidated into the cell-state/interface module.
- GSE274557 Visium, GSE274673 Xenium, H&E morphology, TCGA context and TLS stress testing were consolidated into the independent validation/pathology-context module.

## Boundary Preservation

The revised Results retains the key boundaries:

- lymph-node metastasis remains a five-sample, hypothesis-generating lead;
- mechanism analyses prioritize candidate axes but do not establish causal signaling;
- GSE274557 does not validate lymph-node-specific immune decoupling;
- GSE274673 is targeted-panel cell-resolution support, not whole-transcriptome abundance;
- H&E and TCGA analyses provide exploratory context, not clinical-grade prediction or spatial validation;
- TLS stress testing does not support mature TLS reframing.

## Files

- Revised manuscript source: `results/revision_2026_06_29/manuscript/Manuscript_NatureSubjournal_revised.md`
- Revised manuscript Word export: `results/revision_2026_06_29/manuscript/Manuscript_NatureSubjournal_revised.docx`
