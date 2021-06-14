I am interested in the potential for reproducible papers to be used the peer review process.

A highly-cited paper from 2010, ["Topology of the human and mouse m6A RNA methylomes revealed by m6A-seq"](https://www.nature.com/articles/nature11112), received [13 "contrasting" citations](https://scite.ai/reports/topology-of-the-human-and-WDmMRO?contradicting=true&mentioning=false&page=1&supporting=false) in scite.ai, indicating some groups found different 

- mRNA splicing in METTL3 knockdown cells
- number of peaks
- peak location 
- overall methylation levels

In order to study the possible role of bioinformatic approaches to this problem employed in Dominissini et al., I am proposing a "test of robustness", in which different tools, parameters, references, subsets are explored.

Rather than attempt to faithfully reproduce this paper using its original analysis stack, I have used a couple off-the-shelf pipelines to reanalyze this data:
- https://github.com/eQTL-Catalogue/rnaseq - The workflow processes raw data from FastQ inputs (FastQC, Trim Galore!), aligns the reads (STAR or HiSAT2), generates gene counts (featureCounts, StringTie) and performs extensive quality-control on the results (RSeQC, dupRadar, Preseq, edgeR, MultiQC).
- https://github.com/kingzhuky/meripseqpipe - MeRIP-seq analysis pipeline arranged multiple alignment tools, peakCalling tools, Merge Peaks' methods and methylation analysis methods.


To run:

adjust `max_cpus` and `max_memory` in `meripseqpipe/conf/m6a.config` as needed

```
gsutil -m cp -r gs://truwl-dominissini/SRP012098 .  (or s3://dominissini/SRP012098)
gsutil -m cp -r gs://truwl-dominissini/refs .       (or s3://dominissini/refs)
curl -s https://get.nextflow.io | bash
nextflow run meripseqpipe -profile m6a,docker
```

As a reproducible peer reviewer, you may...

- Manipulating parameters, data subsets
- Swapping out tools, reference data, or even primary data as needed
- Implementing further downstream analyses
- Implement something more similar to the original analysis

...in order to evaluate the robustness of the results presented in this paper!


# Data
Data from SRP012098 (m6a RNA IP and Input for IFNg (200ng/ml) or HGF/SF (10 ng/ml) over night. Stress effects were tested in HepG2 cells by either 30 minutes incubation at 43ºC (heat shock) or UV irradiation)
Data from SRP012096 (METTL3_KD1 RNA-seq and mock controls)

# What has been done so far
- Code to generate manifests for meripseqpipe and rna-seq pipelines (see utils/metautils.py)
- [A RNA-seq comparison of METTL3 knockdown & control HepG2 cells](https://github.com/BarryDigby/GSE37001) by Barry Digby
