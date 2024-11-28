# Supporting Information

Last Updated on 2024-11-28

[![](https://zenodo.org/badge/837544982.svg)](https://zenodo.org/doi/10.5281/zenodo.13358758)

This repository contains the supporting information for: **Comparative
analysis of proteomic expression between oesophageal adenocarcinoma and
normal adjacent tissue**.

This comprises Tables S1-S4 as `csv` files, and the prediction of the
effect of the G95R mutation on NOP58 using DDGun
([1](#ref-montanucci2019)) below.

The column names and contents of the `csv` files in the `tables` folder
are described below.

## Patient Information

S1 Table contains the patient information for the 7 male donors in this
study.

                      File                                  
                    ────────────────────────────────────────
                      S1-Table-OAC-Patient-Information.csv  
                    ────────────────────────────────────────

Column names: File

| Column name          | Description                 |
|----------------------|-----------------------------|
| `donor_id`           | Donor identifier            |
| `age_at_diagnosis`   | Age of patient at diagnosis |
| `sex`                | Sex of patient              |
| `location_of_tumour` | Location of tumour          |
| `treatment_modality` | Treatment modality          |

Patient information

## Peaks normalised Top 3 peptide intensities

Label free quantification using the Peaks Q module of Peaks Studio
([2](#ref-zhang2012),[3](#ref-lin2013)) yielding matrices of protein
identifications as quantified by their normalised top 3 peptide
intensities.

S2 Table contain normalised top 3 peptide intensities.

           File                                                        
         ──────────────────────────────────────────────────────────────
           S2-Table-OAC-Peaks-Normalised-Top3-Peptide-Intensities.csv  
         ──────────────────────────────────────────────────────────────

Column names: File

| Column name | Description |
|-----------------------------------------------|-------------------------|
| `protein` | protein short name |
| `gene` | HGNC gene symbol |
| `sample_id` the donor id or donor id suffixed with `T` for tumour or `N` for NAT samples | Normalised top 3 peptide intensity from Peaks |

Peaks normalised top 3 peptide intensities Table information

## Differential protein expression with DEqMS

The normalised top 3 peptide intensities were filtered to remove any
proteins for which there were more than two missing values across the
samples. Differential protein expression (DEP) was then calculated with
DEqMS using the default steps ([4](#ref-deqms)).

S3 Table contain the output of DEqMS. A negative `logFC` indicates
higher expresssion in OAC tissue and a postive `logFC` indicates higher
expression in NAT.

                         File                            
                       ──────────────────────────────────
                         S3-Table-OAC-DEqMS-Results.csv  
                       ──────────────────────────────────

Column names: File

| Column name    | Description                                    |
|----------------|------------------------------------------------|
| `logFC`        | log2 fold change between two groups            |
| `AveExpr`      | the mean of the log2 ratios across all samples |
| `t`            | Limma t-values                                 |
| `P.Value`      | Limma p-values                                 |
| `adj.P.Val`    | BH method adjusted Limma p-values              |
| `B`            | Limma B values                                 |
| `gene`         | HGNC gene symbol                               |
| `count`        | peptide count values                           |
| `sca.t`        | DEqMS t-statistics                             |
| `sca.P.Value`  | DEqMS p-values                                 |
| `sca.adj.pval` | BH method adjusted DEqMS p-values              |
| `protein`      | protein short name                             |

DEqMS Table information

## Functional analysis with g:Profiler

Functional enrichment analysis used g:Profiler ([5](#ref-kolberg2020))
using default settings for homo sapiens modified to exclude GO
electronic annotations. Protein identifiers were used as inputs for
DEPs.

S4 Table contain the g:Profiler outputs.

                       File                                
                     ──────────────────────────────────────
                       S4-Table-OAC-gProfiler-Results.csv  
                     ──────────────────────────────────────

Column names: File

<table>
<caption>g:Profiler Table information</caption>
<colgroup>
<col style="width: 15%" />
<col style="width: 84%" />
</colgroup>
<thead>
<tr class="header">
<th>Column name</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>query</code></td>
<td>the name of the input query</td>
</tr>
<tr class="even">
<td><code>significant</code></td>
<td>indicator for statistically significant results</td>
</tr>
<tr class="odd">
<td><code>p_value</code></td>
<td>hypergeometric p-value after correction for multiple testing</td>
</tr>
<tr class="even">
<td><code>term_size</code></td>
<td>number of genes that are annotated to the term</td>
</tr>
<tr class="odd">
<td><code>query_size</code></td>
<td>number of genes that were included in the query</td>
</tr>
<tr class="even">
<td><code>intersection_size</code></td>
<td>the number of genes in the input query that are annotated to the
corresponding term</td>
</tr>
<tr class="odd">
<td><code>precision</code></td>
<td>the proportion of genes in the input list that are annotated to the
function (defined as intersection_size/query_size)</td>
</tr>
<tr class="even">
<td><code>recall</code></td>
<td>the proportion of functionally annotated genes that the query
recovers (defined as intersection_size/term_size)</td>
</tr>
<tr class="odd">
<td><code>term_id</code></td>
<td>unique term identifier</td>
</tr>
<tr class="even">
<td><p><code>source</code></p>
<p><code>term_name</code></p></td>
<td>the abbreviation of the data source for the term (e.g. GO:BP)</td>
</tr>
<tr class="odd">
<td><code>effective_domain_size</code></td>
<td>the total number of genes “in the universe” used for the
hypergeometric test</td>
</tr>
<tr class="even">
<td><code>source_order</code></td>
<td>numeric order for the term within its data source</td>
</tr>
<tr class="odd">
<td><code>parents</code></td>
<td>list of term IDs that are hierarchically directly above the term.
For non-hierarchical data sources this points to an artificial root
node.</td>
</tr>
<tr class="even">
<td><code>evidence_codes</code></td>
<td>a lists of all evidence codes for the intersecting genes between
input and the term. The evidences are separated by comma for each
gene.</td>
</tr>
<tr class="odd">
<td><code>intersection</code></td>
<td>a comma separated list of genes from the query that are annotated to
the corresponding term</td>
</tr>
</tbody>
</table>

g:Profiler Table information

## NOP58 G95R mutation DDGun protein stability prediction

The effect of the G95R mutation on Nucleolar protein 58 (NOP58) was
predicted using <a href="#0" style="font-size: 12pt;">DDGun</a>
([1](#ref-montanucci2019)).

The output of the prediction using the amino acid sequence is showing
how G95R decreases the stability of NOP58 is shown below:

    PROTEIN: Protein name.
    VARIANT: Protein mutation.
    CONSERVATION: Frequencies of the wild-type and mutant residues in the mutated positions.
    S_KD: Kyte-Doolittle substitution scores (AAINDEX1:KYTJ820101).
    S_BL: BLOSUM62 substitution score (AAINDEX2:ENS920102).
    S_DDG[SEQ]: Prediction of the DDG for the single amino acid substitutions.
    T_DDG[SEQ]: Global DDG for multiple amino acid substitutions.
    STABILITY[SEQ]: Variation of protein stability

|  |  |  |  |  |  |  |  |  |
|:-----------:|:----:|:-------:|:---:|:---:|:----:|:-------:|:-------:|:---------:|
| PROTEIN | VARIANT | CONSERVATION | S_KD | S_BL | S_PROF | S_DDG\[SEQ\] | T_DDG\[SEQ\] | STABILITY\[SEQ\] |
| Nucleolar protein 58 | G95R | 30.1 | 0.3 | 0.107 | -3.440 | -1.222 | -0.5 | Decrease |

## References

<span class="csl-left-margin">1.
</span><span class="csl-right-inline">Montanucci L, Capriotti E, Frank
Y, Ben-Tal N, Fariselli P. DDGun: an untrained method for the prediction
of protein stability changes upon single and multiple point variations.
BMC Bioinformatics \[Internet\]. 2019 Jul;20(S14). Available from:
<http://dx.doi.org/10.1186/s12859-019-2923-1></span>

<span class="csl-left-margin">2.
</span><span class="csl-right-inline">Zhang J, Xin L, Shan B, Chen W,
Xie M, Yuen D, et al. PEAKS DB: De novo sequencing assisted database
search for sensitive and accurate peptide identification. Molecular &
Cellular Proteomics. 2012;11(4):M111010587. </span>

<span class="csl-left-margin">3.
</span><span class="csl-right-inline">Lin H, He L, Ma B. A combinatorial
approach to the peptide feature matching problem for label-free
quantification. Bioinformatics \[Internet\]. 2013 May 10;29(14):1768–75.
Available from: <http://dx.doi.org/10.1093/bioinformatics/btt274></span>

<span class="csl-left-margin">4.
</span><span class="csl-right-inline">DEqMS \[Internet\]. Available
from: <http://bioconductor.org/packages/DEqMS/></span>

<span class="csl-left-margin">5.
</span><span class="csl-right-inline">Kolberg L, Raudvere U, Kuzmin I,
Vilo J, Peterson H. gprofiler2 – an R package for gene list functional
enrichment analysis and namespace conversion toolset g:Profiler.
F1000Research \[Internet\]. 2020 Nov 17;9:709. Available from:
<http://dx.doi.org/10.12688/f1000research.24956.2></span>
