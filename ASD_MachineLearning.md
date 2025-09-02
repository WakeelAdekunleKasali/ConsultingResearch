Autism Spectrum Disorder (High Dimensional Analysis with Machine Learning) - Project Report
================

2025-04-02



## Introduction

Autism spectrum disorder (ASD) is a heterogeneous neurodevelopmental
condition characterized by a wide range of cognitive, behavioral, and
neurological differences. Despite significant advances in genetic and
molecular research, the precise cellular mechanisms underlying ASD
remain incompletely understood (Charman et al., 2011). ASD arises during
prenatal development, with disruptions occurring in two major
developmental epochs. From the first to third trimesters, alterations in
cell proliferation, neurogenesis, neuronal migration, and fate
determination have been identified, while in the third trimester and
early postnatal life, cortical wiring, synaptogenesis, and neural
network organization are disrupted (Courchesne et al., 2020). Studies of
postmortem ASD brains reveal increased cortical neuron numbers, abnormal
laminar organization, and dysregulated synaptic connectivity, supporting
the hypothesis that ASD results from altered prenatal corticogenesis
(Courchesne et al., 2019). The differences associated with ASD has also
been found to upregulate pathways related to stress, which may also play
a part in the resulting phenotype (Pangrazzi et al. 2020) Large genomic
studies with cohorts of over 35 584 samples have identified 103 ASD risk
genes, many of which play crucial roles in neurodevelopment (Satterstrom
et al., 2021). Amongst these are genes related to chromatin
accessibility, transcriptional regulation and synaptic function where
dysregulation of excitatory and inhibitory signalling appears as a
disease hallmark, also in organoids (De Rubeis et al., 2014). Building
off of these discoveries, induced pluripotent stem cell (iPSC)-derived
brain organoid models have emerged as a transformative approach for
studying the neurodevelopment in ASD. These three-dimensional models
recapitulate key aspects of cortical development, allowing researchers
to investigate patient-specific cellular and molecular phenotypes, as
can be seen in figure 1. Organoid studies have already provided insights
into ASD-related alterations, including changes in progenitor
proliferation, premature neuronal differentiation, and synaptic
dysfunction (Eichmüller et al., 2022).

<figure>
<img src="../images/organoid_development.jpg"
alt="Figure 1 - Comparion of cerebral organoid development compared to an actual human brain, adapted from Kelava et al. (2016)" />
<figcaption aria-hidden="true"><strong>Figure 1</strong> - Comparion of
cerebral organoid development compared to an actual human brain, adapted
from Kelava et al. (2016)</figcaption>
</figure>

Recent research by Han et al. highlighted a downregulation of fatty
acid-binding protein 7 (FABP7) in cerebral organoids derived from
normocephalic non-regressive ASD patients. This reduction in FABP7 led
to changes in Mitogen-Activated Protein Kinase Kinase 1/2 (MEK1/2)
phosphorylation, resulting in premature neuronal differentiation. FABP7
is a crucial regulator of radial glial function and neuroepithelial
maintenance, and its dysregulation may contribute to early neurogenic
defects in ASD (Han et al., 2025). In thread with this, excessive
neurogenesis has been reported in ASD patients with macrocephaly,
whereas normocephalic ASD cases exhibit premature differentiation
instead (Han et al., 2025). Suggesting that different ASD subtypes may
follow different neurodevelopmental pathways.

Despite the application of cerebral organoids in ASD research, a major
gap remains in understanding how intercellular communication changes
over developmental time. Our project aims to address this gap by
analyzing single-cell RNA sequencing (scRNA-seq) data from normocephalic
ASD and control cerebral organoids from Day 0 to Day 100. By utilizing
tools such as CellChatV2, we will investigate alterations in cell-cell
communication networks across developmental stages to provide novel
insights into the molecular and cellular underpinnings of normocephalic
ASD. Cerebral organoid models do also present challenges, including the
formation of necrotic cores due to diffusion limitations and potential
stress responses in deeper cell layers (Eichmüller et al., 2022). These
factors must be carefully considered when interpreting the
transcriptomic data. Nonetheless, by applying advanced computational
analyses and focusing on a well-characterized noncephalic ASD organoid
dataset, this project will contribute to a more nuanced understanding of
ASD neurodevelopment and its cellular interactions across development.

## Research question

Han et al. performed a study on autism spectrum disorder (ASD) organoids
where they found premature neuronal differentiation to be caused by
FABP7. The aim of this project is to investigate how intercellular
communication networks are altered in normocephalic ASD during early
cortical development using single-cell RNA sequencing data from cerebral
organoids. By applying CellChatV2 and related analytical tools, we seek
to identify timepoint-specific differences in neurodevelopmental
signaling between ASD and control organoids from Day 0 to Day 100.

## Methods

### Pre-processing

Before going into analysis, the data need to be pre-processed. To do
this we created Seurat objects from the data retrieved from the publicly
available dataset in the Gene Expression Omnibus (GEO) at
<https://www.ncbi.nlm.nih.gov/geo/> with accession number GSE203201. The
Seurat objects were created by following the pipeline described by
Satijalab in the following protocol
<https://satijalab.org/seurat/articles/pbmc3k_tutorial> (version 5.2.1).
The protocol differed in the normalization and scaling method, where
SCTransform was used instead
<https://satijalab.org/seurat/articles/sctransform_vignette>, and we
used SC-type for clustering
<https://github.com/IanevskiAleksandr/sc-type> followed by conservative
manual curation of annotations with canonical marker genes for neurons
and neural progenitor cells. This pre-processing deviates from the paper
originally published by Han et al.

### Pseudotime

While the original study (Han et al., 2025) processed the trajectory
analysis and pseudotemporal ordering of cells with the monocle2 package
of R (version 3.4.0) to characterize the potential functional changes
between ASD and IMR90-4 (a standard iPSC line), we are interested in
that if similar patterns could be observed when comparing with multiple
standard iPSC lines. Here, we used monocle 3 (version 1.3.7, R version
4.4.3, package guideline
<https://cole-trapnell-lab.github.io/monocle3/docs/getting_started/> )
developed by Cao et al. to identify developmental trajectories of ASD
and multiple control cerebral organoids (NC3-1, NC1B-3, and IMR90-4)
from 6 time points ranging from Day 0 to Day 100. Monocle3 is an
advanced package tailored for the analysis of scRNA-seq data, with
particular proficiency in trajectory inference and pseudotime ordering.
It initially maps cells into a low-dimensional space that encodes
transcriptional states utilizing UMAP. Subsequently, it clusters
mutually analogous cells with the Louvain community discovery technique
and consolidates nearby groups into ‘supergroups’. Ultimately, it
delineates the pathways or trajectories that individual cells may follow
during development, pinpointing the sites of bifurcations and
convergences within each supergroup.

**Preprocessing specified for monocle3**  
In this study, we first preprocessed each cell sample using seurat as
depicted before and merged 12 Seurat objects into two combined seurat
objects: one for ASD samples and one for control controls. And those 2
objects are converted to cell_data_set (CDS) objects which are direct
input for monocle pipeline. We applied principal component analysis
(PCA) and Uniform Manifold Approximation and Projection (UMAP) to reduce
the data dimensionality. The number of principal components was set to
30 (num_dim = 30), based on the balance between computational efficiency
and the capacity to capture transcriptional heterogeneity.

**Trajectory Inference and Pseudotime Ordering**  
A trajectory graph was learned using monocle3’s learn_graph() function
with the parameter of use_partition setting as FALSE, as we aimed to
infer a global developmental trajectory across all clusters without
partitioning. To ensure biological interpretability, the root node was
automatically assigned to cells derived from Day 0 (D0) samples, which
represent the earliest developmental timepoint.

**GO (gene ontology) Enrichment Analysis**  
To elucidate the biological significance of pseudotime-associated
transcriptional changes, we performed GO enrichment analysis on
differentially expressed genes (DEGs) identified along the developmental
trajectory using the clusterProfiler (version 4.14.6) package of R
(version 4.4.3). We did this because we wanted to see whether some
specifical biological pathways/functions would be closely related to
what we would figure out in the trajectory inference.

R script for this analysis can be found at
[`src/pseudotime_and_go_enrichment_analysis.R`](../src/pseudotime_and_go_enrichment_analysis.R).

### CellChat

**Cell-cell communication analysis**  
To investigate intercellular communication patterns, we used the
CellChat R package (v2.1.2). A comprehensive protocol by its developers
Jin et al. can be found at
<https://www.nature.com/articles/s41596-024-01045-4>. Seurat objects
were subset by condition (ASD vs Control) and developmental time point
(D12, D30, D60) and converted into CellChat objects with
createCellChat() and filterCommunication() parameter min.cells set to
15. The code for loading CellChat objects canbe found at
[src/create_cellchat_objects.R](https://github.com/STAT540-UBC-2025/project-amethyst/blob/main/src/create_cellchat_objects.R)
. Signaling pathway inference was performed using the
computeCommunProb() and computeCommunProbPathway() functions after
identifying overexpressed ligands and receptors. Merging CellChat
objects for comparison was done using mergeCellChat(). For each time
point, ASD and Control CellChat objects were merged and compared using
netVisual_bubble() and plotGeneExpression() to identify global
differences in signaling. The code for generating these plots can be
found at
[src/pathway_analysis_script.R](https://github.com/STAT540-UBC-2025/project-amethyst/blob/main/src/pathway_analysis_script.R).

**Differential Expression Analysis**  
To identify ligand-receptor interactions driving differences between ASD
and Control, we used two strategies: Communication Probability
Comparison: netVisual_bubble() was used to compare interaction strength
between groups. Differential Gene Expression (DEA): We ran
identifyOverExpressedGenes() on merged CellChat objects, followed by
netMappingDEG() and subsetCommunication() to extract ligand-receptor
pairs with significant changes in expression in sender or receiver
cells. This was done using both cell type-specific DEA and combined DEA
(group.DE.combined = TRUE).

**Pathway Prioritization**  
Selected pathways were prioritized using: Overall information flow
(summarized signaling strength across all cell-cell interactions) from
rankNet(). Statistical significance (Wilcoxon test) of pathway-level
differences. Biological relevance to neurodevelopment (e.g.,
proliferation, differentiation, synapse formation).

### E/I ratio

To quantify the excitatory/inhibitory (E/I) balance within each sample,
we developed a custom function named compute_ei_ratio(). This function
classifies each cell as excitatory, inhibitory, or unclassified based on
the expression of curated gene markers. Excitatory markers correspond to
projection neurons within the glutamatergic pathway, while inhibitory
markers identify GABAergic interneurons. Marker gene selection was
guided by established literature on cortical development, and a fixed
set of genes was used across all samples to ensure comparability.

Before classification, the function filters the input marker sets to
include only those genes that are present in the Seurat object, ensuring
robustness to missing features and preventing errors during downstream
analyses. Expression values for the retained markers are extracted using
FetchData(), targeting a specified assay and data layer (e.g., “data”
slot of the “RNA” assay).Cells are categorized as follows:

- Excitatory: Cells expressing at least one excitatory marker (non-zero
  expression).
- Inhibitory: Cells expressing no excitatory markers but at least one
  inhibitory marker.
- Unclassified: Cells expressing none of the specified markers.

The function then calculates the E/I ratio as the number of excitatory
cells divided by the number of inhibitory cells within each sample. This
scalar value provides a summary measure of the excitatory/inhibitory
balance based on marker-defined cell identities. The compute_ei_ratio()
function was applied to Seurat objects from Day 30, Day 60, and Day 100
organoid samples, covering both ASD and Control groups. Resulting E/I
ratios were used to assess developmental changes in excitatory and
inhibitory balance, and to evaluate potential shifts between ASD and
Control conditions at each time point. The code can be accessed via
[src/Excitatory_Inhibitory_ratio.R](https://github.com/STAT540-UBC-2025/project-amethyst/blob/main/src/Excitatory_Inhibitory_ratio.R)

## Results

### ASD organoids show different pseudotime trajectory compared to Control

Trajectory inference was effectively generated based on the single-cell
transcriptome profiles of ASD and control cell lines (Figure 1). The
computed pseudotime signifies the relative extent of each cell’s
advancement along the developmental continuum, from undifferentiated
progenitor state to differentiated terminal cell type.

<figure>
<img src="../images/trajectory_ASD.png"
alt="Figure 2A - Trajectory inference result of ASD cells D0-D100. UMAP plot displaying the inferred developmental trajectory of single cells from ASD samples, colored by pseudotime. Black circles indicate trajectory graph nodes. Cells are ordered along the trajectory from early (purple) to late (yellow) pseudotime states. The graph structure represents inferred lineage relationships and potential branching events during cellular development." />
<figcaption aria-hidden="true"><strong>Figure 2A</strong> - Trajectory
inference result of ASD cells D0-D100. UMAP plot displaying the inferred
developmental trajectory of single cells from ASD samples, colored by
pseudotime. Black circles indicate trajectory graph nodes. Cells are
ordered along the trajectory from early (purple) to late (yellow)
pseudotime states. The graph structure represents inferred lineage
relationships and potential branching events during cellular
development.</figcaption>
</figure>

<figure>
<img src="../images/trajectory_CRTL.png"
alt="Figure 2B - rajectory inference result of control cells D0-D100" />
<figcaption aria-hidden="true"><strong>Figure 2B</strong> - rajectory
inference result of control cells D0-D100</figcaption>
</figure>

The control (CRTL) samples displayed compact and organized trajectories,
featuring slightly less branching nodes and a clear progression from Day
0 cells to differentiated stages (Figure 1b). Conversely, the ASD
samples had longer, more scattered branches (Figure 1a). Notably, at the
mid-stage of pseudotime estimation especially, ASD cells send out longer
and more dispersed branches, which seems to indicate earlier
differentiation. Because the further a cell progresses along with its
trajectory, the more its gene expression pattern diverges from its
starting state, and the more dramatically the cell’s function and
morphology change.

These data indicate that ASD brain progenitor cells differentiate
prematurely and diverge from the standard developmental pathway,
aligning with the original study (Han et al., 2025). Besides replicating
the trajectory analytic methodology from the original study, our
trajectory analysis integrated three distinct control iPSC-derived brain
organoid lines, instead of depending on a singular standard iPSC line as
utilized in prior research. This method not only strengthens the
reliability of trajectory comparison between ASD and control but also
offers additional evidence for the early differentiation and modified
developmental trajectories in ASD organoids (Han et al., 2025).

The results of the GO enrichment are summarized in Figure 2. ASD and
control organoids generally exhibited similar biological processes
related to neurodevelopment, including the regulation of axon and
neurite projection development and mitochondrial organization. These
common terms suggest that core developmental programs are preserved in
both cases. What is different is that ASD samples showed significant
enrichment of pathways related to RNA metabolism, organelle
localization, and protein localization to organelles (Figure 2a). But
there’s no report of direct connections between ASD and these GO terms.

<figure>
<img src="../images/GO_enrichment_ASD.png"
alt="Figure 3A - GO Enrichment result of ASD cells over all time points. The top enriched biological processes are shown, with dot size representing the number of genes involved and color indicating the adjusted p-value (Benjamini–Hochberg corrected)." />
<figcaption aria-hidden="true"><strong>Figure 3A</strong> - GO
Enrichment result of ASD cells over all time points. The top enriched
biological processes are shown, with dot size representing the number of
genes involved and color indicating the adjusted p-value
(Benjamini–Hochberg corrected).</figcaption>
</figure>

<figure>
<img src="../images/GO_enrichment_CRTL.png"
alt="Figure 3B - GO Enrichment result of control cells over all time points." />
<figcaption aria-hidden="true"><strong>Figure 3B</strong> - GO
Enrichment result of control cells over all time points.</figcaption>
</figure>

In contrast, the control organoids exhibited a rich array of biological
processes consistent with normal neurosecretion and cell cycle
progression, including mitosis, cell cycle phase transitions,
cytokinesis, and the regulation of neurogenesis, reflecting a more
morphological and more standardized developmental program (figure 2b).

In summary, our pseudotime analysis results indicate that ASD organoids
display abnormal developmental trajectories marked by early and varied
differentiation patterns, potentially linked to the hypothesis of
premature differentiation. But the mechanism cannot be adequately
explained by GO enrichment results.

#### Overview of differential cell-cell communication between ASD and Control

To assess the global differences of intercellular communication between
ASD and control, we analyzed the differential number of interactions
using CellChat’s comparison framework. The heatmap below visualizes
where ASD show increased (red) or decreased (blue) signaling
interactions compared to control, which serves as the reference in this
case. This is shown for three timepoints, D12, D30, and D60. The earlier
timepoints were excluded due to poor data for getting meaningful
CellChat results (too few cell types and poor annotations).

<figure>
<img src="../images/diff_interactions.jpg"
alt="Figure 4 - Differential interactions between ASD and Control for different celltypes across time (D12, D30, and D60)" />
<figcaption aria-hidden="true"><strong>Figure 4</strong> - Differential
interactions between ASD and Control for different celltypes across time
(D12, D30, and D60)</figcaption>
</figure>

In the earliest timepoint we included, we see that the ASD organoids are
exhibiting increased signaling in Neural progenitor cells (NPCs),
particularly to other NPCs. Additionally, we see that there are enhanced
interactions between neuroepithelial cells and NPCs to multiple other
cell types. This could possibly suggest an early dysregulation of
progenitor network coordination, supporting the initial hypothesis about
premature differentiation of cells in ASD.

At day 30, we observe a that the differential signaling has become more
widespread. There is an overall increase in interactions in ASD,
including immature neurons, neuroepithelial cells, and radila glial
cells. Neuroblasts are showing a more heterogenous pattern, with some
upregulation and downregulation, dependig on the target cell type. These
broad changes suggest a sustained divergence in intercellular
communication that could be explained by an accelerated of dysfunctional
differentiation process.

For day 60, the most prominent change is the increase in neuroblast
signaling in ASD. In addition, there is some upregulation of signaling
in GABAergic and immature neurons. This aligns with the hypothesis that
there is an imbalance in excitatory/inhibitory cell types.

When we view these data in a bigger perspective across time. Some of the
patterns we see are that the nature of differential expression changes
over time. It shifts from NPCs in D12, to neuroblasts and neurons at
D60. We also see that progenitor and neuroepithelial cells show
consistent and sustained alterations in communication in ASD. From this
we can generally hypothesize that ASD organoids diverge from
communication trajectories in Control. It also maintains distinct
signaling dynamics as new cell types emerge over time. This gives us an
overview for our further analyses, where we will start off by going into
specific signaling pathways.

#### Altered Signaling Pathways Across Developmental Timepoints

The global changes in signaling pathways between ASD and Control
organoids can be analyzed by looking at the relative information flow
for the most significant pathways in D12, D30, and D60. Pathways have
been ordered by effect size, with the statistically significant
differences calculated with a paired Wilcoxon test. The bars have been
normalized to a sum of 1 per pathway, and colored basd on condition,
where red represents Control, and blue is ASD.

<figure>
<img src="../images/pathways.jpg"
alt="Figure 5 - Relative information flow in pathways for Control (red) and ASD (blue) at D12, D30, and D60. The pathways are ordered by effect size, calculated with paired Wilcoxon test." />
<figcaption aria-hidden="true"><strong>Figure 5</strong> - Relative
information flow in pathways for Control (red) and ASD (blue) at D12,
D30, and D60. The pathways are ordered by effect size, calculated with
paired Wilcoxon test.</figcaption>
</figure>

ASD Organoids at day 12 are exhibiting upregulation of pathways related
to synapses and neurogenesis. Some of those related to neurodevelopment
have been highlighted with purple arrows in the figure below. One of the
highlighted pathways includes NOTCH, which is necessary for neural
progenitor self renewal (Imayoshi, 2010). Glutamate is another one, and
is associated with excitatory synaptic signaling (Willard, 2013). NRXN
is the last one we will highlight for ASD, and it is involved in
synaptic adhesion and network formation (Gomez, 2021). The pathways that
are enriched in Control compared to ASD are WNT which is necessary for
self-renewal and proliferation of neural progenitor cells (Nusse, 2008)
As well as CADM which is involved in synaptic adhesion and network
formation (Fujita, 2012). Together, these findings suggest that ASD
organoids might be shifting towards neurogenic and synaptic adhesion
prematurely compared to Control, which maintains progenitor signaling
and developmental adhesion for longer.

At day 30, differences between ASD and Control are becoming more
evident. The ASD organoids are showing higher WNT signaling, which is
promoting neurogenesis, and possibly accelerates differentiation. GABA-A
is also upregulated in ASD, while GABA-B is downregulated. This keeps
reflecting the E/I imbalance we are suspectinng in ASD. Control on the
other hand has higher FGF and BMP signaling, both of which are key
pathways in supporting neural progenitor renewal (Pera, 2003). NOTCH is
also still increased, reinforcing progenitor maintenance. In other
words, we might be sseing that ASD organoids are transitioning earlier
from progenitor to neural states, also known as premature
differentiation.

For day 60 organoids, we obserev an increase in APP and cholesterol
pathways in ASD. These are linked to neurodegeneration and membrane
remodeling. WNT and NOTCH are also suddenly more expressed in ASD,
remaining active beyond their typical developmental windows. GABA-A is
also still higher in ASD. So there is clearly some dusregulation of
neurodevelopment here.

In the bigger picture, what we see is that ASD organoids demonstrate
premature and potentially dysregulated differentiation. This is
reflected by increased signaling in pathways promoting synaptic function
and neurogenesis, and reduced reliance on classical developmental
regulators. We have chosen some specific pathways at each timepoint to
look further into, and this will be our continued analysis.

#### Altered cell–cell communication via PTN, WNT, and EFNB pathways in ASD organoids

<div style="font-size: 80%; line-height: 1.4;">

<figure>
<img src="../images/pathways_d12_d30_d60.jpg"
alt="Figure 6 - Analysis of ligand–receptor signaling differences between control and ASD cerebral organoids at three developmental stages (D12, D30, D60).(A, C, E) CellChat bubble plots showing predicted cell–cell interactions via the PTN, WNT, and EFNB pathways. Circle size indicates statistical significance, and color denotes communication probability. (B, D, F) Violin plots showing gene expression of ligands and receptors within relevant sender/receiver cell types. Control and ASD conditions are shown in red and teal, respectively." />
<figcaption aria-hidden="true"><strong>Figure 6</strong> - Analysis of
ligand–receptor signaling differences between control and ASD cerebral
organoids at three developmental stages (D12, D30, D60).(A, C, E)
CellChat bubble plots showing predicted cell–cell interactions via the
PTN, WNT, and EFNB pathways. Circle size indicates statistical
significance, and color denotes communication probability. (B, D, F)
Violin plots showing gene expression of ligands and receptors within
relevant sender/receiver cell types. Control and ASD conditions are
shown in red and teal, respectively.</figcaption>
</figure>

</div>

To investigate the emergence of altered intercellular signaling in ASD,
we analyzed ligand–receptor interactions at developmental stages D12,
D30, and D60 using CellChat. Among signaling pathways showing consistent
ASD-specific changes, we focused on PTN, WNT, and EFNB for detailed
analysis (Figure X).

PTN signaling was reduced in ASD organoids, with lower communication
probabilities observed among immature neurons and neural progenitor
cells at D12 (Figure XA). This was accompanied by decreased expression
of PTN and dysregulation of its receptors in ASD compared to control
(Figure XB).In contrast, WNT signaling was elevated in ASD, particularly
from neuroepithelial and radial glial cells at D30 (Figure XC).
Expression of WNT5 ligands was remarkably similar, but the FZD1 and FZD3
receptors showed higher expression in ASD radial glia and
neuroepithelial cells (Figure XD), suggesting aberrant activation of
non-canonical WNT pathways. EFNB–EPHB interactions were also markedly
dysregulated at D60 in ASD organoids (Figure XE), especially among
neuroblasts and neuroepithelial cells. Upregulation of the EFNB3 ligand
and EPHB2/4 receptors in ASD immature neurons and neuroepithelial cells
supports enhanced ephrin signaling in critical neural progenitors
(Figure XF).

Together, these findings highlight a shift in signaling dynamics across
development in ASD organoids, characterized by early loss of
neurotrophic PTN signaling and later-stage dysregulation of ncWNT and
upregulation of EPHB pathways, potentially contributing to aberrant
neurodevelopmental patterning.

### E/I ratio

By history, excitatory neurons in the brain have been linked to a number
of marker genes in many scientific literature. On the list are SLC17A6,
TBR1, SATB2 and FEZF2 which specifically plays role in controlling genes
involved in corticospinal motor neuron identity (Price et al., 2023).
The relative changes in the excitatory and inhibitory neurons during
cortical development are essential for deciphering the growth in ASD
phenotypes. In this study, we apply the concept of excitatory-inhibitory
(E/I ratio) on the cerebral organoid samples derived from induced
pluripotent stem cells. The estimates from this metric are used in
providing results to the research objectives highlighted as follows:

##### Investigate whether normocephalic ASD Organoids exhibit early excitatory dominance during cortical development

Given the well-established excitatory (e.g., SLC17A6, TBR1) and
inhibitory (e.g., GAD1, DLX2) markers used to measure the E/I ratio at
Day 30, 60 and 100, we compare the trajectories of E/I ratio between ASD
and healthy controls groups across the developmental time. From the plot
in (Figure 6), the ASD organoids showed an elevated E/I ratio (~12.4)
compared to controls (~7.9). This supports the early excitatory
expansion in ASD organoids. Such expansion suggests dysregulation of
radial glia and early differentiation as Han et al. (2025) link this to
FABP7 downregulation.

![](../WakeelWorkspace/Excitatory_Inhibitory_ratio_files/figure-gfm/unnamed-chunk-7-1.png)

##### Assess Whether the temporal evolution of E/I balance in ASD differs from controls.

In control samples, E/I ratios gradually dropped over time. This steady
decline is consistent with with GABAergic neuron integration and circuit
stabilization. In contrast, ASD samples, however, dropped sharply
between D30 and D60, then slightly increased at D100. This is a pointer
to disruption in the maturation process. From (Figure 6) it can be
concluded that the E/I trajectories is non-linear between the two groups
across the timepoint.

![](../WakeelWorkspace/Excitatory_Inhibitory_ratio_files/figure-gfm/unnamed-chunk-8-1.png)

##### To quantify the divergence in E/I ratio between the two groups as a marker of disease-specific circuit dysregulation

The divergence plot in (Figure 7), is aimed at capturing information
about ASD-specific shifts in E/I balance. The result showed a large
positive gap at D30 (+4.5), suggesting excess excitation in ASD, a small
negative dip at D60 (−0.2), indicating transient inhibition possibility
or almost convergence in the excitatory-inhibitory signal neuron.
Meanwhile, at D100, divergence returned to positive (~+1.8) revealing an
inconsistent E/I inbalance between the two groups in study.

As a whole, this part of the study highlight a possible inconsistent
behavior in the excitatory-inhibitory balance in the system
neuronal-signal system of the ASD patients during the cortical
development. This analogy further describes the complexity in
understanding the neuronal development common to ASD phenotypes.

## Discussion

The trajectory inference results indicate that ASD brain progenitor
cells differentiate prematurely and diverge from the standard
developmental pathway, aligning with the original study (Han et al.,
2025). Besides replicating the trajectory analytic methodology from the
original study, our trajectory analysis integrated three distinct
control iPSC-derived brain organoid lines, instead of depending on a
singular standard iPSC line as utilized in prior research. This method
not only strengthens the reliability of trajectory comparison between
ASD and control but also offers additional evidence for the early
differentiation and modified developmental trajectories in ASD organoids
(Han et al., 2025).

The CellChat pathway analysis aimed to explore how cell-cell
communication may differ between ASD and control cerebral organoids at
developmental stages D12, D30, and D60. These timepoints were selected
based on cell composition, and them being representative for the neural
development from neural progenitor cells to the formation of mature
neuronal patterns. The results from this analysis highlight potential
disruptions in signaling pathways central to neurodevelopmental timing
and synaptic circuit formation.

We observed in the altered signaling pathways that ASD organoids
exhibited early and sustained upregulation of pathways associated with
neurogenesis and synaptic function (e.g., NOTCH, WNT, GABA-A, NRXN),
alongside reduced activity in classical developmental regulators such as
FGF, BMP, and GABA-B. These shifts reflect a pattern of premature
differentiation, with ASD organoids transitioning earlier from
progenitor to neuronal states. This trajectory is further supported by
changes in cell-type-specific signaling, where progenitor populations
such as NPCs and neuroepithelial cells show altered interaction profiles
in ASD as early as D12, and mature neuronal populations dominate
signaling by D60.

At D12, PTN signaling was notably reduced in ASD, particularly in
immature neurons. PTN regulates neurogenesis by promoting neurite
outgrowth, supporting survival, and organizing radial glial scaffolding
(Wang et al, 2020). Reduced PTN in ASD may contribute to accelerated
depletion of the progenitor pool, aligning with previously observed
premature differentiation phenotypes. At D30, we observed elevated
non-canonical WNT (ncWNT) signaling in ASD, especially among radial glia
and neuroepithelial cells. While ncWNT is less studied than the
canonical β-catenin pathway, it is known to regulate early polarity,
migration, and fate decisions (Suthon et al., 2021). This upregulation
could further reinforce premature developmental transitions or atypical
differentiation. Although this interpretation fits the trajectory of
accelerated development, the role of ncWNT in this specific context
remains less well defined and should be considered
hypothesis-generating. By D60, we observed reduced EPHB signaling in ASD
organoids. Ephrin-B signaling supports neuronal migration, axon
guidance, and synapse formation. Its downregulation may reflect impaired
neural circuit assembly and could reinforce the disorganized
developmental phenotype suggested by early differentiation signals (Lai
et al., 2009).

It is essential to recognize that these interpretations are heavily
reliant on cell type annotations, which were derived through a
combination of SC-type and manual curation. The ambiguity of cell
identity—particularly among transitional cell states like neural
progenitors and neuroblasts—could impact pathway assignment and
downstream inferences. Furthermore, CellChat’s analytical framework is
predictive, not definitive. It assumes that mRNA expression of ligands
and receptors translates directly to functional protein-level signaling,
does not account for post-translational modifications, and relies on a
curated interaction database. Therefore, the pathway changes highlighted
here should be interpreted as hypothesis-generating, guiding further
mechanistic or functional experiments, rather than serving as conclusive
evidence of altered signaling.

Interestingly, although overall EPHB signaling activity was reduced in
ASD according to CellChat’s pathway-level analysis, certain EPHB
receptors appeared selectively expressed in ASD within specific cell
types. This suggests that while the global coordination of EPHB-mediated
interactions may be weakened, isolated ligand or receptor expression may
still occur in ASD, possibly reflecting disrupted or incomplete synaptic
signaling programs. This highlights the importance of interpreting
pathway strength in conjunction with cell-type-specific expression data.

## References

- Cao, J., Spielmann, M., Qiu, X. et al. (2019). The single-cell
  transcriptional landscape of mammalian organogenesis. Nature 566,
  496–502. <https://doi.org/10.1038/s41586-019-0969-x>

- Charman, T., Jones, C. R. G., Pickles, A., Simonoff, E., Baird, G., &
  Happé, F. (2011). Defining the cognitive phenotype of autism. Brain
  Research, 1380, 10–21.
  <https://doi.org/10.1016/j.brainres.2010.10.075>

- Courchesne, E., Pramparo, T., Gazestani, V. H., Lombardo, M. V.,
  Pierce, K., & Lewis, N. E. (2019). The ASD Living Biology: from cell
  proliferation to clinical phenotype. Molecular psychiatry, 24(1),
  88-107.

- Courchesne, E., Gazestani, V. H., & Lewis, N. E. (2020). Prenatal
  origins of ASD: the when, what, and how of ASD development. Trends in
  neurosciences, 43(5), 326-342.

- De Rubeis, S., He, X., Goldberg, A. P., Poultney, C. S., Samocha, K.,
  Ercument Cicek, A., … & Buxbaum, J. D. (2014). Synaptic,
  transcriptional and chromatin genes disrupted in autism. Nature,
  515(7526), 209-215.

- Eichmüller, O. L., & Knoblich, J. A. (2022). Human cerebral
  organoids—a new tool for clinical neurology research. Nature Reviews
  Neurology, 18(11), 661-680.

- Fujita, E., Tanabe, Y., Imhof, B. A., Momoi, M. Y., & Momoi, T.
  (2012). A complex of synaptic adhesion molecule CADM 1, a molecule
  related to autism spectrum disorder, with MUPP 1 in the cerebellum.
  Journal of neurochemistry, 123(5), 886-894.

- Gomez, A. M., Traunmüller, L., & Scheiffele, P. (2021). Neurexins:
  molecular codes for shaping neuronal synapses. Nature Reviews
  Neuroscience, 22(3), 137-151.

- Han, X., He, Y., Wang, Y., Hu, W., Chu, C., Huang, L., Hong, Y., Han,
  L., Zhang, X., Gao, Y., Lin, Y., Ma, H., Shen, H., Ke, X., Liu, Y., &
  Hu, Z. (2025). Deficiency of FABP7 Triggers Premature Neural
  Differentiation in Idiopathic Normocephalic Autism Organoids. Advanced
  Science, 12(2), 2406849. <https://doi.org/10.1002/advs.202406849>

- Imayoshi, I., Sakamoto, M., Yamaguchi, M., Mori, K., & Kageyama, R.
  (2010). Essential roles of Notch signaling in maintenance of neural
  stem cells in developing and adult brains. Journal of Neuroscience,
  30(9), 3489-3498.

- Jin, S., Plikus, M. V., & Nie, Q. (2025). CellChat for systematic
  analysis of cell–cell communication from single-cell transcriptomics.
  Nature Protocols, 20(1), 180-219.

- Lai, K. O., & Ip, N. Y. (2009). Synapse development and plasticity:
  roles of ephrin/Eph receptor signaling. Current opinion in
  neurobiology, 19(3), 275-283.

- Nusse, R., Fuerer, C., Ching, W., Harnish, K., Logan, C., Zeng, A., …
  & Kalani, Y. (2008, January). Wnt signaling and stem cell control. In
  Cold Spring Harbor symposia on quantitative biology (Vol. 73,
  pp. 59-66). Cold Spring Harbor Laboratory Press.

- Pangrazzi, L., Balasco, L., & Bozzi, Y. (2020). Oxidative stress and
  immune system dysfunction in autism spectrum disorders. International
  journal of molecular sciences, 21(9), 3293.

- Pera, E. M., Ikeda, A., Eivers, E., & De Robertis, E. M. (2003).
  Integration of IGF, FGF, and anti-BMP signals via Smad1
  phosphorylation in neural induction. Genes & development, 17(24),
  3023-3028.

- Satterstrom, F. K., Kosmicki, J. A., Wang, J., Breen, M. S., De
  Rubeis, S., An, J. Y., Peng, M., Collins, R., Grove, J., Klei, L.,
  Stevens, C., Reichert, J., Mulhern, M. S., Artomov, M., Gerges, S.,
  Sheppard, B., Xu, X., Bhaduri, A., Norman, U., Brand, H., …
  Buxbaum, J. D. (2020). Large-Scale Exome Sequencing Study Implicates
  Both Developmental and Functional Changes in the Neurobiology of
  Autism. Cell, 180(3), 568–584.e23.
  <https://doi.org/10.1016/j.cell.2019.12.036>

- Suthon, S., Perkins, R. S., Bryja, V., Miranda-Carboni, G. A., &
  Krum, S. A. (2021). WNT5B in Physiology and Disease. Frontiers in cell
  and developmental biology, 9, 667581.

- Wang, X. (2020). Pleiotrophin: Activity and mechanism. Advances in
  clinical chemistry, 98, 51-89.

- Willard, S. S., & Koochekpour, S. (2013). Glutamate, glutamate
  receptors, and downstream signaling pathways. International journal of
  biological sciences, 9(9), 948.
