# Erythropoietin-driven dynamic proteome adaptations during erythropoiesis prevent iron-overload in the developing embryo

## Citation
Chakraborty S and Andrieux G et al., *SUBMITTED*

## Abstract
Erythropoietin (Epo) ensures survival and proliferation of colony-forming unit erythroid progenitor cells (CFU E) and their differentiation to hemoglobin-containing mature erythrocytes. A lack of Epo-induced responses causes embryonic lethality, but mechanisms facilitating the dynamic interplay of cellular alterations to the organismal level remain unresolved. By time-resolved transcriptomics and proteomics, we show that Epo induces in CFU E cells a gradual transition from proliferation signature-proteins to proteins indicative for differentiation including heme-synthesis enzymes. In the absence of the Epo-receptor (EpoR) in embryos, we observe a lack of hemoglobin in CFU-E cells and massive iron-overload of the fetal liver pointing to a miscommunication between the liver and placenta. A reduction of iron-sulfur cluster-containing proteins involved in oxidative phosphorylation leads in these embryos to a metabolic shift towards glycolysis in the fetal liver. This unexpected link connecting erythropoiesis with the regulation of iron-homeostasis and metabolic reprogramming suggests, that balancing these interactions is crucial for protection from iron-intoxication and survival.

## Scripts
* protein_prediction.R predicts protome adaptation from mRNA-to-protein ratio
* prediction_limma.R performs differential analysis on predicted proteome
* prediction_gage.R performs gene-set enrichment analysis on predicted proteome
* aracne_bt.R reconstructs gene regulatory network from predicted proteome
* protein_limma.R performs differential analysis on CFU-E, liver and placenta proteomes
* protein_gage.R performs gene-set enrichment analysis on CFU-E, liver and placenta proteomes

