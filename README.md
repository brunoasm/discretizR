# discretizR

## Purpose
This program is an interactive tool to assist in discretizing measurements for phylogenetic analyses.

It helps in removing effects of size, separating dimorphic characters when the dimorphism is not constant betwen species, checking for character correlation and discretizing characters based on Gaussian mixtures.

## Input

To run the program, we need two kinds of input file (see the `test_data` folder on github for formatting):

1. measurement table

A table with measurements, in csv format. It needs at least one column named species and one column for each measurement. If your species have morphs (e. g. sexes), you need a column with the morphs as well (you can choose whatever column title you want).


2. Character table

A table in csv format with the characters that you want to extract from the measurements table. It needs three columns, with the follwing names:

* primary_variable
* size_variable
* morph_variable

In each row, you will list which variable (i. e. column in the measurement table) is the one you want to use for this character. This is the primary_variable.

If you want to remove the effect of size, you can also provide a size_variable. If you want to calculate the size as the geometric mean of more than one variable in the measurement table, you can enclose all variables with quotes and separate them with commas, as in "scutellum_width,scutellum_length"

Finally, if you want to account for morph, you can provide a morph_variable, which should contain the name of a column in the measurement table with morph information

##Tabs

## Left panel

This panel contains the main input/output functions. One can upload tables, select variables to use for body size, download final matrices or select characters to be visualized in regression plots.

### Regression plots

The graph on top plots the selected variable against corresponding size variable. In the bottom, regression residuals for each species are shown. The use has the option to color points by species, character state, or morph.

If a species is associated with more than one character state, it will contain a polymorphism in the output.

### Gaussian misxture diagnostics

After fitting regressions, discretizR uses the mclust package to fit mixtures of models to the distribution of residuals. The number of components is decided automatically by using the Bayesian Information Criterion. See help of the package mclust (https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html) for more information on the diagnostic plots shown in this screen

### Character correlation

Here the user can select a minimum Pearson correlation coefficient to cluster characters together. After selecting a minimum level, characters will be clustered in groups of correlated characters, and the user can use the selection box to the right to navigate through these groups. 

For each group, discretizR shows all pairwise Pearson correlation coefficients and the associated p values. This is meant to be used as an exploratory tool to help select characters for inclusion. If one decides to remove certain characters from the output, one needs to select these characters and then click on **EXCLUDE SELECTED FROM OUTPUT**.

If no character is selected for exclusion, all characters will be present in the output regardless of their correlation group.

##Output

After selecting characters for exclusion, the user can download character matrices in the formats NEXUS (discretized only) and TNT (continuous or discretized). This is done by clicking in **Download Nexus** or **Download TNT**. For discrete character matrices, one can select in the dialogue whether to include invariant characters and whether to treat characters as ordered. For the continuous TNT matrix, One can select the range to rescale a character.

In both cases, a fully annotated matrix with character names and states explaining their definition is downloaded.

##FAQ

### Why Gaussian mixtures?

If characters are approximatelly normally distributed within species, they will become a mixture of normals as species diverge. This paper used this fact to delimit species, for example:
```Zapata, F. & Jiménez, I. (2012) Species delimitation: Inferring gaps in morphology across geography. Systematic Biology 61, 179–194. https://doi.org/10.1093/sysbio/syr084```

Here we use the R package [mclust] (https://github.com/cran/mclust) to test the distribution of a character (i. e. the residuals explained above) for a mixture of normals with single variance. If the data can be explained by a single normal distribution, it cannot be clustered into discrete characters and is removed from analysis.

### How to test for sexual dimorphism and other polymorphisms?

If sexual dimorphism is constant across species (for example, males are always about 1.5 times larger), we can account for that simply by including morph (i. e. sex) in the regression to remove size information. 

As a first step, discretizR fits a regression with size, morph and species as predictors. If we find a significant interaction between species and morph, it means that sexual dimorphism does not have a constant effect across species.

In that is the case, we do not include morph as a predictor to calculate residuals. Additionally, we create one character for each morph (e. g. male and female) in the output. The user can later decide whether to include only one or both characters to infer a phylogeny.

### How to remove character correlation?
Characters are correlated primarily due to phylogeny itself, but there are other sources of correlation (for example, genetic correlations). Many of these other sources of correlation should be found within species, so a quick test to check for correlated characters is to test for correlation (here, Pearson correlation) in the residuals of a regression in which species is one of the predictors (size and morph are also included, if there is information in the character table). 

Of course, this can only be accomplished if there are multiple individuals measured per species (i. e. degress of freedom). Here we only attempt this method if there are more than 10 degrees of freedom in the regression residuals, which roughly means at least 10 independent observations across species. This is a somewhat low threshold, the more the better. In any case, if your data does not include multiple measurements per species, I do not think there is a decent way to remove correlation between characters that does not remove phylogenetic signal.

### Can I upload the nexus file to MorphoBank?

Yes, the output is fully compatible with MorphoBank!
