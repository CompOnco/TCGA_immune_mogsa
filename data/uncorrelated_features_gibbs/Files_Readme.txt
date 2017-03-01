

Feature selection on gene sets
Gibbs, Feb 20 2017

In order to find a minimally sized set of uncorrelated gene-set-features
for downstream clustering, I have performed a search over sets, that minimizes the
within-set correlation. This is done in a repeated mcmc style of algorithm.

The result files in each folder are:

adaptive_feature_selection.R:
  The script used to generate the results. Parameters are found at the top of
  the file. The most important being the learning rate and repetitions.
  For each repetition, a probability distribution is generated over gene sets.
  In each iteration, we find the within-set correlation, and use that to
  update the probability distribution, scaled by the learning rate. Finally,
  after each repetition that span a number of iterations, we have a list of
  probability distributions. We do one final sampling from each distribution,
  and count how often each set has been selected. This is performed on two
  splits of the data, so results can be compared.


run_solutions_and_parametes.rda:
  The rdata file containing results and run parameters. This includes the
  indices for the data split (idx1, idx2), the table of counts from each
  data split, and the probability distributions for each repetition
  (soln1, soln2 list item Ps).


selected_features_barplot_soln1.pdf
  Barplot showing the number of times a feature was selected across reps.


solution_comparison.pdf
  A plot showing the counts for each feature, between the two data splits.

Soution_Counts.txt
  The table of feature counts between the two data splits, sorted by avg. counts.

  
