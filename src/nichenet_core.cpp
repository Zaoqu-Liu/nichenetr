// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <numeric>
#include <cmath>
#include <unordered_set>
#include <unordered_map>
#include <string>

using namespace Rcpp;
using namespace arma;

//' Compute ROC and PR Metrics
//' 
//' @description
//' Computes Area Under ROC Curve (AUROC), Area Under Precision-Recall Curve (AUPR),
//' and related classification metrics.
//' 
//' @param target_vector Logical vector indicating positive class membership
//' @param scores Numeric vector of prediction scores
//' 
//' @return List containing auroc, aupr, aupr_corrected, sensitivity_roc, specificity_roc
//' 
//' @details
//' AUROC is computed using trapezoidal integration under the ROC curve.
//' AUPR is computed using trapezoidal integration under the Precision-Recall curve.
//' The corrected AUPR subtracts the random baseline (proportion of positives).
//' Optimal sensitivity/specificity are determined by maximizing Youden's J statistic.
//' 
//' @keywords internal
// [[Rcpp::export]]
List compute_roc_pr_metrics(LogicalVector target_vector, NumericVector scores) {
  int n = scores.size();
  int n_pos = sum(target_vector);
  int n_neg = n - n_pos;
  
  // Handle edge cases
  if (n_pos == 0 || n_neg == 0 || n == 0) {
    return List::create(
      Named("auroc") = NA_REAL,
      Named("aupr") = NA_REAL,
      Named("aupr_corrected") = NA_REAL,
      Named("sensitivity_roc") = NA_REAL,
      Named("specificity_roc") = NA_REAL
    );
  }
  
  // Sort indices by score (descending order)
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(),
            [&scores](int i, int j) { return scores[i] > scores[j]; });
  
  // Compute ROC and PR curves via single pass
  int tp = 0, fp = 0;
  double auroc = 0.0;
  double aupr = 0.0;
  double prev_tpr = 0.0;
  double prev_fpr = 0.0;
  double prev_recall = 0.0;
  
  for (int i = 0; i < n; i++) {
    if (target_vector[idx[i]]) {
      tp++;
    } else {
      fp++;
    }
    
    double tpr = static_cast<double>(tp) / n_pos;  // True Positive Rate = Recall
    double fpr = static_cast<double>(fp) / n_neg;  // False Positive Rate
    double precision = static_cast<double>(tp) / (i + 1);
    
    // AUROC: trapezoidal rule integration
    auroc += (fpr - prev_fpr) * ((tpr + prev_tpr) / 2.0);
    
    // AUPR: rectangular approximation (standard for PR curves)
    aupr += (tpr - prev_recall) * precision;
    
    prev_tpr = tpr;
    prev_fpr = fpr;
    prev_recall = tpr;
  }
  
  // Baseline-corrected AUPR
  double random_aupr = static_cast<double>(n_pos) / n;
  double aupr_corrected = aupr - random_aupr;
  
  // Find optimal threshold via Youden's J statistic
  double best_sens = 0.5, best_spec = 0.5;
  double max_youden = -1.0;
  tp = 0; fp = 0;
  
  for (int i = 0; i < n; i++) {
    if (target_vector[idx[i]]) tp++; else fp++;
    double sens = static_cast<double>(tp) / n_pos;
    double spec = 1.0 - static_cast<double>(fp) / n_neg;
    double youden = sens + spec - 1.0;
    
    if (youden > max_youden) {
      max_youden = youden;
      best_sens = sens;
      best_spec = spec;
    }
  }
  
  return List::create(
    Named("auroc") = auroc,
    Named("aupr") = aupr,
    Named("aupr_corrected") = aupr_corrected,
    Named("sensitivity_roc") = best_sens,
    Named("specificity_roc") = best_spec
  );
}

//' Predict Ligand Activities (C++ Implementation)
//' 
//' @description
//' Computes ligand activity metrics for all ligands simultaneously.
//' 
//' @param target_genes Character vector of target gene symbols
//' @param gene_names Character vector of all gene names (matrix rownames)
//' @param ligand_target_matrix Numeric matrix of ligand-target scores (genes x ligands)
//' @param ligand_names Character vector of ligand names (matrix colnames)
//' 
//' @return DataFrame with columns: test_ligand, auroc, aupr, aupr_corrected, pearson
//' 
//' @details
//' For each ligand, computes:
//' - AUROC: Area under ROC curve for target gene prediction
//' - AUPR: Area under Precision-Recall curve
//' - AUPR corrected: AUPR minus random baseline
//' - Pearson: Correlation between ligand-target scores and target membership
//' 
//' @keywords internal
// [[Rcpp::export]]
DataFrame predict_ligand_activities_cpp(CharacterVector target_genes,
                                        CharacterVector gene_names,
                                        NumericMatrix ligand_target_matrix,
                                        CharacterVector ligand_names) {
  
  int n_ligands = ligand_target_matrix.ncol();
  int n_genes = ligand_target_matrix.nrow();
  
  if (n_ligands == 0 || n_genes == 0) {
    return DataFrame::create();
  }
  
  // Build target gene set for O(1) lookup
  std::unordered_set<std::string> target_set;
  for (int j = 0; j < target_genes.size(); j++) {
    target_set.insert(as<std::string>(target_genes[j]));
  }
  
  // Create binary target indicator vector
  LogicalVector is_target(n_genes, false);
  for (int i = 0; i < n_genes; i++) {
    std::string gene_i = as<std::string>(gene_names[i]);
    if (target_set.find(gene_i) != target_set.end()) {
      is_target[i] = true;
    }
  }
  
  // Precompute target statistics for Pearson correlation
  int n_pos = sum(is_target);
  double target_mean = static_cast<double>(n_pos) / n_genes;
  double target_ss = 0.0;
  for (int i = 0; i < n_genes; i++) {
    double diff = (is_target[i] ? 1.0 : 0.0) - target_mean;
    target_ss += diff * diff;
  }
  double target_sd = std::sqrt(target_ss / n_genes);
  
  // Result vectors
  NumericVector auroc_vals(n_ligands);
  NumericVector aupr_vals(n_ligands);
  NumericVector aupr_corrected_vals(n_ligands);
  NumericVector pearson_vals(n_ligands);
  
  // Process each ligand
  for (int lig = 0; lig < n_ligands; lig++) {
    NumericVector scores = ligand_target_matrix(_, lig);
    
    // Compute ROC/PR metrics
    List metrics = compute_roc_pr_metrics(is_target, scores);
    auroc_vals[lig] = metrics["auroc"];
    aupr_vals[lig] = metrics["aupr"];
    aupr_corrected_vals[lig] = metrics["aupr_corrected"];
    
    // Compute Pearson correlation
    double score_mean = mean(scores);
    double score_ss = 0.0;
    double covariance = 0.0;
    
    for (int i = 0; i < n_genes; i++) {
      double score_diff = scores[i] - score_mean;
      double target_val = is_target[i] ? 1.0 : 0.0;
      double target_diff = target_val - target_mean;
      
      score_ss += score_diff * score_diff;
      covariance += score_diff * target_diff;
    }
    
    double score_sd = std::sqrt(score_ss / n_genes);
    
    // Numerical stability check
    if (score_sd > 1e-10 && target_sd > 1e-10) {
      pearson_vals[lig] = covariance / (n_genes * score_sd * target_sd);
    } else {
      pearson_vals[lig] = NA_REAL;
    }
  }
  
  return DataFrame::create(
    Named("test_ligand") = ligand_names,
    Named("auroc") = auroc_vals,
    Named("aupr") = aupr_vals,
    Named("aupr_corrected") = aupr_corrected_vals,
    Named("pearson") = pearson_vals,
    Named("stringsAsFactors") = false
  );
}

//' Compute Pearson Correlations for Multiple Ligands
//' 
//' @param target_vector Logical vector for target gene membership
//' @param ligand_target_matrix Matrix of ligand-target scores
//' 
//' @return Numeric vector of Pearson correlations
//' 
//' @keywords internal
// [[Rcpp::export]]
NumericVector compute_pearson_batch(LogicalVector target_vector, 
                                    NumericMatrix ligand_target_matrix) {
  int n_ligands = ligand_target_matrix.ncol();
  int n_genes = ligand_target_matrix.nrow();
  NumericVector correlations(n_ligands);
  
  // Convert logical to numeric for correlation
  NumericVector target_numeric(n_genes);
  for (int i = 0; i < n_genes; i++) {
    target_numeric[i] = target_vector[i] ? 1.0 : 0.0;
  }
  
  double mean_target = mean(target_numeric);
  
  // Precompute target variance
  double var_target = 0.0;
  for (int i = 0; i < n_genes; i++) {
    double diff = target_numeric[i] - mean_target;
    var_target += diff * diff;
  }
  double sd_target = std::sqrt(var_target / n_genes);
  
  if (sd_target < 1e-10) {
    std::fill(correlations.begin(), correlations.end(), NA_REAL);
    return correlations;
  }
  
  // Compute correlation for each ligand
  for (int lig = 0; lig < n_ligands; lig++) {
    NumericVector scores = ligand_target_matrix(_, lig);
    double mean_scores = mean(scores);
    
    double var_scores = 0.0;
    double covar = 0.0;
    
    for (int i = 0; i < n_genes; i++) {
      double diff_target = target_numeric[i] - mean_target;
      double diff_scores = scores[i] - mean_scores;
      var_scores += diff_scores * diff_scores;
      covar += diff_target * diff_scores;
    }
    
    double sd_scores = std::sqrt(var_scores / n_genes);
    
    if (sd_scores < 1e-10) {
      correlations[lig] = NA_REAL;
    } else {
      correlations[lig] = covar / (n_genes * sd_target * sd_scores);
    }
  }
  
  return correlations;
}

//' Extract Top Target Genes for a Ligand
//' 
//' @param ligand_idx Integer index of the ligand (0-based)
//' @param ligand_target_matrix Matrix of ligand-target scores
//' @param target_genes Character vector of genes in the target set
//' @param gene_names Character vector of all gene names
//' @param ligand_names Character vector of ligand names
//' @param top_n Number of top targets to return
//' 
//' @return DataFrame with ligand-target links
//' 
//' @keywords internal
// [[Rcpp::export]]
DataFrame get_top_targets_cpp(int ligand_idx,
                              NumericMatrix ligand_target_matrix,
                              CharacterVector target_genes,
                              CharacterVector gene_names,
                              CharacterVector ligand_names,
                              int top_n) {
  
  int n_genes = ligand_target_matrix.nrow();
  NumericVector scores = ligand_target_matrix(_, ligand_idx);
  
  // Build hash set of target genes
  std::unordered_set<std::string> target_set;
  for (int j = 0; j < target_genes.size(); j++) {
    target_set.insert(as<std::string>(target_genes[j]));
  }
  
  // Filter to genes in target set
  std::vector<int> geneset_indices;
  std::vector<double> geneset_scores;
  
  for (int i = 0; i < n_genes; i++) {
    std::string gene = as<std::string>(gene_names[i]);
    if (target_set.find(gene) != target_set.end()) {
      geneset_indices.push_back(i);
      geneset_scores.push_back(scores[i]);
    }
  }
  
  if (geneset_indices.empty()) {
    return DataFrame::create(
      Named("ligand") = CharacterVector::create(),
      Named("target") = CharacterVector::create(),
      Named("weight") = NumericVector::create()
    );
  }
  
  // Sort by score (descending)
  std::vector<int> sorted_idx(geneset_indices.size());
  std::iota(sorted_idx.begin(), sorted_idx.end(), 0);
  std::sort(sorted_idx.begin(), sorted_idx.end(),
            [&geneset_scores](int i, int j) { 
              return geneset_scores[i] > geneset_scores[j]; 
            });
  
  int actual_n = std::min(top_n, static_cast<int>(sorted_idx.size()));
  CharacterVector target_names(actual_n);
  NumericVector weights(actual_n);
  CharacterVector ligand_rep(actual_n);
  
  std::string ligand_name = as<std::string>(ligand_names[ligand_idx]);
  for (int i = 0; i < actual_n; i++) {
    int orig_idx = geneset_indices[sorted_idx[i]];
    target_names[i] = gene_names[orig_idx];
    weights[i] = geneset_scores[sorted_idx[i]];
    ligand_rep[i] = ligand_name;
  }
  
  return DataFrame::create(
    Named("ligand") = ligand_rep,
    Named("target") = target_names,
    Named("weight") = weights,
    Named("stringsAsFactors") = false
  );
}

//' Subset Matrix by Gene Names
//' 
//' @param matrix Input matrix
//' @param row_names Current row names
//' @param genes_to_keep Genes to retain
//' 
//' @return Subsetted matrix
//' 
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix subset_matrix_by_genes(NumericMatrix matrix,
                                     CharacterVector row_names,
                                     CharacterVector genes_to_keep) {
  
  int n_rows = matrix.nrow();
  int n_cols = matrix.ncol();
  
  // Build hash set for O(1) lookup
  std::unordered_set<std::string> keep_set;
  for (int j = 0; j < genes_to_keep.size(); j++) {
    keep_set.insert(as<std::string>(genes_to_keep[j]));
  }
  
  // Find indices to keep
  std::vector<int> keep_indices;
  CharacterVector new_row_names;
  
  for (int i = 0; i < n_rows; i++) {
    std::string gene = as<std::string>(row_names[i]);
    if (keep_set.find(gene) != keep_set.end()) {
      keep_indices.push_back(i);
      new_row_names.push_back(row_names[i]);
    }
  }
  
  // Create subsetted matrix
  int new_n_rows = keep_indices.size();
  NumericMatrix result(new_n_rows, n_cols);
  
  for (int i = 0; i < new_n_rows; i++) {
    result(i, _) = matrix(keep_indices[i], _);
  }
  
  rownames(result) = new_row_names;
  colnames(result) = colnames(matrix);
  
  return result;
}

//' Complete Ligand Activity Prediction Pipeline
//' 
//' @description
//' Performs the complete ligand activity prediction workflow:
//' gene intersections, matrix subsetting, and metric computation.
//' 
//' @param geneset Character vector of target genes
//' @param background_genes Character vector of expressed genes
//' @param ligand_target_matrix Ligand-target matrix (genes x ligands)
//' @param ligand_names All ligand names
//' @param potential_ligands Ligands to evaluate
//' 
//' @return DataFrame with ligand activity results sorted by aupr_corrected
//' 
//' @keywords internal
// [[Rcpp::export]]
DataFrame predict_ligand_activities_direct(
    CharacterVector geneset,
    CharacterVector background_genes,
    NumericMatrix ligand_target_matrix,
    CharacterVector ligand_names,
    CharacterVector potential_ligands) {
  
  int n_genes_full = ligand_target_matrix.nrow();
  int n_ligands_full = ligand_target_matrix.ncol();
  
  // Build gene sets
  std::unordered_set<std::string> background_set;
  for (int i = 0; i < background_genes.size(); i++) {
    background_set.insert(as<std::string>(background_genes[i]));
  }
  
  std::unordered_set<std::string> geneset_set;
  for (int i = 0; i < geneset.size(); i++) {
    geneset_set.insert(as<std::string>(geneset[i]));
  }
  
  std::unordered_set<std::string> potential_set;
  for (int i = 0; i < potential_ligands.size(); i++) {
    potential_set.insert(as<std::string>(potential_ligands[i]));
  }
  
  // Find valid gene indices (intersection with background)
  std::vector<int> gene_indices;
  std::vector<std::string> valid_genes;
  CharacterVector matrix_rownames = rownames(ligand_target_matrix);
  
  for (int i = 0; i < n_genes_full; i++) {
    std::string gene = as<std::string>(matrix_rownames[i]);
    if (background_set.find(gene) != background_set.end()) {
      gene_indices.push_back(i);
      valid_genes.push_back(gene);
    }
  }
  
  if (gene_indices.empty()) {
    return DataFrame::create();
  }
  
  // Find valid ligand indices
  std::vector<int> ligand_indices;
  std::vector<std::string> valid_ligands;
  
  for (int j = 0; j < n_ligands_full; j++) {
    std::string ligand = as<std::string>(ligand_names[j]);
    if (potential_set.find(ligand) != potential_set.end()) {
      ligand_indices.push_back(j);
      valid_ligands.push_back(ligand);
    }
  }
  
  if (ligand_indices.empty()) {
    return DataFrame::create();
  }
  
  int n_genes = gene_indices.size();
  int n_ligands = ligand_indices.size();
  
  // Create subsetted matrix
  arma::mat mat_subset(n_genes, n_ligands);
  for (int i = 0; i < n_genes; i++) {
    for (int j = 0; j < n_ligands; j++) {
      mat_subset(i, j) = ligand_target_matrix(gene_indices[i], ligand_indices[j]);
    }
  }
  
  // Create binary target vector
  LogicalVector is_target(n_genes);
  for (int i = 0; i < n_genes; i++) {
    is_target[i] = (geneset_set.find(valid_genes[i]) != geneset_set.end());
  }
  
  int n_pos = sum(is_target);
  if (n_pos == 0) {
    return DataFrame::create();
  }
  
  // Precompute target statistics
  double target_mean = static_cast<double>(n_pos) / n_genes;
  double target_var = 0.0;
  for (int i = 0; i < n_genes; i++) {
    double diff = (is_target[i] ? 1.0 : 0.0) - target_mean;
    target_var += diff * diff;
  }
  double target_sd = std::sqrt(target_var / n_genes);
  
  // Compute metrics for all ligands
  NumericVector auroc_vals(n_ligands);
  NumericVector aupr_vals(n_ligands);
  NumericVector aupr_corrected_vals(n_ligands);
  NumericVector pearson_vals(n_ligands);
  
  for (int lig = 0; lig < n_ligands; lig++) {
    arma::vec scores = mat_subset.col(lig);
    
    // Sort indices by score (descending)
    arma::uvec sorted_idx = arma::sort_index(scores, "descend");
    
    // Compute AUROC and AUPR
    int tp = 0, fp = 0;
    double auroc = 0.0, aupr = 0.0;
    double prev_tpr = 0.0, prev_fpr = 0.0, prev_recall = 0.0;
    
    for (size_t i = 0; i < sorted_idx.n_elem; i++) {
      if (is_target[sorted_idx[i]]) {
        tp++;
      } else {
        fp++;
      }
      
      double tpr = static_cast<double>(tp) / n_pos;
      double fpr = static_cast<double>(fp) / (n_genes - n_pos);
      double precision = static_cast<double>(tp) / (i + 1);
      
      auroc += (fpr - prev_fpr) * ((tpr + prev_tpr) / 2.0);
      aupr += (tpr - prev_recall) * precision;
      
      prev_tpr = tpr;
      prev_fpr = fpr;
      prev_recall = tpr;
    }
    
    auroc_vals[lig] = auroc;
    aupr_vals[lig] = aupr;
    aupr_corrected_vals[lig] = aupr - target_mean;
    
    // Compute Pearson correlation
    double score_mean = arma::mean(scores);
    double score_var = 0.0, covar = 0.0;
    
    for (int i = 0; i < n_genes; i++) {
      double score_diff = scores[i] - score_mean;
      double target_val = is_target[i] ? 1.0 : 0.0;
      double target_diff = target_val - target_mean;
      
      score_var += score_diff * score_diff;
      covar += score_diff * target_diff;
    }
    
    double score_sd = std::sqrt(score_var / n_genes);
    
    if (score_sd > 1e-10 && target_sd > 1e-10) {
      pearson_vals[lig] = covar / (n_genes * score_sd * target_sd);
    } else {
      pearson_vals[lig] = NA_REAL;
    }
  }
  
  // Create result DataFrame
  CharacterVector ligand_names_result(n_ligands);
  for (int i = 0; i < n_ligands; i++) {
    ligand_names_result[i] = valid_ligands[i];
  }
  
  // Sort by aupr_corrected descending
  arma::uvec sort_idx = arma::sort_index(as<arma::vec>(aupr_corrected_vals), "descend");
  
  CharacterVector sorted_ligands(n_ligands);
  NumericVector sorted_auroc(n_ligands);
  NumericVector sorted_aupr(n_ligands);
  NumericVector sorted_aupr_corr(n_ligands);
  NumericVector sorted_pearson(n_ligands);
  
  for (int i = 0; i < n_ligands; i++) {
    int idx = sort_idx[i];
    sorted_ligands[i] = ligand_names_result[idx];
    sorted_auroc[i] = auroc_vals[idx];
    sorted_aupr[i] = aupr_vals[idx];
    sorted_aupr_corr[i] = aupr_corrected_vals[idx];
    sorted_pearson[i] = pearson_vals[idx];
  }
  
  return DataFrame::create(
    Named("test_ligand") = sorted_ligands,
    Named("auroc") = sorted_auroc,
    Named("aupr") = sorted_aupr,
    Named("aupr_corrected") = sorted_aupr_corr,
    Named("pearson") = sorted_pearson,
    Named("stringsAsFactors") = false
  );
}
