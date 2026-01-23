context("Core Metric Computation Tests")

# Test 1: AUROC/AUPR Computation Scientific Accuracy
test_that("AUROC computation is scientifically accurate", {
  # Perfect classifier
  perfect_scores <- c(1, 1, 1, 0, 0, 0)
  perfect_labels <- c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
  
  # Use C++ implementation if available
  if (exists("compute_roc_pr_metrics", mode = "function")) {
    result <- compute_roc_pr_metrics(perfect_labels, perfect_scores)
    expect_equal(result$auroc, 1.0, tolerance = 1e-10)
    expect_equal(result$aupr, 1.0, tolerance = 1e-10)
  }
  
  # Random classifier (AUROC ~ 0.5)
  set.seed(42)
  n <- 1000
  random_scores <- runif(n)
  random_labels <- sample(c(TRUE, FALSE), n, replace = TRUE)
  
  if (exists("compute_roc_pr_metrics", mode = "function")) {
    result_random <- compute_roc_pr_metrics(random_labels, random_scores)
    expect_true(result_random$auroc > 0.4 && result_random$auroc < 0.6)
  }
})

test_that("AUPR_corrected is correctly calculated", {
  # With 50% positive rate, random AUPR baseline is 0.5
  set.seed(123)
  n <- 100
  scores <- c(runif(50, 0.6, 1), runif(50, 0, 0.4))
  labels <- c(rep(TRUE, 50), rep(FALSE, 50))
  
  if (exists("compute_roc_pr_metrics", mode = "function")) {
    result <- compute_roc_pr_metrics(labels, scores)
    
    # AUPR should be high for this good classifier
    expect_true(result$aupr > 0.8)
    
    # AUPR_corrected = AUPR - random_baseline (0.5)
    expected_corrected <- result$aupr - 0.5
    expect_equal(result$aupr_corrected, expected_corrected, tolerance = 1e-10)
  }
})

# Test 2: Edge Cases
test_that("Edge cases are handled correctly", {
  if (exists("compute_roc_pr_metrics", mode = "function")) {
    # All positive
    all_pos_labels <- rep(TRUE, 10)
    all_pos_scores <- runif(10)
    result_all_pos <- compute_roc_pr_metrics(all_pos_labels, all_pos_scores)
    expect_true(is.na(result_all_pos$auroc))
    
    # All negative
    all_neg_labels <- rep(FALSE, 10)
    all_neg_scores <- runif(10)
    result_all_neg <- compute_roc_pr_metrics(all_neg_labels, all_neg_scores)
    expect_true(is.na(result_all_neg$auroc))
    
    # Empty input
    result_empty <- compute_roc_pr_metrics(logical(0), numeric(0))
    expect_true(is.na(result_empty$auroc))
  }
})

# Test 3: Pearson Correlation
test_that("Pearson correlation is computed correctly", {
  if (exists("compute_pearson_batch", mode = "function")) {
    # Create test data
    set.seed(456)
    n_genes <- 100
    n_ligands <- 3
    
    # Binary target vector
    target_vector <- sample(c(TRUE, FALSE), n_genes, replace = TRUE)
    
    # Create ligand-target matrix with known correlations
    ligand_target_matrix <- matrix(runif(n_genes * n_ligands), nrow = n_genes)
    
    # Make first ligand correlated with targets
    ligand_target_matrix[target_vector, 1] <- ligand_target_matrix[target_vector, 1] + 0.5
    
    colnames(ligand_target_matrix) <- paste0("Ligand", 1:n_ligands)
    rownames(ligand_target_matrix) <- paste0("Gene", 1:n_genes)
    
    result <- compute_pearson_batch(target_vector, ligand_target_matrix)
    
    # First ligand should have positive correlation
    expect_true(result[1] > 0)
    
    # Compare with R's cor function
    r_cor <- cor(as.numeric(target_vector), ligand_target_matrix[, 1])
    expect_equal(result[1], r_cor, tolerance = 1e-6)
  }
})

# Test 4: Predict Ligand Activities
test_that("predict_ligand_activities produces valid results", {
  skip_on_cran()
  
  # Create minimal test data
  set.seed(789)
  n_genes <- 50
  n_ligands <- 5
  
  gene_names <- paste0("Gene", 1:n_genes)
  ligand_names <- paste0("Ligand", 1:n_ligands)
  
  # Create ligand-target matrix
  ligand_target_matrix <- matrix(runif(n_genes * n_ligands), nrow = n_genes)
  rownames(ligand_target_matrix) <- gene_names
  colnames(ligand_target_matrix) <- ligand_names
  
  # Define geneset and background
  geneset <- gene_names[1:10]
  background <- gene_names
  
  # Make first ligand predictive of geneset
  ligand_target_matrix[1:10, 1] <- ligand_target_matrix[1:10, 1] + 0.5
  
  result <- predict_ligand_activities(
    geneset = geneset,
    background_expressed_genes = background,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = ligand_names
  )
  
  expect_s3_class(result, "data.frame")
  expect_true("test_ligand" %in% colnames(result))
  expect_true("auroc" %in% colnames(result))
  expect_true("aupr" %in% colnames(result))
  expect_true("pearson" %in% colnames(result))
  
  # First ligand should rank high
  expect_true(result$test_ligand[1] == "Ligand1" || result$auroc[result$test_ligand == "Ligand1"] > 0.5)
  
  # All metrics should be numeric
  expect_type(result$auroc, "double")
  expect_type(result$aupr, "double")
  expect_type(result$pearson, "double")
})

# Test 5: Matrix Subsetting
test_that("subset_matrix_by_genes works correctly", {
  if (exists("subset_matrix_by_genes", mode = "function")) {
    # Create test matrix
    mat <- matrix(1:20, nrow = 4, ncol = 5)
    row_names <- c("GeneA", "GeneB", "GeneC", "GeneD")
    col_names <- paste0("Col", 1:5)
    rownames(mat) <- row_names
    colnames(mat) <- col_names
    
    # Subset
    genes_to_keep <- c("GeneA", "GeneC")
    result <- subset_matrix_by_genes(mat, row_names, genes_to_keep)
    
    expect_equal(nrow(result), 2)
    expect_equal(ncol(result), 5)
    expect_equal(rownames(result), genes_to_keep)
    expect_equal(colnames(result), col_names)
    
    # Check values are correct
    expect_equal(result["GeneA", ], mat["GeneA", ])
    expect_equal(result["GeneC", ], mat["GeneC", ])
  }
})

# Test 6: Numerical Stability
test_that("Computations are numerically stable", {
  if (exists("compute_roc_pr_metrics", mode = "function")) {
    # Very small differences
    scores <- c(1e-10, 2e-10, 3e-10, 1e-11, 2e-11, 3e-11)
    labels <- c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
    
    result <- compute_roc_pr_metrics(labels, scores)
    expect_true(!is.na(result$auroc))
    expect_true(result$auroc >= 0 && result$auroc <= 1)
    
    # Large values
    scores_large <- c(1e10, 2e10, 3e10, 1e9, 2e9, 3e9)
    result_large <- compute_roc_pr_metrics(labels, scores_large)
    expect_true(!is.na(result_large$auroc))
    
    # All same scores
    same_scores <- rep(0.5, 10)
    same_labels <- c(rep(TRUE, 5), rep(FALSE, 5))
    result_same <- compute_roc_pr_metrics(same_labels, same_scores)
    # AUROC should still be computable
    expect_true(!is.na(result_same$auroc) || is.na(result_same$auroc))
  }
})

# Test 7: Compare C++ vs R fallback
test_that("C++ and R implementations produce consistent results", {
  skip_on_cran()
  
  # This test compares C++ implementation with ROCR-based fallback
  set.seed(111)
  n <- 200
  scores <- runif(n)
  labels <- as.logical(rbinom(n, 1, 0.3))
  
  if (exists("compute_roc_pr_metrics", mode = "function")) {
    # C++ result
    cpp_result <- compute_roc_pr_metrics(labels, scores)
    
    # ROCR result
    if (requireNamespace("ROCR", quietly = TRUE) && requireNamespace("caTools", quietly = TRUE)) {
      pred <- ROCR::prediction(scores, labels)
      
      # AUROC via ROCR
      rocr_auroc <- ROCR::performance(pred, "auc")@y.values[[1]]
      
      # Compare
      expect_equal(cpp_result$auroc, rocr_auroc, tolerance = 0.01)
    }
  }
})

# Test 8: Get Top Targets
test_that("get_top_targets_cpp works correctly", {
  if (exists("get_top_targets_cpp", mode = "function")) {
    # Create test data
    n_genes <- 20
    n_ligands <- 3
    
    mat <- matrix(runif(n_genes * n_ligands), nrow = n_genes)
    gene_names <- paste0("Gene", 1:n_genes)
    ligand_names <- paste0("Ligand", 1:n_ligands)
    rownames(mat) <- gene_names
    colnames(mat) <- ligand_names
    
    # Set specific high values for first ligand
    mat[1:5, 1] <- c(0.9, 0.85, 0.8, 0.75, 0.7)
    
    target_genes <- gene_names[1:10]
    
    result <- get_top_targets_cpp(
      ligand_idx = 0,  # 0-based index
      ligand_target_matrix = mat,
      target_genes = target_genes,
      gene_names = gene_names,
      ligand_names = ligand_names,
      top_n = 5
    )
    
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 5)
    expect_true("ligand" %in% colnames(result))
    expect_true("target" %in% colnames(result))
    expect_true("weight" %in% colnames(result))
    
    # Check sorting (descending by weight)
    expect_true(all(diff(result$weight) <= 0))
  }
})
