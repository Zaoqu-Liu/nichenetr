#' Database Loading and Caching System
#'
#' @description
#' Efficient I/O system for NicheNet databases with optional caching.
#'
#' @keywords internal

# Global cache environment
.nichenet_db_cache <- new.env(parent = emptyenv())

#' Load Ligand-Target Matrix with Optional Caching
#'
#' @param file_path Path to RDS or qs file
#' @param use_qs Use qs format if available (default: TRUE)
#' @param cache Store in memory cache (default: TRUE)
#' @param cache_id Custom cache identifier
#' 
#' @return Ligand-target matrix
#' 
#' @keywords internal
load_ligand_target_matrix <- function(file_path, use_qs = TRUE, cache = TRUE, cache_id = NULL) {
  if (is.null(cache_id)) {
    cache_id <- digest::digest(file_path, algo = "xxhash32")
  }

  cache_key <- paste0("ltm_", cache_id)

  # Return cached version if available

if (cache && exists(cache_key, envir = .nichenet_db_cache)) {
    if (getOption("nichenetr.verbose", TRUE)) {
      message("Using cached matrix (", cache_id, ")")
    }
    return(.nichenet_db_cache[[cache_key]])
  }

  # Load from file
  if (use_qs && requireNamespace("qs", quietly = TRUE)) {
    qs_path <- sub("\\.rds$", ".qs", file_path)

    if (file.exists(qs_path)) {
      if (getOption("nichenetr.verbose", TRUE)) {
        message("Loading from qs format...")
      }
      start_time <- Sys.time()
      matrix_data <- qs::qread(qs_path, nthreads = parallel::detectCores())
      load_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      if (getOption("nichenetr.verbose", TRUE)) {
        message(sprintf("  Loaded in %.2f seconds", load_time))
      }
    } else if (file.exists(file_path)) {
      if (getOption("nichenetr.verbose", TRUE)) {
        message("Loading from RDS (will convert to qs for future use)...")
      }
      start_time <- Sys.time()
      matrix_data <- readRDS(file_path)
      load_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      if (getOption("nichenetr.verbose", TRUE)) {
        message(sprintf("  Loaded in %.2f seconds", load_time))
      }

      # Save as qs format for future loading
      tryCatch(
        {
          qs::qsave(matrix_data, qs_path, preset = "fast", nthreads = parallel::detectCores())
          if (getOption("nichenetr.verbose", TRUE)) {
            message("  Saved qs version for future use")
          }
        },
        error = function(e) {
          warning("Could not save qs version: ", e$message)
        }
      )
    } else {
      stop("File not found: ", file_path)
    }
  } else {
    # Fallback to RDS
    if (getOption("nichenetr.verbose", TRUE)) {
      message("Loading from RDS...")
    }
    start_time <- Sys.time()
    matrix_data <- readRDS(file_path)
    load_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    if (getOption("nichenetr.verbose", TRUE)) {
      message(sprintf("  Loaded in %.2f seconds", load_time))
    }
  }

  # Cache if requested
  if (cache) {
    .nichenet_db_cache[[cache_key]] <- matrix_data
    if (getOption("nichenetr.verbose", TRUE)) {
      message(sprintf("  Cached (%.1f MB)", object.size(matrix_data) / 1024^2))
    }
  }

  return(matrix_data)
}

#' Preload NicheNet Databases
#'
#' @param lr_network_path Path to ligand-receptor network
#' @param ligand_target_matrix_path Path to ligand-target matrix
#' @param weighted_networks_path Path to weighted networks (optional)
#' 
#' @return List of loaded databases
#' 
#' @keywords internal
preload_nichenet_databases <- function(lr_network_path,
                                       ligand_target_matrix_path,
                                       weighted_networks_path = NULL) {
  message("Preloading NicheNet databases...")

  start_total <- Sys.time()

  # Use parallel loading on Unix-like systems
  if (requireNamespace("parallel", quietly = TRUE) && Sys.info()["sysname"] != "Windows") {
    n_cores <- min(3, parallel::detectCores())
    
    paths <- list(lr_network_path, ligand_target_matrix_path)
    if (!is.null(weighted_networks_path)) paths <- c(paths, list(weighted_networks_path))

    results <- parallel::mclapply(paths, function(path) {
      load_ligand_target_matrix(path, use_qs = TRUE, cache = TRUE)
    }, mc.cores = n_cores)

    lr_network <- results[[1]]
    ligand_target_matrix <- results[[2]]
    weighted_networks <- if (length(results) > 2) results[[3]] else NULL
  } else {
    # Sequential loading on Windows
    lr_network <- load_ligand_target_matrix(lr_network_path, use_qs = TRUE, cache = TRUE)
    ligand_target_matrix <- load_ligand_target_matrix(ligand_target_matrix_path, use_qs = TRUE, cache = TRUE)
    weighted_networks <- if (!is.null(weighted_networks_path)) {
      load_ligand_target_matrix(weighted_networks_path, use_qs = TRUE, cache = TRUE)
    } else {
      NULL
    }
  }

  total_time <- as.numeric(difftime(Sys.time(), start_total, units = "secs"))
  message(sprintf("Preloading complete in %.2f seconds", total_time))

  return(list(
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    weighted_networks = weighted_networks,
    load_time = total_time
  ))
}

#' Clear Database Cache
#' 
#' @description Removes all cached databases from memory.
#' 
#' @export
clear_database_cache <- function() {
  cache_items <- ls(.nichenet_db_cache)
  if (length(cache_items) > 0) {
    total_size <- sum(sapply(cache_items, function(x) {
      object.size(.nichenet_db_cache[[x]])
    })) / 1024^2
    rm(list = cache_items, envir = .nichenet_db_cache)
    message(sprintf(
      "Cleared %d cached items (%.1f MB freed)",
      length(cache_items), total_size
    ))
  } else {
    message("Cache is empty")
  }
  invisible(TRUE)
}

#' Get Cache Statistics
#' 
#' @description Shows information about cached databases.
#' 
#' @export
get_database_cache_stats <- function() {
  cache_items <- ls(.nichenet_db_cache)

  if (length(cache_items) == 0) {
    message("No databases cached")
    return(invisible(NULL))
  }

  stats <- lapply(cache_items, function(key) {
    obj <- .nichenet_db_cache[[key]]
    size_mb <- object.size(obj) / 1024^2
    dims <- if (is.matrix(obj)) dim(obj) else c(NA, NA)

    data.frame(
      cache_id = key,
      type = class(obj)[1],
      rows = dims[1],
      cols = dims[2],
      size_mb = round(size_mb, 1),
      stringsAsFactors = FALSE
    )
  })

  stats_df <- do.call(rbind, stats)
  print(stats_df)
  cat(sprintf("\nTotal cached: %.1f MB\n", sum(stats_df$size_mb)))
  invisible(stats_df)
}
