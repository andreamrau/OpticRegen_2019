#----------------------------------------------------------------------------
# Function slightly modifying sleuth_prep so that I can change target_ids
#----------------------------------------------------------------------------

my_sleuth_prep <- function(
  sample_to_covariates,
  full_model = NULL,
  filter_fun = sleuth:::basic_filter,
  target_mapping = NULL,
  max_bootstrap = NULL,
  norm_fun_counts = sleuth:::norm_factors,
  norm_fun_tpm = sleuth:::norm_factors,
  aggregation_column = NULL,
  read_bootstrap_tpm = FALSE,
  extra_bootstrap_summary = FALSE,
  transformation_function = sleuth:::log_transform,
  num_cores = max(1L, parallel::detectCores() - 1L),
  ## AR add-on:
  targetIds_to_transcripts = NULL,
  ...) {
  
  ##############################
  # check inputs
  
  # data types
  
  if (!is(sample_to_covariates, "data.frame")) {
    stop(paste0("'", substitute(sample_to_covariates), "' (sample_to_covariates) must be a data.frame"))
  }
  
  if (!is(full_model, "formula") && !is(full_model, "matrix") && !is.null(full_model)) {
    stop(paste0("'", substitute(full_model), "' (full_model) must be a formula or a matrix"))
  }
  
  if (!("sample" %in% colnames(sample_to_covariates))) {
    stop(paste0("'", substitute(sample_to_covariates),
                "' (sample_to_covariates) must contain a column named 'sample'"))
  }
  
  if (!("path" %in% colnames(sample_to_covariates))) {
    stop(paste0("'", substitute(sample_to_covariates)),
         "' (sample_to_covariates) must contain a column named 'path'")
  }
  
  if (!is.null(target_mapping) && !is(target_mapping, 'data.frame')) {
    stop(paste0("'", substitute(target_mapping),
                "' (target_mapping) must be a data.frame or NULL"))
  } else if (is(target_mapping, 'data.frame')){
    if (!("target_id" %in% colnames(target_mapping))) {
      stop(paste0("'", substitute(target_mapping),
                  "' (target_mapping) must contain a column named 'target_id'"))
    }
  }
  
  if (!is.null(max_bootstrap) && max_bootstrap <= 0 ) {
    stop("max_bootstrap must be > 0")
  }
  
  if (any(is.na(sample_to_covariates))) {
    warning("Your 'sample_to_covariance' data.frame contains NA values. This will likely cause issues later.")
  }
  
  if (is(full_model, "matrix") &&
      nrow(full_model) != nrow(sample_to_covariates)) {
    stop("The design matrix number of rows are not equal to the number of rows in the sample_to_covariates argument.")
  }
  
  if (!is(norm_fun_counts, 'function')) {
    stop("norm_fun_counts must be a function")
  }
  
  if (!is(norm_fun_tpm, 'function')) {
    stop("norm_fun_tpm must be a function")
  }
  
  if (!is.null(aggregation_column) && is.null(target_mapping)) {
    stop(paste("You provided a 'aggregation_column' to aggregate by,",
               "but not a 'target_mapping'. Please provided a 'target_mapping'."))
  }
  
  num_cores <- sleuth:::check_num_cores(num_cores)
  
  # TODO: ensure transcripts are in same order -- if not, report warning that
  # kallisto index might be incorrect
  
  # done
  ##############################
  
  sleuth:::msg('reading in kallisto results')
  sample_to_covariates$sample <- as.character(sample_to_covariates$sample)
  
  if(nrow(sample_to_covariates) == 1 && !is.null(full_model)) {
    warning("There is only one sample present, but you also provided a model. ",
            "The model will be set to NULL to prevent downstream errors.\n",
            "The sample can be viewed using sleuth_live after preparation, ",
            "but you need more than one sample to run the other aspects of Sleuth.")
    full_model <- NULL
  }
  
  kal_dirs <- sample_to_covariates$path
  sample_to_covariates$path <- NULL
  
  sleuth:::msg('dropping unused factor levels')
  sample_to_covariates <- droplevels(sample_to_covariates)
  
  nsamp <- 0
  # append sample column to data
  kal_list <- lapply(seq_along(kal_dirs),
                     function(i) {
                       nsamp <- sleuth:::dot(nsamp)
                       path <- kal_dirs[i]
                       suppressMessages({
                         kal <- sleuth:::read_kallisto(path, read_bootstrap = FALSE,
                                              max_bootstrap = max_bootstrap)
                       })
                       kal$abundance <- dplyr::mutate(kal$abundance,
                                                      sample = sample_to_covariates$sample[i])
                       
                       ## AR add-on: match up kal$abundance$target_id with targetIds_to_transcripts
                       if(!is.null(targetIds_to_transcripts)) {
                         targetIds_to_transcripts$target_id <- 
                           as(targetIds_to_transcripts$target_id, 
                              class(kal$abundance$target_id))
                         nr <- nrow(kal$abundance)
                         kal$abundance <- dplyr::left_join(kal$abundance, 
                                                           targetIds_to_transcripts,
                                                 by="target_id")
                         kal$abundance <- dplyr::select(kal$abundance,
                                                 -target_id, -chr, -pos)
                         kal$abundance <- dplyr::select(kal$abundance,
                                                        target_id = transcript,
                                                        est_counts,
                                                        eff_len, len, tpm, sample)
                         if(nr != nrow(kal$abundance)) warning("STOP!")
                       }
                       kal
                     })
  sleuth:::msg('')
  
  check_result <- sleuth:::check_kal_pack(kal_list)
  kal_versions <- check_result$versions
  
  obs_raw <- dplyr::bind_rows(lapply(kal_list, function(k) k$abundance))
  
  design_matrix <- NULL
  if (is(full_model, 'formula')) {
    design_matrix <- model.matrix(full_model, sample_to_covariates)
  } else if (is(full_model, 'matrix')) {
    if (is.null(colnames(full_model))) {
      stop("If matrix is supplied, column names must also be supplied.")
    }
    design_matrix <- full_model
  }
  
  if (!is.null(full_model)) {
    rownames(design_matrix) <- sample_to_covariates$sample
  }
  
  obs_raw <- dplyr::arrange(obs_raw, target_id, sample)
  
  # ###
  # # try to deal with weird ensemble names
  # ###
  # if (!is.null(target_mapping)) {
  #   tmp_names <- data.frame(target_id = kal_list[[1]]$abundance$target_id,
  #                           stringsAsFactors = FALSE)
  #   target_mapping <- sleuth:::check_target_mapping(tmp_names, target_mapping)
  #   rm(tmp_names)
  # }
  
  ret <- list(
    kal = kal_list,
    kal_versions = kal_versions,
    obs_raw = obs_raw,
    sample_to_covariates = sample_to_covariates,
    bootstrap_summary = NA,
    full_formula = full_model,
    design_matrix = design_matrix,
    target_mapping = target_mapping,
    gene_mode = !is.null(aggregation_column),
    gene_column = aggregation_column,
    transform_fun = transformation_function
  )
  
  # TODO: eventually factor this out
  normalize <- TRUE
  if (normalize ) {
    
    sleuth:::msg("normalizing est_counts")
    est_counts_spread <- sleuth:::spread_abundance_by(obs_raw, "est_counts",
                                             sample_to_covariates$sample)
    filter_bool <- apply(est_counts_spread, 1, filter_fun)
    filter_true <- filter_bool[filter_bool]
    
    sleuth:::msg(paste0(sum(filter_bool), ' targets passed the filter'))
    est_counts_sf <- norm_fun_counts(est_counts_spread[filter_bool, , drop = FALSE])
    
    filter_df <- sleuth:::adf(target_id = names(filter_true))
    
    est_counts_norm <- sleuth:::as_df(t(t(est_counts_spread) / est_counts_sf))
    
    est_counts_norm$target_id <- rownames(est_counts_norm)
    est_counts_norm <- tidyr::gather(est_counts_norm, sample, est_counts, -target_id)
    
    obs_norm <- est_counts_norm
    obs_norm$target_id <- as.character(obs_norm$target_id)
    obs_norm$sample <- as.character(obs_norm$sample)
    rm(est_counts_norm)
    
    # deal w/ TPM
    sleuth:::msg("normalizing tpm")
    tpm_spread <- sleuth:::spread_abundance_by(obs_raw, "tpm",
                                      sample_to_covariates$sample)
    tpm_sf <- norm_fun_tpm(tpm_spread[filter_bool, , drop = FALSE])
    tpm_norm <- sleuth:::as_df(t(t(tpm_spread) / tpm_sf))
    tpm_norm$target_id <- rownames(tpm_norm)
    tpm_norm <- tidyr::gather(tpm_norm, sample, tpm, -target_id)
    tpm_norm$sample <- as.character(tpm_norm$sample)
    
    sleuth:::msg('merging in metadata')
    # put everyone in the same order to avoid a slow join
    obs_norm <- dplyr::arrange(obs_norm, target_id, sample)
    tpm_norm <- dplyr::arrange(tpm_norm, target_id, sample)
    
    stopifnot(all.equal(obs_raw$target_id, obs_norm$target_id) &&
                all.equal(obs_raw$sample, obs_norm$sample))
    
    suppressWarnings({
      if (!all.equal(dplyr::select(obs_norm, target_id, sample),
                     dplyr::select(tpm_norm, target_id, sample))) {
        stop('Invalid column rows. In principle, can simply join. Please report error.')
      }
      
      # obs_norm <- dplyr::left_join(obs_norm, data.table::as.data.table(tpm_norm),
      #   by = c('target_id', 'sample'))
      obs_norm <- dplyr::bind_cols(obs_norm, dplyr::select(tpm_norm, tpm))
    })
    
    # add in eff_len and len
    obs_norm <- dplyr::bind_cols(obs_norm, dplyr::select(obs_raw, eff_len, len))
    
    
    obs_norm <- sleuth:::as_df(obs_norm)
    ret$obs_norm <- obs_norm
    ret$est_counts_sf <- est_counts_sf
    ret$filter_bool <- filter_bool
    ret$filter_df <- filter_df
    ret$obs_norm_filt <- dplyr::semi_join(obs_norm, filter_df, by = 'target_id')
    ret$tpm_sf <- tpm_sf
    
    #### This code through the for loop is a candidate for moving to another function
    path <- kal_dirs[1]
    kal_path <- sleuth:::get_kallisto_path(path)
    target_id <- as.character(rhdf5::h5read(kal_path$path, "aux/ids"))
    ## AR add-on: change to transcript IDs
    if(!is.null(targetIds_to_transcripts)) {
      targetIds_to_transcripts$target_id <- 
        as(targetIds_to_transcripts$target_id, 
           class(target_id))
      nr <- length(target_id)
      target_id_df <- data.frame(target_id = target_id, stringsAsFactors = FALSE)
      target_id_df <- dplyr::left_join(target_id_df,
                                       targetIds_to_transcripts,
                                       by="target_id")
      target_id <- dplyr::select(target_id_df, target_id=transcript)$target_id
      if(nr != length(target_id)) warning("STOP!")
    }
    num_transcripts <- length(target_id)
    ret$bs_quants <- list()
    
    which_target_id <- ret$filter_df$target_id
    
    if (ret$gene_mode) {
      sleuth:::msg(paste0("aggregating by column: ", aggregation_column))
      # Get list of IDs to aggregate on (usually genes)
      # Also get the filtered list and update the "filter_df" and "filter_bool"
      # variables for the sleuth object
      target_mapping <- data.table::data.table(target_mapping)
      target_mapping[target_mapping[[aggregation_column]] == "",
                     aggregation_column] <- NA
      agg_id <- unique(target_mapping[, aggregation_column, with = FALSE])
      agg_id <- agg_id[[1]]
      agg_id <- agg_id[!is.na(agg_id)]
      mappings <- dplyr::select_(target_mapping, "target_id", aggregation_column)
      mappings <- data.table::as.data.table(mappings)
      which_tms <- which(mappings$target_id %in% which_target_id)
      which_agg_id <- unique(mappings[which_tms, aggregation_column, with = FALSE])
      which_agg_id <- which_agg_id[[1]]
      which_agg_id <- which_agg_id[!is.na(which_agg_id)]
      filter_df <- sleuth:::adf(target_id = which_agg_id)
      filter_bool <- agg_id %in% which_agg_id
      
      sleuth:::msg(paste0(length(which_agg_id), " genes passed the filter"))
      
      # Taken from gene_summary; scale normalized observed counts to "reads/base"
      norm_by_length <- TRUE
      tmp <- data.table::as.data.table(ret$obs_raw)
      tmp <- merge(tmp, mappings,
                   by = "target_id", all.x = TRUE)
      scale_factor <- tmp[, scale_factor := median(eff_len),
                          by=list(sample,eval(parse(text=aggregation_column)))]
      obs_norm_gene <- sleuth:::reads_per_base_transform(ret$obs_norm,
                                                scale_factor, aggregation_column, mappings, norm_by_length)
      # New code: get gene-level TPM (simple sum of normalized transcript TPM)
      tmp <- data.table::as.data.table(tpm_norm)
      tmp <- merge(tmp, mappings,
                   by = "target_id", all.x = T)
      if (any(is.na(tmp[[aggregation_column]]))) {
        rows_to_remove <- is.na(tmp[[aggregation_column]])
        num_missing <- length(unique(tmp[rows_to_remove, target_id]))
        warning(num_missing, " target_ids are missing annotations for the aggregation_column: ",
                aggregation_column, ".\nThese target_ids will be dropped from the gene-level analysis.",
                "\nIf you did not expect this, check your 'target_mapping' table for missing values.")
        tmp <- tmp[!rows_to_remove]
      }
      tpm_norm_gene <- tmp[, j = list(tpm = sum(tpm)),
                           by = list(sample, eval(parse(text = aggregation_column)))]
      data.table::setnames(tpm_norm_gene, 'parse', 'target_id')
      tpm_norm_gene <- sleuth:::as_df(tpm_norm_gene)
      
      # Same steps as above to add TPM column to "obs_norm" table
      obs_norm_gene <- dplyr::arrange(obs_norm_gene, target_id, sample)
      tpm_norm_gene <- dplyr::arrange(tpm_norm_gene, target_id, sample)
      
      stopifnot(all.equal(dplyr::select(obs_norm_gene, target_id, sample),
                          dplyr::select(tpm_norm_gene, target_id, sample)))
      suppressWarnings({
        if ( !all.equal(dplyr::select(obs_norm_gene, target_id, sample),
                        dplyr::select(tpm_norm_gene, target_id, sample), check.attributes = FALSE) ) {
          stop('Invalid column rows. In principle, can simply join. Please report error.')
        }
        
        # obs_norm <- dplyr::left_join(obs_norm, data.table::as.data.table(tpm_norm),
        #   by = c('target_id', 'sample'))
        obs_norm_gene <- dplyr::bind_cols(obs_norm_gene, dplyr::select(tpm_norm_gene, tpm))
      })
      
      # These are the updated gene-level variables
      ret$filter_df <- sleuth:::adf(target_id = which_agg_id)
      ret$filter_bool <- agg_id %in% which_agg_id
      ret$obs_norm <- obs_norm_gene
      ret$obs_norm_filt <- dplyr::semi_join(obs_norm_gene, filter_df, by = 'target_id')
      
      rm(obs_norm, tpm_norm, obs_norm_gene, tpm_norm_gene)
      
      # This is the gene-level version of the matrix
      all_sample_bootstrap <- matrix(NA_real_,
                                     nrow = length(which_agg_id),
                                     ncol = length(ret$kal))
      which_ids <- which_agg_id
    } else {
      all_sample_bootstrap <- matrix(NA_real_,
                                     nrow = length(which_target_id),
                                     ncol = length(ret$kal))
      which_ids <- which_target_id
    }
    
    sleuth:::msg('summarizing bootstraps')
    apply_function <- if (num_cores == 1) {
      lapply
    } else {
      function(x, y) parallel::mclapply(x, y, mc.cores = num_cores)
    }
    bs_results <- apply_function(seq_along(kal_dirs), function(i) {
      samp_name <- sample_to_covariates$sample[i]
      kal_path <- sleuth:::get_kallisto_path(kal_dirs[i])
      sleuth:::process_bootstrap(i, samp_name, kal_path,
                        num_transcripts, est_counts_sf[[i]],
                        read_bootstrap_tpm, ret$gene_mode,
                        extra_bootstrap_summary,
                        target_id, mappings, which_ids, ret$gene_column,
                        ret$transform_fun)
    })
    
    # if mclapply results in an error (a warning is shown), then print error and stop
    error_status <- sapply(bs_results, function(x) is(x, "try-error"))
    if (any(error_status)) {
      print(attributes(bs_results[error_status])$condition)
      stop("At least one core from mclapply had an error. See the above error message(s) for more details.")
    }
    
    # mclapply is expected to retun the bootstraps in order; this is a sanity check of that
    indices <- sapply(bs_results, function(result) result$index)
    stopifnot(identical(indices, order(indices)))
    
    ## AR add-on: match up bs_results[[*]]$bs_quants[[*]] with targetIds_to_transcripts
    if(!is.null(targetIds_to_transcripts)) {
      for(ii in seq_len(length(bs_results)))
        for(jj in seq_len(length(bs_results[[ii]]$bs_quants))) {
          tmp <- bs_results[[ii]]$bs_quants[[jj]]
          targetIds_to_transcripts$target_id <- 
            as(targetIds_to_transcripts$target_id, 
               class(rownames(tmp)))
          nr <- nrow(tmp)
          tmp_df <- data.frame(target_id = rownames(tmp),
                               tmp, stringsAsFactors=FALSE)
          tmp_df <- dplyr::left_join(tmp_df, targetIds_to_transcripts,
                                            by="target_id")
          rownames(tmp_df) <- tmp_df$transcript
          tmp_df <- dplyr::select(tmp_df, -target_id, -chr, -pos, -transcript)
          if(nr != nrow(tmp_df)) warning("STOP!")
          tmp_mat <- as.matrix(tmp_df)
          bs_results[[ii]]$bs_quants[[jj]] <- tmp_mat
        }
    }
    
    if(read_bootstrap_tpm | extra_bootstrap_summary) {
      ret$bs_quants <- lapply(bs_results, function(result) result$bs_quants)
      names(ret$bs_quants) <- sample_to_covariates$sample
    }
    
    all_sample_bootstrap <- sapply(bs_results, function(result) result$bootstrap_result)
    rownames(all_sample_bootstrap) <- which_ids
    
    # end summarize bootstraps
    sleuth:::msg('')
    
    sigma_q_sq <- rowMeans(all_sample_bootstrap)
    
    # This is the rest of the gene_summary code
    if (ret$gene_mode) {
      names(sigma_q_sq) <- which_agg_id
      obs_counts <- sleuth:::obs_to_matrix(ret, "scaled_reads_per_base")[which_agg_id, , drop = FALSE]
    } else {
      names(sigma_q_sq) <- which_target_id
      obs_counts <- sleuth:::obs_to_matrix(ret, "est_counts")[which_target_id, , drop = FALSE]
    }
    
    sigma_q_sq <- sigma_q_sq[order(names(sigma_q_sq))]
    obs_counts <- ret$transform_fun(obs_counts)
    obs_counts <- obs_counts[order(rownames(obs_counts)),]
    
    ret$bs_summary <- list(obs_counts = obs_counts, sigma_q_sq = sigma_q_sq)
  }
  
  class(ret) <- 'sleuth'
  
  ret
}


#----------------------------------------------------------------------------
# Function slightly modifying my_plot_transcript_heatmap so that the sample
# and row names can both be clustered
#----------------------------------------------------------------------------

my_plot_transcript_heatmap <- function (obj, transcripts, units = "tpm", 
                                        trans = "log", offset = 1,
                                        clustRows = TRUE,
                                        clustCols = TRUE) 
{
  units <- sleuth:::check_quant_mode(obj, units)
  if (!all(transcripts %in% obj$obs_norm$target_id)) {
    stop("Couldn't find the following transcripts: ", paste(transcripts[!(transcripts %in% 
                                                                            obj$obs_norm$target_id)], collapse = ", "), "\n\tIt is highly likely that some of them were filtered out.")
  }
  tabd_df <- obj$obs_norm[obj$obs_norm$target_id %in% transcripts, 
                          ]
  if (units == "tpm") {
    tabd_df <- dplyr::select(tabd_df, target_id, sample, 
                             tpm)
    tabd_df <- reshape2::dcast(tabd_df, target_id ~ sample, 
                               value.var = "tpm")
  } else if (units == "est_counts") {
    tabd_df <- dplyr::select(tabd_df, target_id, sample, 
                             est_counts)
    tabd_df <- reshape2::dcast(tabd_df, target_id ~ sample, 
                               value.var = "est_counts")
  } else {
    stop("Didn't recognize the following unit: ", units)
  }
  rownames(tabd_df) <- tabd_df$target_id
  tabd_df$target_id <- NULL
  p <- NULL
  if (nchar(trans) > 0 && !is.null(trans)) {
    tFunc <- eval(parse(text = trans))
    p <- sleuth:::ggPlotExpression(as.matrix(tFunc(tabd_df + offset)), 
                                   clustRows = clustRows, 
                                   clustCols = clustCols)
  } else {
    p <- sleuth:::ggPlotExpression(as.matrix(tabd_df), 
                                   clustRows = clustRows,
                                   clustCols = clustCols)
  }
  p
}


#----------------------------------------------------------------------------
# Misc
#----------------------------------------------------------------------------
my_plot_bootstrap <- function (obj, target_id, units = "est_counts",
                               color_by = setdiff(colnames(obj$sample_to_covariates),
                                                                   "sample"),
                               x_axis_angle = 50, divide_groups = TRUE)
{
  units <- sleuth:::check_quant_mode(obj, units)
  df <- sleuth:::get_bootstrap_summary(obj, target_id, units)
  p <- ggplot(df, aes(x = sample, ymin = min, lower = lower,
                      middle = mid, upper = upper, ymax = max))
  p <- p + geom_boxplot(stat = "identity", aes_string(fill = color_by))
  p <- p + theme(axis.text.x = element_text(angle = x_axis_angle,
                                            hjust = 1))
  p <- p + ggtitle(target_id)
  p <- p + ylab(units)
  if (divide_groups) {
    p <- p + facet_wrap(color_by, switch = "x", scales = "free_x")
  }
  p
}
my_get_bootstrap_summary <- function (obj, target_id, units = "est_counts")
{
  stopifnot(is(obj, "sleuth"))
  if (units != "est_counts" && units != "tpm" && units != "scaled_reads_per_base") {
    stop(paste0("'", units, "' is invalid for 'units'. please see documentation"))
  }
  if (is.null(obj$bs_quants)) {
    if (units == "est_counts") {
      stop("bootstrap summary missing. rerun sleuth_prep() with argument 'extra_bootstrap_summary = TRUE'")
    } else {
      stop("bootstrap summary missing. rerun sleuth_prep() with argument 'extra_bootstrap_summary = TRUE' and 'read_bootstrap_tpm = TRUE'")
    }
  }
  if (!(target_id %in% rownames(obj$bs_quants[[1]][[units]]))) {
    stop(paste0("couldn't find target_id '", target_id, "'"))
  }
  df <- sleuth:::as_df(do.call(rbind, lapply(obj$bs_quants, function(sample_bs) {
    sample_bs[[units]][target_id, ]
  })))
  df <- dplyr::bind_cols(df, obj$sample_to_covariates)
  df
}


