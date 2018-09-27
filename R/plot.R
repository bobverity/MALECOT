
#------------------------------------------------
# default MALECOT colours
#' @noRd
default_colours <- function(K) {
  
  # generate palette and colours
  raw_cols <- c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4")
  #raw_cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
  my_palette <- colorRampPalette(raw_cols)
  
  # simple case if small K
  if (K<=2) {
    return(my_palette(K))
  }
  
  # some logic to choose a palette size and sequence of colours that is
  # consistent across different values of K
  ncol <- 3
  while(ncol<K) {
    ncol <- ncol+(ncol-1)
  }
  dist_mat <- matrix(1:ncol, ncol, ncol)
  dist_mat <- abs(t(dist_mat)-dist_mat)
  x <- rep(FALSE, ncol)
  
  col_index <- 1
  for (i in 2:K) {
    x[col_index] <- TRUE
    s <- apply(dist_mat[which(x),,drop=FALSE], 2, min)
    next_index <- which.max(s)
    col_index <- c(col_index, next_index)
  }
  col_index
  ret <- my_palette(ncol)[col_index]
  
  return(ret)
}

#------------------------------------------------
#' @title Produce PCA plot from MALECOT data
#'
#' @description TODO
#'
#' @details TODO
#'
#' @param data TODO
#' @param type TODO
#' @param missing_data TODO
#' @param target_group TODO
#'
#' @export
#' @examples
#' # TODO

plot_pca <- function(data, type = "3D", missing_data = -1, target_group = NULL) {
  
  # TODO - check inputs
  assert_in(type, c("2D", "3D"))
  
  # TODO - get data format automatically
  
  # define defaults
  n <- nrow(data)
  target_group <- define_default(target_group, rep(1, n))
  
  # impute missing with mean allele frequencies per locus
  locus_means <- colMeans(data, na.rm=TRUE)
  data <- apply(data, 2, function(x){
    x[x==missing_data] <- locus_means[x==missing_data]
    return(x)})
  
  # compute PCA
  PCA <- prcomp(data)
  
  # compute variance explained
  PCA$var <- PCA$sdev^2/sum(PCA$sdev^2) * 100
  
  # barchart
  plot_bar <- plot_ly(x = colnames(PCA$x)[1:9], y=PCA$var[1:9], type="bar", width=500, height=500)
  
  # 2D or 3D scatter
  col_vec <- default_colours(length(unique(target_group)))
  
  if (type=="3D") {
    plot_scat <- plot_ly(as.data.frame(PCA$x[,1:3]), x = ~PC1, y = ~PC2, z = ~PC3, color = target_group, colors = col_vec, type = "scatter3d", mode = "markers", marker = list(size=5))
  } else {
    plot_scat <- plot_ly(as.data.frame(PCA$x), x = ~PC1, y = ~PC2, color = target_group, colors = col_vec, type = "scatter", mode = "markers", marker = list(size=10))
  }
  print(plot_scat)
  
  # return list of plots invisibly
  plot_list <- list(plot_bar, plot_scat)
  invisible(plot_list)
}

#------------------------------------------------
# ggplot theme with minimal objects
#' @noRd
theme_empty <- function() {
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
}

#------------------------------------------------
# Default plot for class malecot_qmatrix
#' @noRd
plot.malecot_qmatrix <- function(x, y, ...) {
  
  # get data into ggplot format
  m <- unclass(x)
  n <- nrow(m)
  K <- ncol(m)
  df <- data.frame(ind = rep(1:n,each=K), k = as.factor(rep(1:K,times=n)), val = as.vector(t(m)))
  
  # produce basic plot
  plot1 <- ggplot(df) + theme_empty()
  plot1 <- plot1 + geom_bar(aes_(x = ~ind, y = ~val, fill = ~k), width = 1, stat = "identity")
  plot1 <- plot1 + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
  plot1 <- plot1 + xlab("sample") + ylab("probability")
  
  # add legends
  plot1 <- plot1 + scale_fill_manual(values = default_colours(K), name = "group")
  plot1 <- plot1 + scale_colour_manual(values = "white")
  plot1 <- plot1 + guides(colour = FALSE)
  
  # add border
  plot1 <- plot1 + theme(panel.border = element_rect(colour = "black", size = 2, fill = NA))
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Posterior allocation plot
#'
#' @description Produce posterior allocation plot of current active set.
#'
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#' @param K which value of K to produce the plot for
#' @param divide_ind_on whether to add dividing lines between bars
#'
#' @export

plot_structure <- function(project, K = NULL, divide_ind_on = FALSE) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  if (!is.null(K)) {
    assert_pos_int(K)
  }
  assert_single_logical(divide_ind_on)
  
  # get active set and check non-zero
  s <- project$active_set
  if (s==0) {
    stop("no active parameter set")
  }
  
  # set default K to all values with output
  null_output <- mapply(function(x) {is.null(x$summary$qmatrix)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no output for active parameter set")
  }
  K <- define_default(K, which(!null_output))
  
  # check output exists for chosen K
  qmatrix_list <- list()
  for (i in 1:length(K)) {
    qmatrix_list[[i]] <- project$output$single_set[[s]]$single_K[[K[i]]]$summary$qmatrix
    if (is.null(qmatrix_list[[i]])) {
      stop(sprintf("no qmatrix output for K = %s of active set", K[i]))
    }
  }
  n <- nrow(qmatrix_list[[1]])
  
  # get data into ggplot format
  df <- NULL
  for (i in 1:length(K)) {
    m <- unclass(qmatrix_list[[i]])
    df <- rbind(df, data.frame(K = as.numeric(K[i]), ind = rep(1:n,each=K[i]), k = as.factor(rep(1:K[i],times=n)), val = as.vector(t(m))))
  }
  
  # produce basic plot
  plot1 <- ggplot(df) + theme_empty()
  plot1 <- plot1 + geom_bar(aes_(x = ~ind, y = ~val, fill = ~k), width = 1, stat = "identity")
  plot1 <- plot1 + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
  
  # arrange in rows
  if (length(K)==1) {
    plot1 <- plot1 + facet_wrap(~K, ncol = 1)
    plot1 <- plot1 + theme(strip.background = element_blank(), strip.text = element_blank())
    plot1 <- plot1 + xlab("sample") + ylab("probability")
  } else {
    plot1 <- plot1 + facet_wrap(~K, ncol = 1, strip.position = "left")
    plot1 <- plot1 + theme(strip.background = element_blank())
    plot1 <- plot1 + xlab("sample") + ylab("K")
  }
  
  # add legends
  plot1 <- plot1 + scale_fill_manual(values = default_colours(max(K)), name = "group")
  plot1 <- plot1 + scale_colour_manual(values = "white")
  plot1 <- plot1 + guides(colour = FALSE)
  
  # add border
  plot1 <- plot1 + theme(panel.border = element_rect(colour = "black", size = 2, fill = NA))
  
  # optionally add dividing lines
  if (divide_ind_on) {
    plot1 <- plot1 + geom_segment(aes_(x = ~x, y = ~y, xend = ~x, yend = ~y+1, col = "white"), size = 0.3, data = data.frame(x = 1:n-0.5, y = rep(0,n)))
  }
  
  return(plot1)
}

#------------------------------------------------
# Default plot for class malecot_loglike_intervals
#' @noRd
plot.malecot_loglike_intervals <- function(x, y, ...) {
  
  # get data into ggplot format
  df <- as.data.frame(unclass(x))
  n <- nrow(df)
  x_vec <- y
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw()
  plot1 <- plot1 + geom_segment(aes_(x = ~x_vec, y = ~Q2.5, xend = ~x_vec, yend = ~Q97.5))
  plot1 <- plot1 + geom_point(aes_(x = ~x_vec, y = ~Q50))
  plot1 <- plot1 + xlab("rung") + ylab("log-likelihood")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot loglikelihood 95\% credible intervals
#'   
#' @description Plot loglikelihood 95\% credible intervals of current active set
#'   
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#' @param K which value of K to produce the plot for
#' @param axis_type how to format the x-axis. 1 = integer rungs, 2 = values of
#'   beta, 3 = values of beta raised to the GTI power
#' @param connect_points whether to connect points in the middle of intervals
#' @param connect_whiskers whether to connect points at the ends of the whiskers
#'
#' @export

plot_loglike <- function(project, K = NULL, axis_type = 1, connect_points = FALSE, connect_whiskers = FALSE) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  if (!is.null(K)) {
    assert_single_pos_int(K, zero_allowed = FALSE)
  }
  assert_in(axis_type, 1:3)
  assert_single_logical(connect_points)
  assert_single_logical(connect_whiskers)
  
  # get active set and check non-zero
  s <- project$active_set
  if (s==0) {
    stop("no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$summary$loglike_intervals)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no loglike_intervals output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  loglike_intervals <- project$output$single_set[[s]]$single_K[[K]]$summary$loglike_intervals
  if (is.null(loglike_intervals)) {
    stop(sprintf("no loglike_intervals output for K = %s of active set", K))
  }
  
  # produce plot with different axis options
  rungs <- nrow(loglike_intervals)
  if (axis_type==1) {
    x_vec <- 1:rungs
    plot1 <- plot(loglike_intervals, as.factor(x_vec))
    
  } else if (axis_type==2) {
    x_vec <- (1:rungs)/rungs
    plot1 <- plot(loglike_intervals, x_vec)
    plot1 <- plot1 + xlab(parse(text = "beta"))
    plot1 <- plot1 + coord_cartesian(xlim = c(0,1))
    
  } else {
    GTI_pow <- project$output$single_set[[s]]$single_K[[K]]$function_call$args$GTI_pow
    x_vec <- ((1:rungs)/rungs)^GTI_pow
    plot1 <- plot(loglike_intervals, x_vec)
    plot1 <- plot1 + xlab(parse(text = "beta^gamma"))
    plot1 <- plot1 + coord_cartesian(xlim = c(0,1))
  }
  
  # optionally add central line
  if (connect_points) {
    df <- as.data.frame(unclass(loglike_intervals))
    plot1 <- plot1 + geom_line(aes(x = x_vec, y = df$Q50))
  }
  
  # optionally connect whiskers
  if (connect_whiskers) {
    df <- as.data.frame(unclass(loglike_intervals))
    plot1 <- plot1 + geom_line(aes(x = x_vec, y = df$Q2.5), linetype = "dotted") + geom_line(aes(x = x_vec, y = df$Q97.5), linetype = "dotted")
  }
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
# Default plot for class malecot_COI_intervals
#' @noRd
plot.malecot_COI_intervals <- function(x, y, ...) {
  
  # get data into ggplot format
  df <- as.data.frame(unclass(x))
  n <- nrow(df)
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw()
  plot1 <- plot1 + geom_segment(aes_(x = ~1:n, y = ~Q2.5, xend = ~1:n, yend = ~Q97.5))
  plot1 <- plot1 + geom_point(aes_(x = ~1:n, y = ~Q50))
  plot1 <- plot1 + scale_y_continuous(limits = c(0, max(df$Q97.5)*1.1), expand = c(0,0))
  plot1 <- plot1 + xlab("sample") + ylab("COI")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot COI 95\% credible intervals
#'
#' @description Plot COI 95\% credible intervals of current active set
#'
#' @details TODO
#'
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#' @param K which value of K to produce the plot for
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

plot_COI <- function(project, K = NULL, ...) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  if (!is.null(K)) {
    assert_single_pos_int(K, zero_allowed = FALSE)
  }
  
  # get active set and check non-zero
  s <- project$active_set
  if (s == 0) {
    stop("no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$summary$COI_intervals)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no COI_intervals output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  COI_intervals <- project$output$single_set[[s]]$single_K[[K]]$summary$COI_intervals
  if (is.null(COI_intervals)) {
    stop(sprintf("no COI_intervals output for K = %s of active set", K))
  }
  
  # produce quantile plot
  plot1 <- plot(project$output$single_set[[s]]$single_K[[K]]$summary$COI_intervals)
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
# Default plot for class malecot_p_intervals
#' @noRd
plot.malecot_p_intervals <- function(x, y, ...) {
  
  # get maximum number of alleles over all loci
  max_alleles <- max(mapply(nrow, x))
  
  # split plotting method for bi-allelic vs. multi-allelic estimates
  if (max_alleles == 1) {
    plot1 <- plot_p_biallelic(x)
  } else {
    plot1 <- plot_p_multiallelic(x)
  }
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
# plot bi-allelic allele frequency estimates
#' @noRd
plot_p_biallelic <- function(p) {
  
  # get data into ggplot format
  df <- as.data.frame(t(mapply(function(x){x}, p)))
  names(df) <- c("Q2.5", "Q50", "Q97.5")
  n <- nrow(df)
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw()
  plot1 <- plot1 + geom_segment(aes_(x = ~1:n, y = ~Q2.5, xend = ~1:n, yend = ~Q97.5))
  plot1 <- plot1 + geom_point(aes_(x = ~1:n, y = ~Q50))
  plot1 <- plot1 + xlab("locus") + ylab("allele frequency")
  plot1 <- plot1 + ylim(0,1)
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
# plot multi-allelic allele frequency estimates
#' @noRd
plot_p_multiallelic <- function(x) {
  
  # get maximum number of alleles at any locus
  L <- length(x)
  max_alleles <- max(mapply(nrow, x))
  
  # add empty rows so same number of alleles at every locus
  x_expanded <- mapply(function(y) {
    ret <- y
    if (nrow(ret) < max_alleles) {
      ret <- rbind(ret, matrix(0, max_alleles - nrow(ret), ncol(ret)))
    }
    cbind(allele = 1:nrow(ret), ret)
  }, x, SIMPLIFY = FALSE)
  
  # get into dataframe
  df <- do.call(rbind.data.frame, x_expanded)
  df$locus <- rep(1:L, each = max_alleles)
  
  # bar colours
  col_vec <- brewer.pal(max_alleles+3, "Blues")[-(1:3)]
  
  # create basic plot
  plot1 <- ggplot(df, aes_(x = ~as.factor(locus), y = ~Q50, fill = ~as.factor(allele)))  + theme_bw() + theme(panel.grid.major.x = element_blank())
  
  # add grouped barpot
  plot1 <- plot1 + geom_bar(stat = "identity", position = "dodge")
  
  # add error bars
  plot1 <- plot1 + geom_errorbar(aes_(ymin = ~Q2.5, ymax = ~Q97.5), position = "dodge")
  
  # scales, titles, legends etc.
  plot1 <- plot1 + scale_fill_manual(values = col_vec, name = "allele")
  plot1 <- plot1 + scale_y_continuous(limits = c(0,1), expand = c(0,0))
  plot1 <- plot1 + geom_vline(xintercept = 1:L + 0.5, col = grey(0.9))
  plot1 <- plot1 + xlab("locus") + ylab("allele frequency")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot allele frequency 95\% credible intervals
#'
#' @description Plot allele frequency 95\% credible intervals of current active set
#'
#' @details TODO
#'
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#' @param K TODO
#' @param deme TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

plot_p <- function(project, K = NULL, deme = 1, ...) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  if (!is.null(K)) {
    assert_single_pos_int(K, zero_allowed = FALSE)
  }
  assert_single_pos_int(deme)
  
  # get active set and check non-zero
  s <- project$active_set
  if (s == 0) {
    stop("no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$summary$p_intervals)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no p_intervals output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  p_intervals <- project$output$single_set[[s]]$single_K[[K]]$summary$p_intervals
  if (is.null(p_intervals)) {
    stop(sprintf("no p_intervals output for K = %s of active set", K))
  }
  
  # check that selected deme is within limits
  assert_leq(deme, length(p_intervals))
  
  # produce quantile plot
  plot1 <- plot(project$output$single_set[[s]]$single_K[[K]]$summary$p_intervals[[deme]])
  
  # add deme title
  plot1 <- plot1 + ggtitle(sprintf("deme%s", deme))
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
# Default plot for class malecot_e_intervals
#' @noRd
plot.malecot_e_intervals <- function(x, y, ...) {
  
  # check for NULL
  if (is.null(x)) {
    message("no e_intervals to plot")
  }
  
  # get data into ggplot format
  df <- as.data.frame(unclass(x))
  xvals <- c("e1", "e2")
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw()
  plot1 <- plot1 + geom_segment(aes_(x = xvals, y = ~Q2.5, xend = xvals, yend = ~Q97.5))
  plot1 <- plot1 + geom_point(aes_(x = xvals, y = ~Q50))
  plot1 <- plot1 + xlab("") + ylab("posterior error")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot error rate 95\% credible intervals
#'
#' @description Plot error rate 95\% credible intervals of current active set
#'
#' @details TODO
#'
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#' @param K TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

plot_e <- function(project, K = NULL, ...) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  if (!is.null(K)) {
    assert_single_pos_int(K, zero_allowed = FALSE)
  }
  
  # get active set and check non-zero
  s <- project$active_set
  if (s == 0) {
    stop("no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$summary$e_intervals)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no e_intervals output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  e_intervals <- project$output$single_set[[s]]$single_K[[K]]$summary$e_intervals
  if (is.null(e_intervals)) {
    stop(sprintf("no e_intervals output for K = %s of active set", K))
  }
  
  # produce quantile plot
  plot1 <- plot(project$output$single_set[[s]]$single_K[[K]]$summary$e_intervals)
  
  # set limits based on e1_max and e2_max
  e1_max <- project$parameter_sets[[s]]$e1_max
  e2_max <- project$parameter_sets[[s]]$e2_max
  plot1 <- plot1 + scale_y_continuous(limits = c(0,max(e1_max, e2_max)), expand = c(0,0))
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
# Default plot for class malecot_COI_mean_intervals
#' @noRd
plot.malecot_COI_mean_intervals <- function(x, y, ...) {
  
  # check for NULL
  if (is.null(x)) {
    message("no COI_mean_intervals to plot")
  }
  
  # get data into ggplot format
  df <- as.data.frame(unclass(x))
  xvals <- rownames(df)
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw()
  plot1 <- plot1 + geom_segment(aes_(x = factor(xvals, levels = xvals), y = ~Q2.5, xend = xvals, yend = ~Q97.5))
  plot1 <- plot1 + geom_point(aes_(x = xvals, y = ~Q50))
  plot1 <- plot1 + xlab("") + ylab("posterior mean COI")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot COI_mean 95\% credible intervals
#'
#' @description Plot COI_mean 95\% credible intervals of current active set
#'
#' @details TODO
#'
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#' @param K TODO
#' @param deme_order the order in which to plot demes. Defaults to increasing order
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

plot_COI_mean <- function(project, K = NULL, deme_order = NULL, ...) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  if (!is.null(K)) {
    assert_single_pos_int(K, zero_allowed = FALSE)
  }
  if (!is.null(deme_order)) {
    assert_pos_int(deme_order, zero_allowed = FALSE)
    assert_vector(deme_order)
  }
  
  # get active set and check non-zero
  s <- project$active_set
  if (s == 0) {
    stop("no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$summary$COI_mean_intervals)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no COI_mean_intervals output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  COI_mean_intervals <- project$output$single_set[[s]]$single_K[[K]]$summary$COI_mean_intervals
  if (is.null(COI_mean_intervals)) {
    stop(sprintf("no COI_mean_intervals output for K = %s of active set", K))
  }
  
  # get quantile object in chosen order
  df <- project$output$single_set[[s]]$single_K[[K]]$summary$COI_mean_intervals
  deme_order <- define_default(deme_order, 1:nrow(df))
  assert_eq(length(deme_order), nrow(df))
  df <- df[deme_order, , drop = FALSE]
  class(df) <- "malecot_COI_mean_intervals"
  
  # produce quantile plot
  plot1 <- plot(df)
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
# Default plot for class maverick_GTI_path
#' @noRd
plot.malecot_GTI_path <- function(x, y, ...) {
  
  # check inputs
  assert_in(y, 1:2)
  
  # get data into ggplot format
  df <- as.data.frame(unclass(x))
  n <- nrow(df)
  df <- rbind(data.frame(mean = 0, SE = 0), df)
  
  # get quantiles
  df$q_min <- df$mean - 1.96*df$SE
  df$q_mid <- df$mean
  df$q_max <- df$mean + 1.96*df$SE
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw()
  if (y==1) {
    q_x <- 1:(n+1)
    width <- 0.1
    plot1 <- plot1 + geom_line(aes_(x = ~as.factor(0:n), y = ~q_mid, group = 1))
    plot1 <- plot1 + xlab("rung")
  } else {
    q_x <- seq(0,1,l=n+1)
    width <- 0.01
    plot1 <- plot1 + geom_line(aes_(x = ~q_x, y = ~q_mid))
    plot1 <- plot1 + xlab(parse(text = "beta"))
  }
  
  # continue building plot
  plot1 <- plot1 + geom_area(aes_(x = ~q_x, y = ~q_mid, fill = "col1", colour = "col1", alpha = 0.5))
  plot1 <- plot1 + geom_segment(aes_(x = ~q_x, y = ~q_min, xend = ~q_x, yend = ~q_max))
  plot1 <- plot1 + geom_segment(aes_(x = ~q_x-width, y = ~q_min, xend = ~q_x+width, yend = ~q_min))
  plot1 <- plot1 + geom_segment(aes_(x = ~q_x-width, y = ~q_max, xend = ~q_x+width, yend = ~q_max))
  
  plot1 <- plot1 + scale_fill_manual(values = "#4575B4")
  plot1 <- plot1 + scale_colour_manual(values = "black")
  plot1 <- plot1 + guides(fill = FALSE, colour = FALSE, alpha = FALSE)
  plot1 <- plot1 + ylab("weighted log-likelihood")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot GTI path of current active set
#'
#' @description Plot GTI path of current active set
#'
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#' @param K which value of K to produce the plot for
#' @param axis_type how to format the x-axis. 1 = integer rungs, 2 = values of
#'   beta
#'
#' @export

plot_GTI_path <- function(project, K = NULL, axis_type = 1) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  if (!is.null(K)) {
    assert_single_pos_int(K, zero_allowed = FALSE)
  }
  assert_in(axis_type, 1:2)
  
  # get active set and check non-zero
  s <- project$active_set
  if (s==0) {
    stop("no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$summary$GTI_path)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  GTI_path <- project$output$single_set[[s]]$single_K[[K]]$summary$GTI_path
  if (is.null(GTI_path)) {
    stop(sprintf("no GTI_path output for K = %s of active set", K))
  }
  
  # produce plot with different axis options
  plot1 <- plot(GTI_path, axis_type)
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
# Default plot for class malecot_GTI_logevidence
#' @noRd
plot.malecot_GTI_logevidence <- function(x, y, ...) {
  
  # get data into ggplot format
  df <- as.data.frame(unclass(x))
  n <- nrow(df)
  
  # get quantiles
  df$q_min <- df$mean - 1.96*df$SE
  df$q_mid <- df$mean
  df$q_max <- df$mean + 1.96*df$SE
  q_x <- 1:n
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw()
  width <- 0.1
  plot1 <- plot1 + geom_point(aes_(x = ~as.factor(K), y = ~q_mid), na.rm = TRUE)
  plot1 <- plot1 + geom_segment(aes_(x = ~q_x, y = ~q_min, xend = ~q_x, yend = ~q_max), na.rm = TRUE)
  plot1 <- plot1 + geom_segment(aes_(x = ~q_x-width, y = ~q_min, xend = ~q_x+width, yend = ~q_min), na.rm = TRUE)
  plot1 <- plot1 + geom_segment(aes_(x = ~q_x-width, y = ~q_max, xend = ~q_x+width, yend = ~q_max), na.rm = TRUE)
  plot1 <- plot1 + xlab("K") + ylab("log-evidence")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot log-evidence estimates over K
#'
#' @description Plot log-evidence estimates over K
#'
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#'
#' @export

plot_logevidence_K <- function(project) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  
  # get active set and check non-zero
  s <- project$active_set
  if (s==0) {
    stop("no active parameter set")
  }
  
  # check output exists for chosen K
  GTI_logevidence <- project$output$single_set[[s]]$all_K$GTI_logevidence
  if (is.null(GTI_logevidence)) {
    stop("no GTI_logevidence output for active set")
  }
  
  # produce plot
  plot1 <- plot(GTI_logevidence)
  
  # return plot object
  return(plot1)
}


#------------------------------------------------
# Default plot for class malecot_GTI_posterior
#' @noRd
plot.malecot_GTI_posterior <- function(x, y, ...) {
  
  # get data into ggplot format
  df <- as.data.frame(unclass(x))
  n <- nrow(df)
  width <- 0.1
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw()
  plot1 <- plot1 + geom_bar(aes_(x = ~K, y = ~Q50, fill = "blue"), stat = "identity", na.rm = TRUE)
  plot1 <- plot1 + geom_segment(aes_(x = ~K, y = ~Q2.5, xend = ~K, yend = ~Q97.5), na.rm = TRUE)
  plot1 <- plot1 + geom_segment(aes_(x = ~K-width, y = ~Q2.5, xend = ~K+width, yend = ~Q2.5), na.rm = TRUE)
  plot1 <- plot1 + geom_segment(aes_(x = ~K-width, y = ~Q97.5, xend = ~K+width, yend = ~Q97.5), na.rm = TRUE)
  
  # add legends
  plot1 <- plot1 + scale_fill_manual(values = "#4575B4")
  plot1 <- plot1 + guides(fill = FALSE)
  
  # modify scales etc.
  plot1 <- plot1 + coord_cartesian(ylim = c(-0.05,1.05))
  plot1 <- plot1 + scale_y_continuous(expand = c(0,0))
  plot1 <- plot1 + xlab("K") + ylab("probability")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot posterior K
#'
#' @description Plot posterior K
#'
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#'
#' @export

plot_posterior_K <- function(project) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  
  # get active set and check non-zero
  s <- project$active_set
  if (s==0) {
    stop("no active parameter set")
  }
  
  # check output exists for chosen K
  GTI_posterior <- project$output$single_set[[s]]$all_K$GTI_posterior
  if (is.null(GTI_posterior)) {
    stop("no GTI_posterior output for active set")
  }
  
  # produce plot
  plot1 <- plot(GTI_posterior)
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
# Default plot for class malecot_GTI_logevidence_model
#' @noRd
plot.malecot_GTI_logevidence_model <- function(x, y, ...) {
  
  # get data into ggplot format
  df <- as.data.frame(unclass(x))
  n <- nrow(df)
  
  # get quantiles
  df$q_min <- df$mean - 1.96*df$SE
  df$q_mid <- df$mean
  df$q_max <- df$mean + 1.96*df$SE
  q_x <- 1:n
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw()
  width <- 0.1
  plot1 <- plot1 + geom_point(aes_(x = ~as.factor(set), y = ~q_mid), na.rm = TRUE)
  plot1 <- plot1 + geom_segment(aes_(x = ~q_x, y = ~q_min, xend = ~q_x, yend = ~q_max), na.rm = TRUE)
  plot1 <- plot1 + geom_segment(aes_(x = ~q_x-width, y = ~q_min, xend = ~q_x+width, yend = ~q_min), na.rm = TRUE)
  plot1 <- plot1 + geom_segment(aes_(x = ~q_x-width, y = ~q_max, xend = ~q_x+width, yend = ~q_max), na.rm = TRUE)
  plot1 <- plot1 + scale_x_discrete(labels = df$name, breaks = 1:n)
  plot1 <- plot1 + xlab("model") + ylab("log-evidence")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot log-evidence estimates over parameter sets
#'
#' @description Plot log-evidence estimates over parameter sets
#'
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#'
#' @export

plot_logevidence_model <- function(project) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  
  # check output exists
  GTI_logevidence_model <- project$output$all_sets$GTI_logevidence_model
  if (is.null(GTI_logevidence_model)) {
    stop("no GTI_logevidence_model output")
  }
  
  # produce plot
  plot1 <- plot(GTI_logevidence_model)
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
# Default plot for class malecot_GTI_posterior_model
#' @noRd
plot.malecot_GTI_posterior_model <- function(x, y, ...) {
  
  # get data into ggplot format
  df <- as.data.frame(unclass(x))
  n <- nrow(df)
  width <- 0.1
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw()
  plot1 <- plot1 + geom_bar(aes_(x = ~as.factor(set), y = ~Q50, fill = "blue"), stat = "identity", na.rm = TRUE)
  plot1 <- plot1 + geom_segment(aes_(x = ~set, y = ~Q2.5, xend = ~set, yend = ~Q97.5), na.rm = TRUE)
  plot1 <- plot1 + geom_segment(aes_(x = ~set-width, y = ~Q2.5, xend = ~set+width, yend = ~Q2.5), na.rm = TRUE)
  plot1 <- plot1 + geom_segment(aes_(x = ~set-width, y = ~Q97.5, xend = ~set+width, yend = ~Q97.5), na.rm = TRUE)
  
  # add legends
  plot1 <- plot1 + scale_fill_manual(values = "#4575B4")
  plot1 <- plot1 + guides(fill = FALSE)
  
  # modify scales etc.
  plot1 <- plot1 + scale_x_discrete(labels = df$name, breaks = 1:n)
  plot1 <- plot1 + coord_cartesian(ylim = c(-0.05,1.05))
  plot1 <- plot1 + scale_y_continuous(expand = c(0,0))
  plot1 <- plot1 + xlab("model") + ylab("probability")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot posterior model
#'
#' @description Plot posterior model
#'
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#'
#' @export

plot_posterior_model <- function(project) {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  
  # check output exists
  GTI_posterior_model <- project$output$all_sets$GTI_posterior_model
  if (is.null(GTI_posterior_model)) {
    stop("no GTI_posterior_model output")
  }
  
  # produce plot
  plot1 <- plot(GTI_posterior_model)
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Produce MCMC trace plot
#'   
#' @description Produce MCMC trace plot of the log-likelihood at each iteration.
#'   
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#' @param K which value of K to plot
#' @param rung which rung to plot. Defaults to the cold chain
#' @param col colour of the trace
#'   
#' @export

plot_trace <- function(project, K = NULL, rung = NULL, col = "black") {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  if (!is.null(K)) {
    assert_single_pos_int(K, zero_allowed = FALSE)
  }
  if (!is.null(rung)) {
    assert_single_pos_int(rung)
  }
  
  # get active set and check non-zero
  s <- project$active_set
  if (s==0) {
    stop("no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$raw$loglike_sampling)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no loglike_sampling output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  loglike_sampling <- project$output$single_set[[s]]$single_K[[K]]$raw$loglike_sampling
  if (is.null(loglike_sampling)) {
    stop(sprintf("no loglike_sampling output for K = %s of active set", K))
  }
  
  # use cold rung by default
  rungs <- ncol(loglike_sampling)
  rung <- define_default(rung, rungs)
  assert_leq(rung, rungs)
  loglike <- as.vector(loglike_sampling[,rung])
  
  # get into ggplot format
  df <- data.frame(x = 1:length(loglike), y = loglike)
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw() + ylab("log-likelihood")
  
  # complete plot
  plot1 <- plot1 + geom_line(aes_(x = ~x, y = ~y, colour = "col1"))
  plot1 <- plot1 + coord_cartesian(xlim = c(0,nrow(df)))
  plot1 <- plot1 + scale_x_continuous(expand = c(0,0))
  plot1 <- plot1 + scale_colour_manual(values = col)
  plot1 <- plot1 + guides(colour = FALSE)
  plot1 <- plot1 + xlab("iteration")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Produce MCMC autocorrelation plot
#'
#' @description Produce MCMC autocorrelation plot of the log-likelihood
#'
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#' @param K which value of K to plot
#' @param rung which rung to plot. Defaults to the cold chain
#' @param col colour of the trace
#'
#' @export

plot_acf <- function(project, K = NULL, rung = NULL, col = "black") {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  if (!is.null(K)) {
    assert_single_pos_int(K, zero_allowed = FALSE)
  }
  if (!is.null(rung)) {
    assert_single_pos_int(rung)
  }
  
  # get active set and check non-zero
  s <- project$active_set
  if (s==0) {
    stop("no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$raw$loglike_sampling)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no loglike_sampling output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  loglike_sampling <- project$output$single_set[[s]]$single_K[[K]]$raw$loglike_sampling
  if (is.null(loglike_sampling)) {
    stop(sprintf("no loglike_sampling output for K = %s of active set", K))
  }
  
  # use cold rung by default
  rungs <- ncol(loglike_sampling)
  rung <- define_default(rung, rungs)
  assert_leq(rung, rungs)
  loglike <- as.vector(loglike_sampling[,rung])
  
  # store variable to plot
  v <- loglike
  
  # get autocorrelation
  lag_max <- round(3*length(v)/effectiveSize(v))
  lag_max <- max(lag_max, 20)
  lag_max <- min(lag_max, length(v))
  
  # get into ggplot format
  a <- acf(v, lag.max = lag_max, plot = FALSE)
  acf <- as.vector(a$acf)
  df <- data.frame(lag = (1:length(acf))-1, ACF = acf)
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw()
  plot1 <- plot1 + geom_segment(aes_(x = ~lag, y = 0, xend = ~lag, yend = ~ACF, colour = "col1"))
  plot1 <- plot1 + scale_colour_manual(values = col)
  plot1 <- plot1 + guides(colour = FALSE)
  plot1 <- plot1 + xlab("lag") + ylab("ACF")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Produce MCMC density plot
#'
#' @description Produce MCMC density plot of the log-likelihood
#'
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#' @param K which value of K to plot
#' @param rung which rung to plot. Defaults to the cold chain
#' @param col colour of the trace
#'
#' @export

plot_density <- function(project, K = NULL, rung = NULL, col = "black") {
  
  # check inputs
  assert_custom_class(project, "malecot_project")
  if (!is.null(K)) {
    assert_single_pos_int(K, zero_allowed = FALSE)
  }
  if (!is.null(rung)) {
    assert_single_pos_int(rung)
  }
  
  # get active set and check non-zero
  s <- project$active_set
  if (s==0) {
    stop("no active parameter set")
  }
  
  # set default K to first value with output
  null_output <- mapply(function(x) {is.null(x$raw$loglike_sampling)}, project$output$single_set[[s]]$single_K)
  if (all(null_output)) {
    stop("no loglike_sampling output for active parameter set")
  }
  if (is.null(K)) {
    K <- which(!null_output)[1]
    message(sprintf("using K = %s by default", K))
  }
  
  # check output exists for chosen K
  loglike_sampling <- project$output$single_set[[s]]$single_K[[K]]$raw$loglike_sampling
  if (is.null(loglike_sampling)) {
    stop(sprintf("no loglike_sampling output for K = %s of active set", K))
  }
  
  # use cold rung by default
  rungs <- ncol(loglike_sampling)
  rung <- define_default(rung, rungs)
  assert_leq(rung, rungs)
  loglike <- as.vector(loglike_sampling[,rung])
  
  # get into ggplot format
  df <- data.frame(v = loglike)
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw() + xlab("log-likelihood")
  
  # produce plot
  #plot1 <- ggplot(df) + theme_bw()
  plot1 <- plot1 + geom_histogram(aes_(x = ~v, y = ~..density.., fill = "col1"), bins = 50)
  plot1 <- plot1 + scale_fill_manual(values = col)
  plot1 <- plot1 + guides(fill = FALSE)
  plot1 <- plot1 + ylab("density")
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Produce diagnostic plots of log-likelihood
#'
#' @description Produce diagnostic plots of the log-likelihood.
#'
#' @param project a MALECOT project, as produced by the function 
#'   \code{malecot_project()}
#' @param K which value of K to plot
#' @param rung which rung to plot. Defaults to the cold chain
#' @param col colour of the trace
#'
#' @export

plot_loglike_dignostic <- function(project, K = NULL, rung = NULL, col = "black") {
  
  # produce individual diagnostic plots and add features
  plot1 <- plot_trace(project, K = K, rung = rung, col = col)
  plot1 <- plot1 + ggtitle("MCMC trace")
  
  plot2 <- plot_acf(project, K = K, rung = rung, col = col)
  plot2 <- plot2 + ggtitle("autocorrelation")
  
  plot3 <- plot_density(project, K = K, rung = rung, col = col)
  plot3 <- plot3 + ggtitle("density")
  
  # produce grid of plots
  ret <- grid.arrange(plot1, plot2, plot3, layout_matrix = rbind(c(1,1), c(2,3)))
}

#------------------------------------------------
#' @title Plot prior on COI
#'
#' @description Produce plot of the prior on COI for given parameters. Options
#'   include the uniform distribution, and a modified form of Poisson and
#'   negative binomial distribution (see details).
#'
#' @details  The prior on COI can be uniform, Poisson, or negative binomial. In 
#'   the uniform case there is an equal chance of any given sample having a COI 
#'   between 1 and \code{COI_max} (inclusive). In the Poisson and negative 
#'   binomial cases it is important to note that the distribution is over 
#'   (COI-1), rather than over COI. This is because both Poisson and negative 
#'   binomial distributions allow for 0 values, which cannot be the case here 
#'   because observed samples must contain at least 1 genotype. Poisson and
#'   negative binomial distributions are also truncated at \code{COI_max}.
#'
#'   The full probability mass distribution for the Poisson case with 
#'   \code{COI_mean}\eqn{ = \mu} and \code{COI_max}\eqn{ = M} can be written
#'   \deqn{ Pr(COI = n) = z (\mu-1)^(n-1) exp(-(\mu-1)) / (n-1)! } where \eqn{z}
#'   is a normalising constant that ensures the distribution sums to unity, and
#'   is defined as: \deqn{1/z = \sum_{i=1}^M (\mu-1)^(i-1) exp(-(\mu-1)) /
#'   (i-1)! }
#'
#'   The mean of this distribution will generally be very close to \eqn{\mu},
#'   and the variance will be close to \eqn{\mu-1} (strictly it will approach
#'   these values as \eqn{M} tends to infinity).
#'
#'   The full probability mass distribution for the negative binomial case with 
#'   \code{COI_mean}\eqn{ = \mu}, \code{COI_dispersion}\eqn{ = v/\mu} and 
#'   \code{COI_max}\eqn{ = M} can be written \deqn{ Pr(COI = n) = z 
#'   \Gamma(n-1+N)/( \Gamma(N)(n-1)! ) p^N (1-p)^(n-1) } where \eqn{N =
#'   (\mu-1)^2/(v-\mu+1)}, \eqn{p = (\mu-1)/v}, and \eqn{z} is a normalising 
#'   constant that ensures the distribution sums to unity, and is defined as: 
#'   \deqn{1/z = \sum_{i=1}^M \Gamma(i-1+N)/( \Gamma(N)(i-1)! ) p^N (1-p)^(i-1) 
#'   }
#'
#'   The mean of this distribution will generally be very close to \eqn{\mu} and
#'   the variance will be close to \eqn{v} (strictly it will approach these 
#'   values as \eqn{M} tends to infinity).
#'
#' @param COI_model the type of prior on COI. Must be one of "uniform", 
#'   "poisson", or "nb" (negative binomial)
#' @param COI_mean the prior mean (before truncating at \code{COI_max}). Note 
#'   that this parameter only applies under the "poisson" and "nb" models
#' @param COI_dispersion  the ratio of the variance to the mean of the prior on
#'   COI. Only applies under the negative binomial model. Must be >1
#' @param COI_max the maximum COI allowed. Distributions are truncated at this
#'   value
#'
#' @export

plot_prior_COI <- function(COI_model = "poisson", COI_mean = 3, COI_dispersion = 2, COI_max = 20) {
  
  # check inputs
  assert_single_string(COI_model)
  assert_in(COI_model, c("uniform", "poisson", "nb"))
  assert_single_pos(COI_mean)
  assert_gr(COI_mean, 1)
  assert_single_pos(COI_dispersion)
  assert_gr(COI_dispersion, 1)
  assert_single_pos_int(COI_max, zero_allowed = FALSE)
  
  # create prior distribution
  COI <- 1:COI_max
  y <- switch(COI_model,
              "uniform" = {
                rep(1/COI_max, COI_max)
              },
              "poisson" = {
                ret <- dpois(COI-1, COI_mean-1)
                ret/sum(ret)
              },
              "nb" = {
                mu <- COI_mean
                v <- mu*COI_dispersion
                ret <- dnbinom(COI-1, size = (mu-1)^2/(v-mu+1), prob = (mu-1)/v)
                ret/sum(ret)
              })
  
  # produce plot
  plot1 <- ggplot(data.frame(COI = COI, probability = y)) + theme_bw()
  plot1 <- plot1 + geom_bar(aes_(x = ~COI, y = ~probability), width = 0.9, stat = "identity")
  plot1 <- plot1 + scale_x_continuous(limits = c(0, COI_max+1), expand = c(0,0)) + scale_y_continuous(limits = c(0, max(y)*1.1), expand = c(0,0))
  
  return(plot1)
}

#------------------------------------------------
#' @title Plot prior on allele frequencies
#'
#' @description Produce plot of the prior on COI for given parameters. This
#'   prior is Beta in the bi-allelic case, and Dirichlet in the multi-allelic
#'   case.
#'
#' @param lambda shape parameter(s) of the Beta or Dirichlet distribution. Can
#'   be a single scalar value, in which case the dimensionality is given by the
#'   number of \code{alleles}, or a vector of values specifying the shape
#'   parameter for each allele
#' @param alleles the dimensionality of the prior. Defaults to the length of
#'   \code{lambda}, or to 2 of \code{lambda} is a scalar
#'
#' @export

plot_prior_p <- function(lambda = 1, alleles = NULL) {
  
  # check inputs
  assert_pos(lambda)
  if (length(lambda) == 1) {
    if (is.null(alleles)) {
      alleles <- 2
    }
    assert_single_pos_int(alleles)
    lambda <- rep(lambda, alleles)
  }
  alleles <- length(lambda)
  
  # split between Beta and Dirichlet plots
  if (alleles == 2) {
    
    # create distribution
    x <- seq(0,1,l=1001)
    y <- dbeta(x, lambda[1], lambda[2])
    
    #produce plot
    plot1 <- ggplot(data.frame(x = x, y = y)) + theme_bw()
    plot1 <- plot1 + geom_area(aes_(x = ~x, y = ~y), colour = "black", fill = "dodgerblue3", alpha = 0.5)
    plot1 <- plot1 + scale_x_continuous(limits = c(0, 1), expand = c(0,0)) + scale_y_continuous(limits = c(0, max(y[is.finite(y)])*1.1), expand = c(0,0))
    plot1 <- plot1 + xlab("reference allele frequency") + ylab("probability density")
    return(plot1)
    
  } else {
    
    # random Dirichlet draws
    n <- 1e3
    dirichlet_draws <- mapply(function(x) {
                                rdirichlet(lambda)
                              }, x = rep(1,n), SIMPLIFY = FALSE)
    df <- data.frame(p = unlist(dirichlet_draws),
                     allele = rep(1:alleles, times = n),
                     sim = rep(1:n, each = alleles))
    
    #produce plot
    plot1 <- ggplot(df) + theme_bw()
    plot1 <- plot1 + geom_line(aes_(x = ~as.factor(allele), y = ~p, group = ~sim),
                               colour = "dodgerblue3", alpha = 0.1)
    plot1 <- plot1 + geom_vline(xintercept = 1:alleles, col = grey(0.6))
    plot1 <- plot1 + scale_y_continuous(limits = c(0, 1), expand = c(0,0))
    plot1 <- plot1 + xlab("allele") + ylab("probability density")
    return(plot1)
    
  }
}
