#' ggvariant: Tidy, ggplot2-Native Visualization for Genomic Variants
#'
#' @docType package
#' @name ggvariant-package
#' @aliases ggvariant
#'
#' @importFrom stats aggregate ave
#' @importFrom utils tail
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggplot aes geom_col geom_point geom_segment
#'   geom_rect geom_text annotate scale_colour_manual scale_fill_manual
#'   scale_x_continuous scale_y_continuous labs coord_flip facet_wrap
#'   theme_minimal theme element_text element_line element_blank
#'   expansion margin .data
#' @importFrom scales comma percent_format
#' @importFrom cli cli_inform cli_warn cli_abort cli_progress_step cli_progress_done
"_PACKAGE"
