#' ggvariant: Tidy, ggplot2-Native Visualization for Genomic Variants
#'
#' @docType package
#' @name ggvariant-package
#' @aliases ggvariant
#'
#' @references
#' Danecek P, Auton A, Abecasis G, et al.; 1000 Genomes Project Analysis
#' Group (2011). The variant call format and VCFtools. *Bioinformatics*,
#' 27(15), 2156-2158. \doi{10.1093/bioinformatics/btr330}
#'
#' Alexandrov LB, Kim J, Haradhvala NJ, et al.; PCAWG Consortium (2020).
#' The repertoire of mutational signatures in human cancer. *Nature*,
#' 578(7793), 94-101. \doi{10.1038/s41586-020-1943-3}
#'
#' @importFrom stats aggregate ave setNames na.omit
#' @importFrom utils tail head
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggplot aes geom_col geom_point geom_segment geom_tile
#'   geom_rect geom_text annotate scale_colour_manual scale_fill_manual
#'   scale_x_continuous scale_x_discrete scale_y_continuous labs coord_flip
#'   facet_wrap theme_minimal theme element_text element_line element_blank
#'   expansion margin .data
#' @importFrom scales comma percent_format
#' @importFrom cli cli_inform cli_warn cli_abort cli_progress_step cli_progress_done
"_PACKAGE"
