#' Example Structural Neuroimaging and Genetic Data for Spatial Model.
#'
#' Simulated data with 632 subjects, 486 SNPs from 24 structural neuroimaging measures.
#'
#' @docType data
#'
#' @usage data(sp_bgsmtr_example_data)
#'
#' @format A list with three components: "SNP_data", "SNP_groups", "BrainMeasures".
#' SNP_data is a 486-by-632 matrix containing minor allele counts for 632 subjects and 486 SNPs.
#' neighbourhood_structure is a 12 by 12 first order neighbourhood structure matrix.
#' BrainMeasures is a 24-by-632 matrix containing simulated volumetric and cortical thickness measures for 24 regions of interest.
#'
#'
#' @keywords datasets
#'
#' @examples
#' data(sp_bgsmtr_example_data)
#' names(sp_bgsmtr_example_data)
#' dim(sp_bgsmtr_example_data$SNP_data)
#' dim(sp_bgsmtr_example_data$BrainMeasures)
#' dim(sp_bgsmtr_example_data$neighbourhood_structure)
"sp_bgsmtr_example_data"
