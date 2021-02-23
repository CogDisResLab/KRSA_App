#' TODO
#'
#' TODO
#'
#' @param data krsa data_modeled$scaled
#' @param samples sample names
#' @param peptides vector of peptides
#'
#' @return vector
#'
#' @import dplyr
#' @import rlang
#'
#' @export
#'
#' @examples
#' TRUE

krsa_violin_plot_2 <- function(data, peptides,facet = T, facet_factor,samples = NULL, groups = NULL) {
  
  data %>% dplyr::filter(Peptide %in% peptides) %>%
    {if(!is.null(samples)) filter(.,SampleName %in% samples) else .} %>%
    {if(!is.null(groups)) filter(.,Group %in% groups) else .} %>%
    ggplot(aes(SampleName, slope)) + 
    geom_violin(aes(fill = Group), show.legend = F, trim = F, width=0.4) +
    geom_boxplot(width=0.1, fill="white") +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=1, alpha = 1/2) +
    geom_line(aes(group = Peptide), alpha = 1/2) +
    labs(
      x = "",
      y = "Signal Intensity"
    ) + theme(
      axis.text.x = element_text(size = 6)
      
    ) + theme_bw() -> gg
  
  {if(facet == T) facet_wrap(facet_factor, scales = "free", nrow = 2, ncol = 2)}
  if(facet == T) {gg + facet_wrap(facet_factor, scales = "free", nrow = 2, ncol = 2)}
  else {gg}
  
}
