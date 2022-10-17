#' @export
# PARAMETERS:
# x: list, containing the binding matrices (wide binding result objects) to be merged

MergeBindingResults <- function(x) {
  reshape2::acast(data = reshape2::melt(x), Var1 ~ Var2, value.var = "value", fun.aggregate = sum)
}
MergeBindingResults()
