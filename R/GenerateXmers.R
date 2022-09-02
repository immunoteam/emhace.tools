#' @export
GenerateXmers <- function(x, l, string_input = TRUE) {
  if(string_input) {x <- unlist(strsplit(x, ""), use.names = FALSE)}
  coordsmtx <- cbind(seq(1, length(x)), seq(1, length(x)) + l - 1)
  coordsmtx <- coordsmtx[coordsmtx[, 1] <= length(x), , drop = FALSE]
  coordsmtx <- coordsmtx[coordsmtx[, 2] <= length(x), , drop = FALSE]
  tmp <- apply(apply(coordsmtx, 1, function(y) {x[y[1]:y[2]]}), 2, function(y) {paste0(y, collapse = "")})
  tmp[nchar(tmp) == l]
}
