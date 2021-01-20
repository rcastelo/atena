## pretty-print names (private)
.pprintnames <- function(x) {
  y <- x
  if (length(x) > 2)
    y <- c(y[1], "...", y[length(y)])
  y <- paste(y, collapse=", ")
  y
}
