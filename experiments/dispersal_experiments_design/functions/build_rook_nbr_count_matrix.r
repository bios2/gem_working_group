# build_rook_nbr_count_matrix.r
# Returns an n_row x n_col integer matrix where each value is the number of
# rook-adjacent (shared-boundary) neighbours for that cell.
# Interior cells = 4, edge cells = 3, corner cells = 2.

build_rook_nbr_count_matrix <- function(n_row, n_col) {
  m          <- matrix(4L, nrow = n_row, ncol = n_col)
  m[1, ]     <- m[1, ]     - 1L
  m[n_row, ] <- m[n_row, ] - 1L
  m[, 1]     <- m[, 1]     - 1L
  m[, n_col] <- m[, n_col] - 1L
  m
}
