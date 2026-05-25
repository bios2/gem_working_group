# build_rook_nbr_count_matrix.r
# Returns an n_row x n_col integer matrix where each value is the number of
# rook-adjacent (shared-boundary) neighbours for that cell.
# Interior cells = 4, edge cells = 3, corner cells = 2.

build_rook_nbr_count_matrix <- function(n_row, n_col) {
  # start every cell at 4 (N, S, W, E all available for interior cells)
  m          <- matrix(4L, nrow = n_row, ncol = n_col)
  # top and bottom rows have no northern / southern neighbour respectively
  m[1, ]     <- m[1, ]     - 1L
  m[n_row, ] <- m[n_row, ] - 1L
  # left and right columns have no western / eastern neighbour respectively
  m[, 1]     <- m[, 1]     - 1L
  m[, n_col] <- m[, n_col] - 1L
  return(m)
}
