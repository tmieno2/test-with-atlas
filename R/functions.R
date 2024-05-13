#' test st_buffer()
#'
#' test
#' 
#' @param sf sf
#' @param radius radius
#' @returns sf
#' @import sf
#' @import data.table
#' @export
#' @examples
#' 
#' data(nc)
#' buffer <- make_buffer(sf::st_centroid(nc), 5)

make_buffer <- function(sf, radius) {
  sf::st_buffer(sf, radius)
}

#' test st_make_grid()
#' 
#' test
#'
#' @param sf sf
#' @param x_length long length
#' @param y_length lat length
#' @returns sf
#' @import sf
#' @export
#' @examples
#' 
#' data(nc)
#' buf <- make_buffer(sf::st_centroid(nc), 100)
#' grids <- make_grid(buf[1, ], 10, 20)
#' 
 
make_grid <- function(sf, x_length, y_length) {
  sf::st_make_grid(sf, cellsize = c(x_length, y_length))
}
