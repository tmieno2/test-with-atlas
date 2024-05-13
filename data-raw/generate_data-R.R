## code to prepare `generate_data.R` dataset goes here

library(sf)
nc <-
  st_read(system.file("shape/nc.shp", package = "sf")) %>%
  st_transform(26917)

st_crs(nc)$wkt <- gsub("°|º", "", st_crs(nc)$wkt)

usethis::use_data(nc, overwrite = TRUE)

data(nc)
buf <- make_buffer(sf::st_centroid(nc), 100)
plot(buf)

make_grid(buf[1, ], 10, 20)

#++++++++++++++++++++++++++++++++++++
#+
#++++++++++++++++++++++++++++++++++++
library(ofpetrial)

data(td_single_input)

usethis::use_data(td_single_input, overwrite = TRUE)
