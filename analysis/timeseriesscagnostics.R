library(tidyverse)

#make file paths
startpath <- here::here("data-raw/timeseries/comp-engine-export-datapoints-{code}.csv")
timepaths <- glue::glue(startpath, code = c("1b6cc", "664b7", "c02e5"))

#load them
ts1 <- read_csv(timepaths[1])
ts2 <- read_csv(timepaths[2])
ts3 <- read_csv(timepaths[3])

#combine into a dataset
