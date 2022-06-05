plasma_data <- read.table("/Users/karekv/GDrive/Research/finite/Data/littmann_data.txt",
           na.strings = c("", " ", "NA"),
           header = TRUE,
           strip.white = FALSE,
           sep = "\t",
           colClasses = c("integer", "factor", "numeric", "factor"))


plasma_data$lpa_cat <- cut(plasma_data$lpa, breaks = c(0, 10, 30, 120, Inf), include.lowest = T)
