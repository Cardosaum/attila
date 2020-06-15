library(tidyverse)

files <- dir("./data/parsed/")

all <- vector()
for (i in seq_along(files)){
    if (i == 1){
        print(files[i])
        tmp <- read_csv(file.path("./data/parsed", files[i]))
        all[i] <- tmp
        print(file.path("./data/parsed/", files[i]))
    }
}


all[1] <- read_csv(file.path("./data/parsed", files[i]))
