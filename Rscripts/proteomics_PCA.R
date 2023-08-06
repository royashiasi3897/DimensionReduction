pacman::p_load(conflicted,
               wrappedtools,
               tidyverse,
               ggfortify, GGally,
               PCAtools, # bioconductor
               FactoMineR)

# conflict_scout()
conflict_prefer('slice','dplyr')
conflict_prefer("filter", "dplyr")
conflict_prefer('screeplot','stats')
conflict_prefer('biplot','stats')

##read
rawdata <- readRDS('Data/cervical.RData')
predvars <- FindVars(c('-'))
#rawdata <- mutate(rawdata,
#                  across(all_of(predvars$names),
#                         .fns = ~log(.x+0.01)))

pca_out <- prcomp(rawdata |> select(predvars$names),
                  center = T,scale. = T)
summary(pca_out)

screeplot(pca_out,npcs = 5)

pca_out$rotation[1:10,1:5]