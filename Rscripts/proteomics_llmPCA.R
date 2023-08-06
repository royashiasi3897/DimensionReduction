pacman::p_load(conflicted,
               wrappedtools,
               tidyverse,
               MLSeq,
               ggfortify, GGally,
               PCAtools, # bioconductor
               FactoMineR)

# conflict_scout()
#conflicts_prefer(dplyr::filter())
conflict_prefer('slice','dplyr')
conflict_prefer("filter", "dplyr")
conflict_prefer('screeplot','stats')
conflict_prefer('biplot','stats')
#pacman::p_load(lmms.filter.lines)



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

## install timeOmics
BiocManager::install('timeOmics', force = TRUE)

install.packages("devtools")
# then load
library(devtools)
install_github("abodein/timeOmics")

library(timeOmics)
library(tidyverse)


library(readxl)
rawdata <- read_csv("F:/bioinformatic_weiterbildung_final/internship/stage_2/simulated_data/example_data.csv")



rawdata <- rawdata |>
  pivot_wider(names_from = metabolites,
              values_from = values)

#view(rawdata)
#saveRDS(rawdata,'Data/cervical.RData')


predvars <- FindVars(c('metabolite'))

###Data preprocessing

remove.low.cv <- function(X, cutoff = 0.5){
  # var.coef
  cv <- unlist(lapply(as.data.frame(X),
                      function(x) abs(sd(x)/mean(x))))
  return(X[,cv > cutoff])
}

data.filtered <- remove.low.cv(rawdata |> select(predvars$names), 0.5)
conflicts_prefer(dplyr::select)
devtools::install_github("cran/lmms", force = TRUE)
library(lmms)
ls("package:lmms")


time <- rawdata$time
head(time)

lmms.output <- lmms::lmmSpline(data = data.filtered, time = time,
                               sampleID = rownames(data.filtered), deri = FALSE,
                               basis = "p-spline", numCores = 4, timePredict = 1:9,
                               keepModels = TRUE)
modelled.data <- t(slot(lmms.output, 'predSpline'))
#lmms.output <- timeOmics.simdata$lmms.output
#modelled.data <- timeOmics.simdata$modelled

data.gathered <- modelled.data %>% as.data.frame() %>%
  rownames_to_column("time") %>%
  mutate(time = as.numeric(time)) %>%
  pivot_longer(names_to="feature", values_to = 'value', -time)

# plot profiles
ggplot(data.gathered, aes(x = time, y = value, color = feature)) + geom_line() +
  theme_bw() + ggtitle("lmms profiles") + ylab("Feature expression") +
  xlab("Time")

filter.res <- lmms.filter.lines(data = data.filtered,
                                lmms.obj = lmms.output, time = time)
profile.filtered <- filter.res$filtered

###pca
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mixOmics")
library(mixOmics)

pca.res <- pca(X = profile.filtered, ncomp = 5, scale=FALSE, center=FALSE)

##########andreas
#pca_out <- prcomp(rawdata |> select(predvars$names),
#                  center = T,scale. = T)
#summary(pca_out)
#screeplot(pca_out,npcs = 5)




################################################################
