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

data("timeOmics.simdata")
sim.data <- timeOmics.simdata$sim

dim(sim.data)


remove.low.cv <- function(X, cutoff = 0.5){
  # var.coef
  cv <- unlist(lapply(as.data.frame(X),
                      function(x) abs(sd(x)/mean(x))))
  return(X[,cv > cutoff])
}

data.filtered <- remove.low.cv(sim.data, 0.5)

devtools::install_github("cran/lmms", force = TRUE)
library(lmms)


time <- timeOmics.simdata$time
head(time)

lmms.output <- lmms::lmmSpline(data = data.filtered, time = time,
                               sampleID = rownames(data.filtered), deri = FALSE,
                               basis = "p-spline", numCores = 4, timePredict = 1:9,
                               keepModels = TRUE)
modelled.data <- t(slot(lmms.output, 'predSpline'))

lmms.output <- timeOmics.simdata$lmms.output
modelled.data <- timeOmics.simdata$modelled

# gather data
data.gathered <- modelled.data %>% as.data.frame() %>%
  rownames_to_column("time") %>%
  mutate(time = as.numeric(time)) %>%
  pivot_longer(names_to="feature", values_to = 'value', -time)

# plot profiles
ggplot(data.gathered, aes(x = time, y = value, color = feature)) + geom_line() +
  theme_bw() + ggtitle("`lmms` profiles") + ylab("Feature expression") +
  xlab("Time")


filter.res <- lmms.filter.lines(data = data.filtered,
                                lmms.obj = lmms.output, time = time)
profile.filtered <- filter.res$filtered


conflicts_prefer(mixOmics::pca)
# run pca
pca.res <- pca(X = profile.filtered, ncomp = 5, scale=FALSE, center=FALSE)

# tuning ncomp
pca.ncomp <- getNcomp(pca.res, max.ncomp = 5, X = profile.filtered,
                      scale = FALSE, center=FALSE)

pca.ncomp$choice.ncomp
plot(pca.ncomp)

# final model
pca.res <- pca(X = profile.filtered, ncomp = 2, scale = FALSE, center=FALSE)

# extract clhttp://127.0.0.1:19297/graphics/plot_zoom_png?width=1229&height=641uster
pca.cluster <- getCluster(pca.res)
head(pca.cluster)


summary(pca.res)

pca.res$rotation[1:20,1:2]

install.packages("ggplot2")
library(ggplot2)
library(mixOmics)
plotIndiv(pca.res)
plotVar(pca.res)
plotLoadings(pca.res)
