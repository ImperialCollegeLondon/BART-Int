library(ggplot2)
require(gridExtra)
library(ggrepel)
library(plotrix)
`%--%` <- function(x, y) 
  # from stack exchange:
  # https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}

getTreeNumLeaf <- function(sampler, chainNum, sampleNum, treeNum)
  
  #'Build tree
  #' 
  #'@description The function builds tree and fills attributes to each node of the tree
  #' 
  #'@param sampler List; Bart model
  #'@param chainNum Integer; The index of the chain 
  #'@param sampleNum Integer; The index of the posterior sample draw
  #'@param treenNum Integer; The index of the tree sample extracted
  #'
  #'@return A tree sample in BART fitting at the position specifized
#'by chainNum, sampleNUm and treeNum, with attributeds filled into
#'each node of the tree

{
  cutPoints <- dbarts:::createCutPoints(sampler)
  if (sampler$control@keepTrees) {
    treeString <- sampler$state[[chainNum]]@savedTrees[treeNum, sampleNum]
    treeFits <- sampler$state[[chainNum]]@savedTreeFits[, treeNum, sampleNum]
  }
  else {
    treeString <- sampler$state[[chainNum]]@trees[treeNum]
    treeFits <- sampler$state[[chainNum]]@treeFits[,treeNum]
  }                           
  
  tree <- dbarts:::buildTree(strsplit(gsub("\\.", "\\. ", treeString),
                                      " ", fixed = TRUE)[[1]])
  tree$remainder <- NULL
  
  tree$indices <- seq_len(length(sampler$data@y))
  tree <- dbarts:::fillObservationsForNode(tree, sampler, cutPoints)
  
  tree <- dbarts:::fillPlotInfoForNode(tree, sampler, treeFits)
  maxDepth <- dbarts:::getMaxDepth(tree)
  
  tree <- dbarts:::fillPlotCoordinatesForNode(tree, maxDepth, 1L, 1L)
  
  return (tree)
}

.simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

# pdf("Figures/2.pdf" %--% c(genzFunctionName), width=11, height=2.8)
pdf("Figures/2.pdf" %--% c(genzFunctionName), width=5, height=5)
# par(pty="s", mfrow=c(1, 3), mai = c(1, 0.4, 0.4, 0.4), cex=1.3)
par(pty="s", cex=1.3)
for (whichGenz in c(7)) {
  # Find Genz function
  if (whichGenz == 1) { genzFunctionName <-  "cont" }
  if (whichGenz == 2) { genzFunctionName <-  "copeak" }
  if (whichGenz == 3) { genzFunctionName <-  "disc" }
  if (whichGenz == 4) { genzFunctionName <-  "gaussian" }
  if (whichGenz == 5) { genzFunctionName <-  "oscil" }
  if (whichGenz == 6) { genzFunctionName <-  "prpeak" }
  if (whichGenz == 7) { genzFunctionName <-  "step" }
  
  load("results/genz/%s/%sDim1Uniform_1.RData" %--% c(whichGenz, genzFunctionName))
  num_leaves_1 <- c()
  for (t in 1:50) {
    num_leaves_1 <- c(num_leaves_1, getTreeNumLeaf(results_models$BART$model$fit, 1, 300, t)$index-1)
  }
  load("results/genz/%s/%sDim10Uniform_1.RData" %--% c(whichGenz, genzFunctionName))
  num_leaves_10 <- c()
  for (t in 1:50) {
    num_leaves_10 <- c(num_leaves_10, getTreeNumLeaf(results_models$BART$model$fit, 1, 300, t)$index-1)
  }
  df_num_leaves <- data.frame(num_leaves_1=num_leaves_1, num_leaves_10=num_leaves_10) 
  
  multhist(list(df_num_leaves$num_leaves_1, df_num_leaves$num_leaves_10), 
           breaks=seq(0.5, 8.5), names.arg=1:8, col=c("orangered", "orangered4"), 
           xlab="K", ylab="Frequency", 
           main=.simpleCap(genzFunctionName),
           ylim=c(0,50), xlim=c(1,21))
  if (whichGenz %in% c(1, 7)) {
    legend("topright", c("d=1", "d=10"), fill=c("orangered", "orangered4"), bty="n")
  }
}
dev.off()

