#' @importFrom tidytree as_tibble
#' @importFrom tidytree child
#' @importFrom tidytree offspring
PairingFunction <- function(x) {
  return(.5*(x[1,]+x[2,])*(x[1,]+x[2,]+1)+x[2,])
}

mychild <- function(myphylo, mynode) {unlist(child(as_tibble(myphylo), mynode)[, 2])}
myoffspring <- function(myphylo, mynode) {unlist(offspring(as_tibble(myphylo), mynode)[, 2])}


#' MASPG metric function
#'
#' Comparison of two phylogeographic trees using the Maximum Agreement Sub-Phylogeography (MASPG) incompatibility metric.
#'
#' @author Benjamin John Singer \email{benjamin.singer@@bnc.ox.ac.uk}
#' @author Antonello Di Nardo \email{antonello.dinardo@@pirbright.ac.uk}
#' @author Luca Ferretti \email{luca.ferretti@@gmail.com}
#'
#' @param tree.a an object of the class \code{treedata}
#' @param tree.b an object of the class \code{treedata} (with the same tip labels as tree.a).
#' @param coord_name (optional) a string specifying the attribute name used to indicate the geographic trait variable in the \code{treedata} object.
#' @param discretization (optional) I_MASPG is a discrete measure which compares if locations are the same or different, and which does not use geographic distances. To use continuous distances a discretization must be specified. This is a numeric measure of distance under which two points are 'the same'. If discretization is used, 2-dimensional co-ordinates are assumed, with names coord_name1 and coord_name2.
#' @param geography (optional) If false, then the Maximum Agreement Sub-Tree (MAST) metric is computed instead of the MASPG metric
#' 
#' @return The distance between the two trees according to the metric
#' 
#' @importFrom tidytree rootnode
#'
#' @export
maspgDist <- function(tree.a, tree.b, coord_name="Location", discretization=FALSE, geography=TRUE) {
  #Vector of node numbering
  NNODES <- tree.a@phylo$Nnode*2+1
  nodes1 <- as.numeric(unlist(tree.a@data$node))
  nodes2 <- as.numeric(unlist(tree.b@data$node))
  #Vector of tip labels
  tips1 <- tree.a@phylo$tip.label
  tips2 <- tree.b@phylo$tip.label
  #Vector of locations
  if (!discretization) {
    if ((paste(coord_name, "1", sep='') %in% colnames(tree.a@data)) & !(coord_name %in% colnames(tree.a@data))) {
      stop("Continuous spatial co-ordinates require a discretization threshold.")
    }
    locs1 <- unlist(tree.a@data[coord_name])
    locs2 <- unlist(tree.b@data[coord_name])
  } else {
    locs1x <- unlist(tree.a@data[paste(coord_name, "1", sep='')])
    locs1y <- unlist(tree.a@data[paste(coord_name, "2", sep='')])
    locs2x <- unlist(tree.b@data[paste(coord_name, "1", sep='')])
    locs2y <- unlist(tree.b@data[paste(coord_name, "2", sep='')])
    locs1 <- array(c(locs1x, locs1y), dim=c(length(locs1x), 2))
    locs2 <- array(c(locs2x, locs2y), dim=c(length(locs2x), 2))
  }
  #Find roots
  root1 <- as.vector(unlist(rootnode(as_tibble(tree.a@phylo))[, 2]))
  root2 <- as.vector(unlist(rootnode(as_tibble(tree.b@phylo))[, 2]))
  
  #Get number of nodes
  N <- length(unlist(tree.a@data$node))
  #and number of tips
  NTIPS <- length(unlist(tree.a@phylo$tip.label))
  
  #Initialise score matrix
  score.matrix <- array(rep(0, N), dim=c(N, N))-1

  #Perform node-node comparisons
  score.matrix[nodes1<=NTIPS, nodes2<=NTIPS] <- array(rep(tips1[nodes1[nodes1<=NTIPS]], NTIPS),
                                                      dim=c(NTIPS, NTIPS))==aperm(array(rep(tips2[nodes2[nodes2<=NTIPS]], NTIPS),
                                                                                        dim=c(NTIPS, NTIPS)), c(2, 1))
  
  #z just helps the iteration work using floor/modulus method
  z <- c(NNODES, 1:(NNODES-1))
  
  #Perform the calculations! The way the array of comparisons to make is set up,
  #everything is calculated in the required order just by going through
  #the array in order
  if (!discretization) {
    for (i in 1:(NNODES^2)) {
      node1 <- nodes1[ceiling(i/NNODES)]
      node2 <- nodes2[z[i%%NNODES+1]]
      if ((node1<=NTIPS) & (node2<=NTIPS)) {
        next
      } else if (node1<=NTIPS) {
        score.matrix[nodes1==node1, nodes2==node2] <- ((nodes1[nodes1<=NTIPS])[tips2==tips1[(nodes1[nodes1<=NTIPS])==node1]] %in%
                                                        myoffspring(tree.b@phylo, node2))
      } else if (node2<=NTIPS) {
        score.matrix[nodes1==node1, nodes2==node2] <- ((nodes2[nodes2<=NTIPS])[tips1==tips2[(nodes2[nodes2<=NTIPS])==node2]] %in%
                                                        myoffspring(tree.a@phylo, node1))
      } else if ((locs1[nodes1==node1] != locs2[nodes2==node2]) & geography) {
        score.matrix[nodes1==node1,nodes2==node2] <- max(score.matrix[nodes1==node1, nodes2==mychild(tree.b@phylo, node2)[1]],
                                                         score.matrix[nodes1==node1, nodes2==mychild(tree.b@phylo, node2)[2]],
                                                         score.matrix[nodes1==mychild(tree.a@phylo, node1)[1], nodes2==node2],
                                                         score.matrix[nodes1==mychild(tree.a@phylo, node1)[2], nodes2==node2])
        if (is.na(score.matrix[nodes1==node1, nodes2==node2])) {
          warning("Non-numerical entry in the score matrix")
        }
      } else {
        mychild.score1 <- score.matrix[nodes1==mychild(tree.a@phylo, node1)[1], nodes2==mychild(tree.b@phylo, node2)[2]] +
                          score.matrix[nodes1==mychild(tree.a@phylo, node1)[2], nodes2==mychild(tree.b@phylo, node2)[1]]
        mychild.score2 <- score.matrix[nodes1==mychild(tree.a@phylo, node1)[1], nodes2==mychild(tree.b@phylo, node2)[1]] +
                          score.matrix[nodes1==mychild(tree.a@phylo, node1)[2], nodes2==mychild(tree.b@phylo, node2)[2]]
        score.matrix[nodes1==node1, nodes2==node2] <- max(mychild.score1, mychild.score2,
                                                          score.matrix[nodes1==node1, nodes2==mychild(tree.b@phylo,node2)[1]],
                                                          score.matrix[nodes1==node1, nodes2==mychild(tree.b@phylo,node2)[2]],
                                                          score.matrix[nodes1==mychild(tree.a@phylo, node1)[1], nodes2==node2],
                                                          score.matrix[nodes1==mychild(tree.a@phylo, node1)[2], nodes2==node2])
      }
    }
  } else {
    for (i in 1:(NNODES^2)) {
      node1 <- nodes1[ceiling(i/NNODES)]
      node2 <- nodes2[z[i%%NNODES+1]]
      if ((node1<=NTIPS) & (node2<=NTIPS)) {
        next
      } else if (node1 <= NTIPS) {
        score.matrix[nodes1==node1, nodes2==node2] <- ((nodes1[nodes1<=NTIPS])[tips2==tips1[(nodes1[nodes1<=NTIPS])==node1]] %in%
                                                        myoffspring(tree.b@phylo, node2))
      } else if (node2 <= NTIPS) {
        score.matrix[nodes1==node1, nodes2==node2] <- ((nodes2[nodes2<=NTIPS])[tips1==tips2[(nodes2[nodes2<=NTIPS])==node2]] %in%
                                                        myoffspring(tree.a@phylo, node1))
      } else if ((dist(rbind(locs1[nodes1==node1, ], locs2[nodes2==node2, ])) > discretization) & geography) {
        score.matrix[nodes1==node1, nodes2==node2] <- max(score.matrix[nodes1==node1, nodes2==mychild(tree.b@phylo, node2)[1]],
                                                          score.matrix[nodes1==node1, nodes2==mychild(tree.b@phylo, node2)[2]],
                                                          score.matrix[nodes1==mychild(tree.a@phylo, node1)[1], nodes2==node2],
                                                          score.matrix[nodes1==mychild(tree.a@phylo, node1)[2], nodes2==node2])
        if (is.na(score.matrix[nodes1==node1, nodes2==node2])) {
          warning("Non-numerical entry in the score matrix")
        }
      } else {
        mychild.score1 <- score.matrix[nodes1==mychild(tree.a@phylo, node1)[1], nodes2==mychild(tree.b@phylo, node2)[2]] +
                          score.matrix[nodes1==mychild(tree.a@phylo, node1)[2], nodes2==mychild(tree.b@phylo, node2)[1]]
        mychild.score2 <- score.matrix[nodes1==mychild(tree.a@phylo, node1)[1], nodes2==mychild(tree.b@phylo, node2)[1]] +
                          score.matrix[nodes1==mychild(tree.a@phylo, node1)[2], nodes2==mychild(tree.b@phylo, node2)[2]]
                          score.matrix[nodes1==node1, nodes2==node2] <- max(mychild.score1, mychild.score2,
                                                                            score.matrix[nodes1==node1, nodes2==mychild(tree.b@phylo, node2)[1]],
                                                                            score.matrix[nodes1==node1, nodes2==mychild(tree.b@phylo, node2)[2]],
                                                                            score.matrix[nodes1==mychild(tree.a@phylo, node1)[1], nodes2==node2],
                                                                            score.matrix[nodes1==mychild(tree.a@phylo, node1)[2], nodes2==node2])
      }
    }
  }
  
  return((NTIPS-score.matrix[nodes1==root1, nodes2==root2])/NTIPS)
}


#' Metric function for \code{multiPhylo} input
#'
#' Comparison of a list of phylogeographic trees using the Maximum Agreement Sub-Phylogeography (MASPG) incompatibility metric. Output is given as a pairwise distance matrix.
#'
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#' @author Benjamin John Singer \email{benjamin.singer@@bnc.ox.ac.uk}
#' @author Antonello Di Nardo \email{antonello.dinardo@@pirbright.ac.uk}
#' @author Luca Ferretti \email{luca.ferretti@@gmail.com}
#'
#' @param trees an object of the class \code{multiPhylo} containing the trees to be compared
#' @param coord_name (optional) a string specifying the attribute name used to indicate the geographic trait variable in the \code{treedata} object.
#' @param discretization (optional) I_MASPG is a discrete measure which compares if locations are the same or different, and which does not use geographic distances. To use continuous distances a discretization must be specified. This is a numeric measure of distance under which two points are 'the same'. If discretization is used, 2-dimensional co-ordinates are assumed, with names coord_name1 and coord_name2.
#' @param geography (optional) If false, then the Maximum Agreement Sub-Tree (MAST) metric is computed instead of the MASPG metric.
#' 
#' @return The pairwise tree incompatibility matrix
#'
#'
#' @export
multiMaspgDist <- function(trees, coord_name="Location", discretization=FALSE, geography=TRUE) {
  if (!inherits(trees, "multiPhylo")) stop("trees should be a multiphylo object")
  num_trees <- length(trees) 
  if (num_trees<2) {
    stop("multiMaspgDist expects at least two trees")
  }
  
  # make name labels well defined
  if (is.null(names(trees))) names(trees) <- 1:num_trees 
  else if (length(unique(names(trees))) != num_trees) {
    warning("duplicates detected in tree labels - using generic names")
    names(trees) <- 1:num_trees
  }
  lab <- names(trees)
  
  # check all trees have same tip labels
  for (i in 1:num_trees) {
    if (!setequal(trees[[i]]@phylo$tip.label, trees[[1]]@phylo$tip.label)) {
      stop(paste0("Tree ", lab[[i]], " has different tip labels from the first tree."))
    } 
  }
  
  distances <- matrix(0.0, num_trees, num_trees)
  
  sapply(1:(num_trees-1), function(i) {
    sapply((i+1):num_trees, function(j) {
      distances[i, j] <<- distances[j, i] <<- maspgDist(trees[[i]], trees[[j]], coord_name, discretization, geography)
    })
  })
  
  return(as.dist(distances))
}
