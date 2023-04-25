############################
# create some test objects
############################

tree_a <- read.beast("mtGeo.MCC.txt")
tree_b <- read.beast("yGeo.MCC.txt")
tree_c <- read.beast("autoGeo.MCC.txt")

trees <- c(tree_a,tree_b,tree_c)
class(trees) <- "multiPhylo"

############################
# MRCA tests
############################

test_that("mrcaDist of mtGeo and yGeo is correct", {
	expect_equal(mrcaDist(tree_a,tree_b, coord_name="geo", continuous=TRUE)[1],6.118718447700729)
	})

test_that("mrcaDist not sensitive to permuting inputs", {
	expect_equal(mrcaDist(tree_a,tree_b, coord_name="geo", continuous=TRUE),mrcaDist(tree_b,tree_a, coord_name="geo", continuous=TRUE))
	})

test_that("mrcaDist same when computed from mrcaVec", {
	expect_equal(mrcaDist(tree_a,tree_b, coord_name="geo", continuous=TRUE),mrcaDist(mrcaVec(tree_a, coord_name="geo", continuous=TRUE),mrcaVec(tree_b, coord_name="geo", continuous=TRUE), coord_name="geo", continuous=TRUE))
	})

test_that("save_memory version of multiMrcaDist equals normal multiMrcaDist", {
	expect_equal(multiMrcaDist(trees, coord_name="geo", continuous=TRUE),multiMrcaDist(trees, coord_name="geo", continuous=TRUE,save.memory=TRUE))
	})

test_that("multiMrcaDist of mtGeo, yGeo, and autoGeo is correct", {
	expect_equal(multiMrcaDist(trees,coord_name="geo",continuous=TRUE)[1:3],c(6.118718447700729257122,4.894385071111685192591,4.555614265727615297408))
	})

############################
# MASPG tests
############################

test_that("maspgDist of mtGeo and yGeo is correct", {
	expect_equal(maspgDist(tree_a,tree_b, coord_name="geo", discretization=1)[1],0.875)
	})

test_that("maspgDist not sensitive to permuting inputs", {
	expect_equal(maspgDist(tree_a,tree_b, coord_name="geo", discretization=1),maspgDist(tree_b,tree_a, coord_name="geo", discretization=1))
	})

test_that("multiMaspgDist of mtGeo, yGeo, and autoGeo is correct", {
	expect_equal(multiMaspgDist(trees, coord_name="geo", discretization=1)[1:3],c(0.875,0.750,0.875))
	})

############################
# test for errors and warmings
############################

test_that("error is given if input is not of class phylo / multiphylo", {
  expect_error(mrcaVec(trees))
  expect_error(mrcaDist(trees))
  expect_error(multiMrcaDist(tree_a))
  expect_error(maspgVec(trees))
  expect_error(maspgDist(trees))
  expect_error(multiMaspgDist(tree_a))
  })

test_that("error is given if trees have different tip labels", {
  tree_d <- tree_a
  tree_d@phylo$tip.label <- 1:16
  expect_error(mrcaDist(tree_a,tree_d))
  expect_error(maspgDist(tree_a,tree_d))
  })