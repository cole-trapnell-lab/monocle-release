library(HSMMSingleCell)
library(monocle)
context("orderCells")

test_that("orderCells() correctly orders cells by State and Pseudotime", {
  set.seed(123)
  small_set <- load_HSMM_markers()
  
  small_set <- setOrderingFilter(small_set, row.names(fData(small_set)))
  
  small_set <- suppressMessages(reduceDimension(small_set, use_irlba=FALSE))
  
  small_set <- suppressMessages(orderCells(small_set))
  
  expect_false(is.null(pData(small_set)$Pseudotime))
  expect_false(is.null(pData(small_set)$State))
  
  expect_equivalent(levels(pData(small_set)$State), "1")
  
  expect_less_than(abs(max(pData(small_set)$Pseudotime) -  124.8243), 1e-4)
})

test_that("orderCells() properly validates its input",{
  
  # small_set <- load_HSMM_markers()
  
  # small_set <- setOrderingFilter(small_set, row.names(fData(small_set)))
  
  # results <- evaluate_promise(orderCells(small_set))
  
  # expect_equal(results$error, "Error: reducedDimS is NULL. Did you call reducedDimension()?")
  
  # small_set <- setOrderingFilter(small_set, row.names(fData(small_set)))
  # small_set <- suppressMessages(reduceDimension(small_set, use_irlba=FALSE))
  
  # results <- evaluate_promise(orderCells(small_set, num_paths=1000))
  
  # expect_equal(results$error, "Error: num_paths must be fewer than the number of cells")
  
  # results <- evaluate_promise(orderCells(small_set, num_paths=0))
  
  # expect_equal(results$error, "Error: num_paths must be > 0")
  
  # results <- evaluate_promise(orderCells(small_set, root_cell ="ARGHYBLARGH"))
  
  # expect_equal(results$error, "Error: root_cell must be the name of a cell in the CellDataSet")
  
  A <- c(2.0000, -81.5590)
  B <- c(58.0000, 161.9340)
  p <- c(42.0000, 41.8890)
  q <- project_point_to_line_segment(p, cbind(A, B))
  #project_point_to_line_segment()
  
  test <- Project2MST(HSMM_myo, project_point_to_line_segment)
  test <- orderCells(test, num_paths=1, root_state = NULL)
  plot_spanning_tree(test, color_by="Time", cell_size=2)
  plot_genes_in_pseudotime(test[46319, ], cell_size = 3)
  plot_spanning_tree(test, color_by="State", cell_size=2)
  unique(pData(test)$Pseudotime)

  Z <- test@reducedDimS
  P <- test@reducedDimK
  
  #   Z <- cds@reducedDimS
  #   P <- cds@reducedDimK
  Z <- as.data.frame(t(Z)); P <- as.data.frame(t(P))
  Z <- cbind(Z, cell = row.names(Z), type = 'Z')
  P <- cbind(P, cell = row.names(Z), type = 'P')
  
  
  ZP_df <- rbind(Z, P)
  ggplot(aes(x = V1, y = V2, group = cell, color = I(type)), 
         data = ZP_df) + geom_point(size = 3) + geom_line(aes(color = "blue")) + 
    xlab("C1") + ylab("C2") +  monocle_theme_opts()
  
  #new function
  test2 <- orderCells(HSMM_myo, num_paths=1, root_state = NULL)
  plot_spanning_tree(test2, color_by="Time", cell_size=2)
  plot_genes_in_pseudotime(test2[46319, ], cell_size = 3)
  plot_spanning_tree(test2, color_by="State", cell_size=2)
  unique(pData(test2)$Pseudotime)
  
  Z2 <- test2@reducedDimS
  P2 <- test2@reducedDimK
  
  #   Z <- cds@reducedDimS
  #   P <- cds@reducedDimK
  Z2 <- as.data.frame(t(Z2)); P2 <- as.data.frame(t(P2))
  Z2 <- cbind(Z2, cell = row.names(Z2), type = 'Z')
  P2 <- cbind(P2, cell = row.names(Z2), type = 'P')
  
  
  ZP2_df <- rbind(Z2, P2)
  ggplot(aes(x = V1, y = V2, group = cell, color = I(type)), 
         data = ZP2_df) + geom_point(size = 3) + geom_line(aes(color = "blue")) + 
    xlab("C1") + ylab("C2") +  monocle_theme_opts()
  
})