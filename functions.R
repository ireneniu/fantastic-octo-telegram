process_image <- function(image_file_name, k_list){
  ## process_image(image_file_name, k_list) conducts k-means clustering 
  ## on an image and produces clustering information including tidied clusters, 
  ## associated RGB values and DMC colour information and the original data frame in a list
  ##
  ## Input:
  ## - image_file_name : The name of a PNG or JPEG image. Example: "x.jpeg"
  ## - k_list : The list of numbers of centres in the clustering. 
  ##   Example: c(2:10), or c(2,4,6)
  ##
  ## Output:
  ## - cluster_info: A list of information derived from the k_means 
  ##  that includes clustering information and the original data frame. 
  ##  The first element is a tibble of clustering information containing the original 
  ##  output of the kclust calls, the tidied clusters, their associated RGB values and 
  ##  their nearest DMC thread colour.
  ##  The second element is the original data frame for further use in functions
  ## 
  ## Example:
  ##   process_image("liza.jpg",c(2,4,6))
  ##
  #install packages 
  library(imager)
  library(dplyr)
  library(tidymodels)
  
  #devtools::install_github("sharlagelfand/dmc")
  library(dmc)
  
  im <- imager::load.image(image_file_name)
  plot(im)
  tidy_dat <- as.data.frame(im, wide="c") %>% rename(R = c.1, G = c.2, B = c.3)
  
  dat <- select(tidy_dat,c(-x,-y))
  kclusts <-
    tibble(k = k_list) %>%
    mutate(
      kclust = map(k, ~kmeans(x = dat , centers = .x, nstart=4)),
      glanced = map(kclust, glance),
      tidied = map(kclust, tidy)
      #augmented = map(kclust, augment, tidy_dat)
    )
  
  #Loop through all the tidied datasets in kclusts and attach rgb information using rgb function
  
  
  for (i in (1:length(kclusts$tidied))){
    #create empty tibble
    dmc <- tibble(
      dmc = character(),
      name = character(),
      hex = character(),
      red = numeric(),
      green = numeric(),
      blue = numeric()
    )
    kclusts$tidied[[i]] <- kclusts$tidied[[i]] %>% mutate(col = rgb(R,G,B))
    for (j in (1:length(kclusts$tidied[[i]]$col))){
      dmc_col <- dmc(kclusts$tidied[[i]]$col[[j]])
      dmc <- dmc %>% add_row(dmc_col)
    } 
    kclusts$tidied[[i]] <- kclusts$tidied[[i]] %>% mutate(dmc)  
    
  }
  
  #keep original data frame information 
  kclusts <- list(kclusts,tidy_dat)
  return(kclusts)
}


#=========================================================================================
#=========================================================================================

scree_plot <- function(cluster_info){
  ## scree_plot(cluster_info) produces and plots a ratio scree plot based on 
  ## the cluster_info output from process_image function
  ##
  ## Input:
  ## - cluster_info: cluster information output from process_image function 
  ## containing glanced dataset
  ##
  ## Output:
  ## - A ratio scree plot with k (number of cluster centres) on the x-axis 
  ## and ratio on the y-axis
  ## 
  ## Example:
  ##   screeplot(cluster_info)
  ##
  clusterings <-
    cluster_info %>%
    unnest(cols = c(glanced))
  
  ggplot(clusterings, aes(k, tot.withinss)) +
    geom_line() +
    geom_point() 
    
  
  #ratio version
  nclust = length(clusterings$k)
  ratio = rep(NA, nclust-1)
  for (kk in 2:nclust) {
    ratio[kk-1] = clusterings$tot.withinss[kk]/clusterings$tot.withinss[kk-1]
  }
  plot_data <- data.frame(k = clusterings$k[2:nclust],ratio)
  ggplot(plot_data, aes(x=k, y = ratio)) + geom_line() + labs(caption="Figure 1. Scree plot")
  
}


#=========================================================================================
#=========================================================================================
colour_strips <- function(cluster_info){
  ## colour_strips(cluster_info) produces colour strips with the DMC colour
  ## closest to the cluster centre colour for every number of cluster centres in the k_list
  ##
  ## Input:
  ## - cluster_info: cluster information output from process_image function containing tidied dataset
  ##
  ## Output:
  ## - A number of colour strips for every element in k_list, containing DMC colours closest 
  ##   to the cluster centre colours
  ## 
  ## Example:
  ##   colour_strips(cluster_info)
  ##
  library(scales)
  par(mfrow = c(3, 3))
  for (i in 1:length(cluster_info$tidied)){
    show_col(cluster_info$tidied[[i]]$hex)
  }
}
#=========================================================================================
#=========================================================================================
change_resolution <- function(image_df, x_size)
{
  ## change_resolution(image_df, x_size) subsamples an image to produce
  ## a lower resolution image. Any non-coordinate columns in the data
  ## frame are summarized with their most common value in the larger
  ## grid cell.
  ##
  ## Input:
  ## - image_df: A data frame in wide format. The x-coordinate column MUST
  ##             be named 'x' and the y-coordinate column MUST be named 'y'.
  ##             Further columns have no naming restrictions.
  ## - x_size:   The number of cells in the x-direction. The number of cells
  ##             in the vertical direction will be computed to maintain the 
  ##             perspective. There is no guarantee that the exact number
  ##             of cells in the x-direction is x_size
  ##
  ## Output:
  ## - A data frame with the same column names as image_df, but with fewer 
  ##   entries that corresponds to the reduced resolution image.
  ##
  ## Example:
  ##   library(imager)
  ##   library(dplyr)
  ##   fpath <- system.file('extdata/Leonardo_Birds.jpg',package='imager') 
  ##   im <- load.image(fpath)
  ##   im_dat<- as.data.frame(im,wide = "c") %>% rename(R = c.1, G = c.2, B = c.3) %>%
  ##            select(x,y,R,G,B)
  ##   agg_image <- change_resolution(im_dat, 50)
  
  if(!require(sp)) {
    stop("The sp packages must be installed. Run install.packages(\"sp\") and then try again.")
  }
  if(!require(dplyr)) {
    stop("The dplyr packages must be installed. Run install.packages(\"dplyr\") and then try again.")
  }
  
  sp_dat <- image_df 
  gridded(sp_dat) = ~x+y
  
  persp = (gridparameters(sp_dat)$cells.dim[2]/gridparameters(sp_dat)$cells.dim[1])
  y_size = floor(x_size*persp)
  orig_x_size = gridparameters(sp_dat)$cells.dim[1]
  orig_y_size = gridparameters(sp_dat)$cells.dim[2]
  
  x_res = ceiling(orig_x_size/x_size)
  y_res = ceiling(orig_y_size/y_size)
  
  gt = GridTopology(c(0.5,0.5), c(x_res, y_res),
                    c(floor(orig_x_size/x_res), floor(orig_y_size/y_res)))
  SG = SpatialGrid(gt)
  agg = aggregate(sp_dat, SG, function(x) names(which.max(table(x)))[1] )
  agg@grid@cellsize <- c(1,1)
  df <- agg %>% as.data.frame %>% rename(x = s1, y = s2)  %>% select(colnames(image_df))
  
  return(df)
  
}

#=========================================================================================
#=========================================================================================
make_pattern <- function (cluster_info, k, x_size, black_white = FALSE, background_colour = NULL){
  ## make_pattern(cluster_info, k, x_size, black_white = FALSE, background_colour = NULL) 
  ## plots a cross-stitch pattern of the original image
  ## Input:
  ## - cluster_info: Cluster information output from process_image function
  ## - k: The chosen cluster size by user
  ## - x_size: The (approximate) total number of possible stitches in the horizontal 
  ## direction
  ## - black_white: (logical) Print the pattern in black and white (TRUE) or 
  ## colour(FALSE, default)
  ## - background_colour: The colour of the background, which should not be stitched 
  ## in the pattern. (Default is to not have a colour)
  ## Output:
  ## - A cross-stitch pattern that can be followed, complete with a legend that has 
  ## thread colour, and a guide grid
  ## 
  ## Example:
  ##   make_pattern(cluster_info, k=6, x_size= 50, black_white=FALSE, background_colour=NULL)
  ##
  library(imager)
  library(dplyr)
  library(cowplot)
  #augment the initial data with clusters 
  augmented <- augment(cluster_info[[1]]$kclust[[k-1]], cluster_info[[2]]) %>% rename(cluster = .cluster)
  augmented <- augmented %>% mutate(col = rgb(R,G,B))
  
  #change resolution 
  resized <- change_resolution(augmented, x_size)
  
  #plot
  
  ggplot(resized, aes(x=x, y = y,colour = cluster, shape=cluster)) +
    geom_point()+
    scale_colour_manual(values= cluster_info[[1]]$tidied[[k-1]]$hex) +
    scale_y_reverse() + theme_void()+ background_grid() + labs(caption="Cross Stitch Image")
}
