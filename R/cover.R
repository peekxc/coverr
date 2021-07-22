#' interval_cover
#' @param x point cloud to construct a cover of.
#' @param num_intervals number of intervals in the cover, per dimension. 
#' @param overlap proportion of overlap between intervals.
#' @description This function constructs an interval over on the unit hypercube [0,1]^d.
#' @useDynLib coverr, .registration = TRUE
#' @import methods Rcpp
#' @export
interval_cover <- function(x, num_intervals, overlap, bounded = FALSE){
	if (is.vector(x) && is.null(dim(x))){ x <- as.matrix(x) }
	stopifnot(!is.null(dim(x)))
	d <- ncol(x)
	stopifnot(length(num_intervals) %in% c(1L, d), length(overlap) %in% c(1L, d))
	if (length(num_intervals) != d){ num_intervals <- rep(num_intervals, d) }
	if (length(overlap) != d){ overlap <- rep(overlap, d) }
	R <- apply(x, 2, range)
	
	## If bounded, used fixed interval cover 
	if (!bounded){
		cover <- local({
			base_width <- apply(R, 2, diff)/num_intervals
			interval_width <- base_width + (base_width*overlap)/(1.0 - overlap)
			eps <- (interval_width/2.0) + sqrt(.Machine$double.eps) ## ensures each point is in the cover
		  cart_prod <- arrayInd(seq(prod(num_intervals)), .dim = num_intervals)
		  set_bounds <- t(apply(cart_prod, 1, function(idx){
		    centroid <- R[1L,] + ((as.integer(idx)-1L)*base_width) + base_width/2.0
		    c(centroid - eps, centroid + eps)
		  }))
		  constructIsoAlignedLevelSets(x, set_bounds)
		})
	} else {
		cover <- local({
			interval_width <- apply(R, 2, diff)/(num_intervals - overlap*(num_intervals - 1L))
	  	base_width <- interval_width * (1 - overlap)
	  	cart_prod <- arrayInd(seq(prod(num_intervals)), .dim = num_intervals)
	    set_bounds <- t(apply(cart_prod, 1, function(idx){
	      set_lb <- R[1,] + (as.integer(idx)-1L) * base_width
	      c(set_lb, set_lb + interval_width)
	    }))
    	constructIsoAlignedLevelSets(x, set_bounds)
		})
	}
	subset_sizes <- sapply(cover, length)
	cover_sm <- Matrix::sparseMatrix(
		i = unlist(cover), 
		j = rep(seq_along(cover), times = subset_sizes), 
		x = rep(TRUE, sum(subset_sizes))
	)
	return(cover_sm)
}

#' Ball cover 
#' @description Constructs a ball cover over a point set at a given radius. 
#' @param x point cloud, as a (n x d) matrix representing `n` points in `d` dimensions. 
#' @param centers ball centers, as a (j x d) matrix representing `j` points in `d` dimensions.
#' @param radii either the radius of each ball, or a single radius for all the balls
#' @details This function determines the set of points that lie within the union of closed balls 
#' of some supplied `radius`.
#' @export
ball_cover <- function(x, centers, radii, compress=FALSE){
	stopifnot(is.matrix(x), is.matrix(centers), ncol(x) == ncol(centers))
	radii <- rep(radii, length.out = nrow(centers))
	cover <- Matrix::which(t(t(proxy::dist(x, centers)) <= radii), arr.ind = TRUE)
	if (compress){
		cover <- as.data.frame(cover[,c(2,1),drop=FALSE])
		colnames(cover) <- c("set", "element")
		cover <- RcppGreedySetCover::greedySetCover(cover, data.table = FALSE)
		cover <- as.matrix(cover[,c(2,1)])
	}
	cover_sm <- Matrix::sparseMatrix(i = cover[,1L], j = cover[,2L], x = TRUE, dims = c(nrow(x), length(unique(cover[,2L]))))
	
	## Check if the whole data set is covered
	cover_check <- vector(mode = "logical", length = nrow(x))
	cover_check[cover_sm@i+1L] <- TRUE
	if (any(!cover_check)){
		warning(sprintf("The union of balls of radius %g only covers %g\\% of the data set!", radius, sum(cover_check)/length(cover_check)*100))
	}
	return(cover_sm)
}

#' Landmark cover
#' @description Creates a cover based on landmarks
#' @export
landmark_cover <- function(x, n_sets=25L, radius=NULL, ...){
	if (!missing(radius) && is.numeric(radius)){
		point_indices <- landmark::landmarks_maxmin(x, radius = radius)
	} else if (missing(radius)){
		stopifnot(is.numeric(n_sets), n_sets == floor(n_sets))
		point_indices <- landmark::landmarks_maxmin(x, num = n_sets)
		radius <- min(proxy::dist(x = x[point_indices[-n_sets],,drop=FALSE], y = x[point_indices[n_sets],,drop=FALSE]))
	}
	cover_sm <- ball_cover(x, centers = x[point_indices,,drop=FALSE], radii = radius, ...)
	return(cover_sm)
}


## TODO: figure out how to scale the balls to guarantee the result is a cover  

#' K-Nearest Neighbor cover
#' @description This function determines the set of points that lie within the union of the top `n_sets`
#' densest balls according to a `k`-NN density estimate. The radii of the these balls is proportional
#' to distance to their k-th nearest neighbors.
#' @export
knn_cover <- function(x, k, n_sets=25L){
	knn_dist <- apply(RANN::nn2(x, x, k = k)$nn.dist, 1, max)
	centers <- order(knn_dist)[seq(n_sets)]
	return(ball_cover(x, centers = x[centers,,drop=FALSE], radii = knn_dist[centers]))
}

