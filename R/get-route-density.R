#' get_routes
#' 
#' Calculates the probability density of routes between two points.
#' 
#' @param bbox the bounding box defining the road network. If not available, it
#' will be downloaded using the osmdatar package.
#' @param start the starting point of the route
#' @param end the ending point of the route
#' @export

library (osmdatar);
library (methods); # needed when using Rscript to run
library (Rcpp);

get_routes <- function (bbox, start, end)
{
    message ("Getting roads...")
    roads <- osmdatar::get_lines (bbox, key='highway');
    
    return (roads)
}

### For testing
start <- function ()
{
    bbox <- matrix (c (-0.11, 51.51, -0.10, 51.52), nrow=2, ncol=2)
    roads <- get_routes (bbox, NULL, NULL)
    summary (roads)
}

start ();
