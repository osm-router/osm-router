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

get_routes <- function (bbox, start, end)
{
    roads <- osmdatar::get_lines (bbox, key='highway');
}
