/***************************************************************************
 *  Project:    osm-router
 *  File:       Utils.h
 *  Language:   C++
 *
 *  osm-router is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation, either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  osm-router is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 *  more details.
 *
 *  You should have received a copy of the GNU General Public License along with
 *  NeutralClusters.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright   Mark Padgham December 2015
 *  Author:     Mark Padgham
 *  E-Mail:     mark.padgham@email.com
 *
 *  Description:    Routing engine for OSM based on boost::graph
 *
 *  Limitations:
 *
 *  Dependencies:       libboost
 *
 *  Compiler Options:   -std=c++11 -lboost_program_options 
 ***************************************************************************/

#include <stdlib.h> // has abs function
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <vector>
#include <string>
#include <iomanip> // for setfill
#include <sys/ioctl.h> // for console width: Linux only!
#include <ctype.h>
#include <fstream>
#include <assert.h>

#include <boost/config.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <boost/program_options.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>

#ifndef UTILS_H
#define UTILS_H

#define PI 3.1415926535897932384626433832795

typedef boost::numeric::ublas::vector <int> ivec;
typedef boost::numeric::ublas::matrix <int> imat;
typedef boost::numeric::ublas::vector <double> dvec;
typedef boost::numeric::ublas::matrix <double> dmat;
typedef boost::numeric::ublas::vector <bool> bvec;
typedef boost::numeric::ublas::matrix <bool> bmat;
typedef boost::numeric::ublas::zero_matrix <double> zmat_d;
typedef boost::numeric::ublas::zero_matrix <int> zmat_i;

const double DOUBLE_MAX = std::numeric_limits<double>::max (),
    DOUBLE_MIN = -DOUBLE_MAX,
    FLOAT_MAX = std::numeric_limits <float>::max ();

struct myTime{
    int hh, mm;
    float ss;	};

struct DistStruct{
    double dx, dy, d;	};


struct RegrResults {
    double r2, cov, slope, intercept, SS, tval;      };

std::string standardise (std::string sin);
std::string substituteNames (bool tube, std::string str);
double calc_angle (double x, double y);
DistStruct getdists (double xa, double ya, double xb, double yb);
DistStruct convert_distance (double dist, double midx, double midy);
void timeout (double tseconds);
RegrResults regression (std::vector <double> x, std::vector <double> y);
double calcMI (std::vector <double> x, std::vector <double> y);
void progLine (double progress);

#endif
