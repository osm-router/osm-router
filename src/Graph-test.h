/***************************************************************************
 *  Project:    osm-router
 *  File:       Graph-test.h
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
 *  osm-router.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Author:     Mark Padgham 
 *  E-Mail:     mark.padgham@email.com 
 *
 *  Description:    C++ implementation of OSM router using boost::graph.
 *                  Designed to work in a designated area, and so reads data
 *                  from a planet.osm file. Hard-coded at present to read data
 *                  for greater London and greater NYC, and to route between
 *                  points given in ./data/routing-points-(city).txt using the
 *                  profile given in profile.cfg
 *
 *  Limitations:
 *
 *  Dependencies:       libboost
 *
 *  Compiler Options:   -std=c++11 -lboost_program_options 
 ***************************************************************************/


#include "xml-parser.h"

#include <boost/unordered_set.hpp>
#include <boost/graph/connected_components.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>

#define PI 3.1415926535897932384626433832795

typedef boost::unordered_map <long long, int> umapInt; // for node int IDs
typedef boost::unordered_map <long long, int>::iterator umapInt_Itr;

typedef std::vector <Way>::iterator Ways_Itr;
typedef std::pair <std::string, float> ProfilePair;
typedef float Weight;

// See http://theboostcpplibraries.com/boost.unordered
std::size_t hash_value(const int &i)
{
    std::size_t seed = 0;
    boost::hash_combine(seed, i);
    return seed;
}

// edge and vertex bundles for the boost::graph
struct bundled_edge_type
{ 
    long long id;
    std::string name, type;
    // key-val pairs not implemented yet
    //std::vector <std::pair <std::string, std::string>> key_val;
    Weight weight; // weighted distance
    float dist;
};

struct bundled_vertex_type
{
    long long id;
    std::string name;
    float lat, lon;
    int component;
};



class Graph: public Xml
{
    // http://stackoverflow.com/questions/18791319/calculate-number-of-in-and-out-edges-in-a-boostgraph-vertex
    // say only "out_degree" is available for directedS, while bidiretionalS
    // enables also "in_degree"
    //using Graph_t = boost::adjacency_list< boost::vecS, boost::vecS, 
    //      boost::bidirectionalS, bundled_vertex_type, bundled_edge_type >;
    //using Graph_t = boost::adjacency_list< boost::vecS, boost::vecS, 
    //      boost::directedS, bundled_vertex_type, bundled_edge_type >;
    using Graph_t = boost::adjacency_list< boost::vecS, boost::vecS, 
          boost::undirectedS, bundled_vertex_type, bundled_edge_type >;

    using Vertex = boost::graph_traits<Graph_t>::vertex_descriptor;

    private:
        Graph_t gFull;

    protected:

    public:
        umapInt nodeIndx;
        umapPair_Itr node_itr;


    Graph (std::string file, float lonmin, float latmin, float lonmax, float latmax)
        : Xml (file, lonmin, latmin, lonmax, latmax)
    {
        // Construction fills vectors of "ways" and "nodes" from Xml class
        // These then need to be used to fill the boost::graph
        std::cout << "Parsed " << nodes.size () << " nodes and " << 
            ways.size () << " ways" << std::endl;

        nodeIndx.clear ();
        makeFullGraph ();
        dumpGraph ();
    }
    ~Graph ()
    {
        nodeIndx.clear ();
    }

    void dumpGraph ();
    void makeFullGraph ();
    float calcDist (std::vector <float> x, std::vector <float> y);
};



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                        FUNCTION::DUMPGRAPH                         **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void Graph::dumpGraph ()
{
    // Actually dumps the raw ways from which the graph is made. Useful for
    // visually inspecting connected components.
    long long ni;
    std::vector <float> lats, lons;
    umapPair_Itr umapitr;
    typedef std::vector <long long>::iterator ll_Itr;
    std::ofstream out_file;

    out_file.open ("junk.txt", std::ofstream::out);

    for (Ways_Itr wi = ways.begin(); wi != ways.end(); ++wi)
    {
        // Set up first origin node
        ni = (*wi).nodes.front ();
        assert ((umapitr = nodes.find (ni)) != nodes.end ());
        lats.resize (0);
        lons.resize (0);
        lons.push_back ((*umapitr).second.first);
        lats.push_back ((*umapitr).second.second);

        // Then iterate over the remaining nodes of that way
        for (ll_Itr it = std::next ((*wi).nodes.begin ());
                it != (*wi).nodes.end (); it++)
        {
            assert ((umapitr = nodes.find (*it)) != nodes.end ());
            lons.push_back ((*umapitr).second.first);
            lats.push_back ((*umapitr).second.second);

            assert (lons.size () == lats.size ()); // can't ever fail
            if (lons.size () == 2) // then add edge
            {
                out_file << lons [0] << "," << lats [0] << "," <<
                    lons [1] << "," << lats [1] << std::endl;
                lons.erase (lons.begin ());
                lats.erase (lats.begin ());
            }
        }
    }
    out_file.close ();

    // R code to examine output:
    /*
    plotgraph <- function ()
    {
        fname <- "./build/junk.txt"
        dat <- read.csv (fname, sep=",", header=FALSE)
        xlims <- range (c (dat [,1], dat [,3]))
        ylims <- range (c (dat [,2], dat [,4]))
        #xlims <- mean (xlims) + c (-0.05, 0.05) * diff (xlims)
        #ylims <- min (ylims) + 0.35 * diff (ylims) + c (-0.05, 0.05) * diff (ylims)
        plot.new ()
        par (mar=rep (0, 4))
        plot (NULL, NULL, xlim=xlims, ylim=ylims,
              xaxt="n", yaxt="n", xlab="", ylab="")
        junk <- apply (dat, 1, function (i)
                       lines (c (i [1], i [3]), c (i [2], i [4])))

        x <- c (dat [,1], dat [,3])
        y <- c (dat [,2], dat [,4])
        dat2 <- cbind (x, y)
        dat3 <- dat2 [which (!duplicated (dat2)),]
        points (dat2, pch=1)

        # zoom function (close window to stop)
        flag <- FALSE
        while (!flag)
        {
            lims <- par ("usr")
            loc <- locator (n=1)
            xlims <- loc$x + c (-0.25, 0.25) * diff (lims [1:2])
            xlims [1] <- max (c (xlims [1], min (x)))
            xlims [2] <- min (c (xlims [2], max (x)))
            ylims <- loc$y + c (-0.25, 0.25) * diff (lims [3:4])
            ylims [1] <- max (c (ylims [1], min (y)))
            ylims [2] <- min (c (ylims [2], max (y)))
            plot (NULL, NULL, xlim=xlims, ylim=ylims,
                  xaxt="n", yaxt="n", xlab="", ylab="")
            junk <- apply (dat, 1, function (i)
                           lines (c (i [1], i [3]), c (i [2], i [4])))
            points (dat2, pch=1)
        }
    }
    */

} // end function Graph::dumpGraph


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                      FUNCTION::MAKEFULLGRAPH                       **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void Graph::makeFullGraph ()
{
    int nodeCount = 0, tempi [2];
    float lon1, lat1, lon2, lat2;
    std::vector <float> lats, lons;
    long long ni;
    typedef std::vector <long long>::iterator ll_Itr;
    bundled_vertex_type oneVert;
    bundled_edge_type oneEdge;
    umapPair_Itr umapitr;

    // The first of these is to check whether edges exist; the second to count
    // them
    typedef boost::graph_traits <Graph_t>::edge_descriptor edge_t;
    boost::graph_traits <Graph_t>::out_edge_iterator ei, ei_end;

    // Vertices are added as new ones appear in ways. boost::graph numbers all
    // vertices with sequential integers indexed here in nodeIndx. The tempi
    // values hold the indices for add the edges.
    for (Ways_Itr wi = ways.begin(); wi != ways.end(); ++wi)
    {
        oneEdge.id = (*wi).id;
        oneEdge.name = (*wi).name;
        oneEdge.type = (*wi).type;

        // Set up first origin node
        ni = (*wi).nodes.front ();
        assert ((umapitr = nodes.find (ni)) != nodes.end ());
        // Add to nodes if not already there
        if (nodeIndx.find (ni) == nodeIndx.end())
        {
            oneVert.id = (*umapitr).first;
            oneVert.lon = (*umapitr).second.first;
            oneVert.lat = (*umapitr).second.second;
            boost::add_vertex (oneVert, gFull);
            tempi [0] = nodeCount;
            nodeIndx [oneVert.id] = nodeCount++;
        } else
        {
            tempi [0] = (*nodeIndx.find (ni)).second; // int index of ni
        }

        lats.resize (0);
        lons.resize (0);
        lons.push_back ((*umapitr).second.first);
        lats.push_back ((*umapitr).second.second);

        // Then iterate over the remaining nodes of that way
        for (ll_Itr it = std::next ((*wi).nodes.begin ());
                it != (*wi).nodes.end (); it++)
        {
            assert ((umapitr = nodes.find (*it)) != nodes.end ());
            // Add to nodes if not already there
            if (nodeIndx.find (*it) == nodeIndx.end ())
            {
                oneVert.id = (*umapitr).first;
                oneVert.lon = (*umapitr).second.first;
                oneVert.lat = (*umapitr).second.second;
                boost::add_vertex (oneVert, gFull);
                tempi [1] = nodeCount;
                nodeIndx [oneVert.id] = nodeCount++;
            } else
            {
                tempi [1] = (*nodeIndx.find (*it)).second;
            }
            
            lons.push_back ((*umapitr).second.first);
            lats.push_back ((*umapitr).second.second);

            // Check whether edge exists:
            std::pair <edge_t, bool> p = boost::edge (tempi [0], tempi [1], gFull);
            if (p.second) 
            {
                edge_t e = *boost::out_edges (tempi [0], gFull).first;
                if (gFull [e].type != (*wi).type)
                    std::cout << "Repeated edge of different type: " <<
                std::cout << "edge exists: " << p.first << ": (N=" <<
                    gFull [e].name << "; ID=" << gFull [e].id << "; type=" <<
                    gFull [e].type << " -> " << (*wi).type << std::endl;
            }
            // or count them:
            /*
            boost::tie (ei, ei_end) = out_edges (tempi [0], gFull);
            int parallel_count = 0;
            for ( ; ei != ei_end; ++ei)
                if (boost::target (*ei, gFull) == tempi [1])
                    parallel_count++;
            if (parallel_count > 0)
                std::cout << "---" << parallel_count << "---" << std::endl;
            */

            assert (lons.size () == lats.size ()); // can't ever fail
            oneEdge.dist = calcDist (lons, lats);
            boost::add_edge (tempi [0], tempi [1], oneEdge, gFull);
            tempi [0] = tempi [1];
            lons.erase (lons.begin ());
            lats.erase (lats.begin ());
        }
    }

    std::vector <int> compvec (num_vertices (gFull));
    tempi [0] = boost::connected_components (gFull, &compvec [0]);
    std::cout << "Graph has " << num_vertices (gFull) << " vertices and " <<
        tempi [0] << " connected components" << std::endl;
} // end function Graph::makeFullGraph



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         FUNCTION::CALCDIST                         **
 **                                                                    **
 ************************************************************************
 ************************************************************************/


float Graph::calcDist (std::vector <float> x, std::vector <float> y)
{
    float d, dsum = 0.0, x0, y0, x1 = x[0], y1 = y[0], xd, yd;
    assert (x.size () == y.size ());

    for (int i=1; i<x.size (); i++)
    {
        x0 = x1;
        y0 = y1;
        x1 = x[i];
        y1 = y[i];
        
        xd = (x1 - x0) * PI / 180.0;
        yd = (y1 - y0) * PI / 180.0;

        d = sin (yd / 2.0) * sin (yd / 2.0) + cos (y1 * PI / 180.0) *
            cos (y0 * PI / 180.0) * sin (xd / 2.0) * sin (xd / 2.0);
        d = 2.0 * atan2 (sqrt (d), sqrt (1.0 - d));

        dsum += d * 6371.0;
    }

    return dsum; // in kilometres!
} // end function Graph::calcDist

