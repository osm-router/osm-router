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
 *  Author:     Mark Padgham / Andreas Petutschnig
 *  E-Mail:     mark.padgham@email.com / andras@petutschnig.de
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
    bool oneway;
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
    // stackoverflow: "calculate-number-of-in-and-out-edges-in-a-boostgraph-vertex"
    // say only "out_degree" is available for directedS, while bidiretionalS
    // also enables "in_degree"
    //using Graph_t = boost::adjacency_list< boost::vecS, boost::vecS, 
    //      boost::bidirectionalS, bundled_vertex_type, bundled_edge_type >;
    using Graph_t = boost::adjacency_list< boost::vecS, boost::vecS, 
          boost::directedS, bundled_vertex_type, bundled_edge_type >;
    //using Graph_t = boost::adjacency_list< boost::vecS, boost::vecS, 
    //      boost::undirectedS, bundled_vertex_type, bundled_edge_type >;

    using Vertex = boost::graph_traits<Graph_t>::vertex_descriptor;

    private:
        bool _compact_graph, obey_oneway;
        std::string _profile_file;

    protected:
        std::vector <ProfilePair> profile;
        Graph_t gFull, gCompact;

    public:
        umapInt nodeIndxFull, nodeIndxCompact;
        umapPair_Itr node_itr;


    Graph (bool compact, std::string xml_file, std::string profile_file,
            float lonmin, float latmin, float lonmax, float latmax)
        : _compact_graph (compact), _profile_file (profile_file),
            Xml (xml_file, lonmin, latmin, lonmax, latmax)
    {
        // Construction fills vectors of "ways" and "nodes" from Xml class
        // These then need to be used to fill the boost::graph
        // #nodes = nodes.size ()
        // #ways = ways.size ()

        boost::filesystem::path p (_profile_file);
        if (boost::filesystem::exists (p))
            setProfile (_profile_file.c_str (), &profile);
        // setProfile also sets obey_oneway

        nodeIndxFull.clear ();
        nodeIndxCompact.clear ();
        if (_compact_graph)
            makeCompactGraph ();
        else
            makeFullGraph ();
        //dumpWays ();
    }
    ~Graph ()
    {
        nodeIndxFull.clear ();
        nodeIndxCompact.clear ();
        gFull.clear ();
        gCompact.clear ();
    }

    std::string getProfileFile () { return _profile_file;   }
    bool getIsCompact () { return _compact_graph;   }

    void dumpWays ();
    void dumpGraph (Graph_t* g);
    void makeFullGraph ();
    void makeCompactGraph ();
    float calcDist (std::vector <float> x, std::vector <float> y);
    void add_edge (int v1, int v2, Graph_t* g,
            bundled_edge_type* oneEdge, std::string new_highway_type);
    int getConnected (Graph_t *g);
    void setProfile (const std::string& profileName, 
            std::vector <ProfilePair>* profile);
    float getProfileWeight (std::vector <ProfilePair>* profile,
            std::string highway_type);
};



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                        FUNCTION::DUMPWAYS                          **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void Graph::dumpWays ()
{
    // All ways parsed from XML. Useful for visually inspecting connected
    // components and comparing with results of dumpGraph.
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

} // end function Graph::dumpWays


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                        FUNCTION::DUMPGRAPH                         **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void Graph::dumpGraph (Graph_t* g)
{
    // TODO: Write this!
}


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
    float wt, lon1, lat1, lon2, lat2;
    std::vector <float> lats, lons;
    long long ni;
    typedef std::vector <long long>::iterator ll_Itr;
    bundled_vertex_type oneVert;
    bundled_edge_type oneEdge;
    umapPair_Itr umapitr;

    // For edge counting
    //boost::graph_traits <Graph_t>::out_edge_iterator ei, ei_end;

    // Vertices are added as new ones appear in ways. boost::graph numbers all
    // vertices with sequential integers indexed here in nodeIndxFull. The tempi
    // values hold the indices for add the edges.
    for (Ways_Itr wi = ways.begin(); wi != ways.end(); ++wi)
    {
        wt = getProfileWeight (&profile, (*wi).type);
        // only proceed if highway type is in profile
        if (wt > 0.0)
        {
            oneEdge.id = (*wi).id;
            oneEdge.name = (*wi).name;
            oneEdge.type = (*wi).type;

            // Set up first origin node
            ni = (*wi).nodes.front ();
            assert ((umapitr = nodes.find (ni)) != nodes.end ());
            // Add to nodes if not already there
            if (nodeIndxFull.find (ni) == nodeIndxFull.end())
            {
                oneVert.id = (*umapitr).first;
                oneVert.lon = (*umapitr).second.first;
                oneVert.lat = (*umapitr).second.second;
                boost::add_vertex (oneVert, gFull);
                tempi [0] = nodeCount;
                nodeIndxFull [oneVert.id] = nodeCount++;
            } else
            {
                tempi [0] = (*nodeIndxFull.find (ni)).second; // int index of ni
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
                if (nodeIndxFull.find (*it) == nodeIndxFull.end ())
                {
                    oneVert.id = (*umapitr).first;
                    oneVert.lon = (*umapitr).second.first;
                    oneVert.lat = (*umapitr).second.second;
                    boost::add_vertex (oneVert, gFull);
                    tempi [1] = nodeCount;
                    nodeIndxFull [oneVert.id] = nodeCount++;
                } else
                {
                    tempi [1] = (*nodeIndxFull.find (*it)).second;
                }
                
                lons.push_back ((*umapitr).second.first);
                lats.push_back ((*umapitr).second.second);
                oneEdge.dist = calcDist (lons, lats);
                oneEdge.weight = oneEdge.dist / wt;

                // count edges (can also count in_edges for bidirectionalS)
                /*
                boost::tie (ei, ei_end) = out_edges (tempi [0], gFull);
                int parallel_count = 0;
                for ( ; ei != ei_end; ++ei)
                    if (boost::target (*ei, gFull) == tempi [1])
                        parallel_count++;
                */

                add_edge (tempi [0], tempi [1], &gFull, &oneEdge, (*wi).type);
                if (!(*wi).oneway || !obey_oneway)
                    add_edge (tempi [1], tempi [0], &gFull, &oneEdge, (*wi).type);

                assert (lons.size () == lats.size ()); // can't ever fail
                tempi [0] = tempi [1];
                lons.erase (lons.begin ());
                lats.erase (lats.begin ());
            } // end for ll_Itr over nodes
        } // end if wt > 0.0
    } // end for Ways_itr over ways

    std::cout << "Full graph has " << num_vertices (gFull) << " vertices and " <<
        getConnected (&gFull) << " connected components" << std::endl;
} // end function Graph::makeFullGraph


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         FUNCTION::ADD_EDGE                         **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void Graph::add_edge (int v1, int v2, Graph_t* g,
        bundled_edge_type* oneEdge, std::string new_highway_type)
{
    typedef boost::graph_traits <Graph_t>::edge_descriptor edge_t;
    std::pair <edge_t, bool> edge_p = boost::edge (v1, v2, *g);
    if (edge_p.second) // edge already exists
    {
        edge_t e = *boost::out_edges (v1, *g).first;
        if ((*g)[e].type != new_highway_type)
            std::cout << "Repeated edge of different type: " <<
        std::cout << "edge exists: " << edge_p.first << ": (N=" <<
            (*g)[e].name << "; ID=" << (*g)[e].id << "; type=" <<
            (*g)[e].type << " -> " << new_highway_type << std::endl;
    } else
    {
        boost::add_edge (v1, v2, *oneEdge, *g);
    }
} // end function Graph::add_edge


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                     FUNCTION::MAKECOMPACTGRAPH                     **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void Graph::makeCompactGraph ()
{
    int nodeCount = 0, tempi [2];
    float wt, lon1, lat1, lon2, lat2;
    std::vector <float> lats, lons;
    long long ni;
    typedef std::vector <long long>::iterator ll_Itr;
    bundled_vertex_type oneVert;
    bundled_edge_type oneEdge;
    umapPair_Itr umapitr;

    // First count all occurrences of each node so compact graph includes only
    // those nodes with counts > 1. To start with, however, all start and end
    // nodes of each way are given counts of 2 to ensure they're included.
    umapInt nodeCounts;
    for (Ways_Itr wi = ways.begin(); wi != ways.end(); ++wi)
    {
        wt = getProfileWeight (&profile, (*wi).type);
        // only proceed if highway type is in profile
        if (wt > 0.0)
        {
            nodeCounts [(*wi).nodes.front ()] = 2;
            nodeCounts [(*wi).nodes.back ()] = 2;
            for (ll_Itr it = std::next ((*wi).nodes.begin ());
                    it != (*wi).nodes.end (); it++)
            {
                assert ((umapitr = nodes.find (*it)) != nodes.end ());
                if (nodeCounts.find (*it) == nodeCounts.end ())
                    nodeCounts [*it] = 1;
                else
                    nodeCounts [*it]++;
            }
        }
    }

    // Vertices are added as new ones appear in ways. boost::graph numbers all
    // vertices with sequential integers indexed here in nodeIndxCompact. The
    // tempi values hold the indices for add the edges.
    for (Ways_Itr wi = ways.begin(); wi != ways.end(); ++wi)
    {
        wt = getProfileWeight (&profile, (*wi).type);
        // only proceed if highway type is in profile
        if (wt > 0.0)
        {
            oneEdge.id = (*wi).id;
            oneEdge.name = (*wi).name;
            oneEdge.type = (*wi).type;

            // Set up first origin node
            ni = (*wi).nodes.front ();
            assert ((umapitr = nodes.find (ni)) != nodes.end ());
            // Add to nodes if not already there (and first nodes are always added)
            if (nodeIndxCompact.find (ni) == nodeIndxCompact.end())
            {
                oneVert.id = (*umapitr).first;
                oneVert.lon = (*umapitr).second.first;
                oneVert.lat = (*umapitr).second.second;
                boost::add_vertex (oneVert, gCompact);
                tempi [0] = nodeCount;
                nodeIndxCompact [oneVert.id] = nodeCount++;
            } else
            {
                tempi [0] = (*nodeIndxCompact.find (ni)).second; // int index of ni
            }

            lats.resize (0);
            lons.resize (0);
            lons.push_back ((*umapitr).second.first);
            lats.push_back ((*umapitr).second.second);

            // Then iterate over the remaining nodes of that way
            oneEdge.dist = 0.0;
            for (ll_Itr it = std::next ((*wi).nodes.begin ());
                    it != (*wi).nodes.end (); it++)
            {
                assert ((umapitr = nodes.find (*it)) != nodes.end ());
                lons.push_back ((*umapitr).second.first);
                lats.push_back ((*umapitr).second.second);
                oneEdge.dist += calcDist (lons, lats);
                oneEdge.weight = oneEdge.dist / wt;

                // Nodes are only potentially added here if nodeCounts > 1,
                // otherwise edge distance is simply accumulated.
                if (nodeCounts [*it] > 1)
                {
                    // Add to nodes if not already there *AND* nodeCounts > 1
                    if (nodeIndxCompact.find (*it) == nodeIndxCompact.end ())
                    {
                        oneVert.id = (*umapitr).first;
                        oneVert.lon = (*umapitr).second.first;
                        oneVert.lat = (*umapitr).second.second;
                        boost::add_vertex (oneVert, gCompact);
                        tempi [1] = nodeCount;
                        nodeIndxCompact [oneVert.id] = nodeCount++;
                    } else
                    {
                        tempi [1] = (*nodeIndxCompact.find (*it)).second;
                    }
                    add_edge (tempi [0], tempi [1], &gCompact, &oneEdge, (*wi).type);
                    if (!(*wi).oneway || !obey_oneway)
                        add_edge (tempi [1], tempi [0], &gCompact, &oneEdge, (*wi).type);

                    tempi [0] = tempi [1];
                    oneEdge.dist = 0.0;
                }
                assert (lons.size () == lats.size ()); // can't ever fail
                lons.erase (lons.begin ());
                lats.erase (lats.begin ());
            } // end for ll_Itr over nodes
        } // end if wt > 0.0
    } // end for Ways_itr over ways
    nodeCounts.clear ();

    std::cout << "Compact graph has " << num_vertices (gCompact) << 
        " vertices and " << getConnected (&gCompact) << 
        " connected components" << std::endl;
} // end function Graph::makeCompactGraph


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



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                       FUNCTION::GETCONNECTED                       **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Graph::getConnected (Graph_t *g)
{
    std::vector <int> compvec (num_vertices (*g));
    int num = boost::connected_components (*g, &compvec[0]);

    // Then store component info in vertices
    /*
    typedef boost::graph_traits <Graph_t>::vertex_iterator viter;
    std::pair <viter, viter> vp;
    for (vp = vertices (*g); vp.first != vp.second; ++vp.first)
        vertex_component [*vp.first] = compvec [*vp.first];
     */

    // Alternative:
    boost::property_map< Graph_t, int bundled_vertex_type::* >::type 
        vertex_component = boost::get (&bundled_vertex_type::component, (*g));
    auto vs = boost::vertices (*g);
    for (auto vit = vs.first; vit != vs.second; ++vit)
        vertex_component [*vit] = compvec [*vit];

    // Optional filtering of component = 0:
    /*
    in_component_0 <VertMap> filter (boost::get 
        (&bundled_vertex_type::component, (*g)));
    boost::filtered_graph <Graph, in_component_0 <VertMap> > fg ((*g), filter);
    // Next lines reveal filtering does not actually reduce the graph
    auto vsfg = boost::vertices (fg);
    for (auto vit = vsfg.first; vit != vsfg.second; ++vit)
        std::cout << *vit << std::endl;
     */

    return num;
} // end function getConnected


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                        FUNCTION::SETPROFILE                        **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void Graph::setProfile (const std::string& profileName, 
        std::vector <ProfilePair>* profile)
{
    int ipos;
    float value;
    std::string line, field;
    std::ifstream in_file;

    // Profile can also include an "obey" key. If this is 0, then oneway
    // highways will *not* be obeyed, and will thus be added to the graph.
    obey_oneway = true;

    profile->resize (0);

    in_file.open (profileName.c_str (), std::ifstream::in);
    assert (!in_file.fail ());

    while (!in_file.eof ())
    {
        getline (in_file, line, '\n');
        if (!in_file.eof ())
        {
            ipos = line.find (',');
            field = line.substr (0, ipos);
            line = line.substr (ipos + 1, line.length () - ipos - 1);
            value = atof (line.c_str ());

            if (field.substr (0, 4) == "obey" && value == 0.0)
                obey_oneway = false;

            profile->push_back (std::make_pair (field, value));
        }
    }
    in_file.close ();
};


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                     FUNCTION::GETPROFILEWEIGHT                     **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

float Graph::getProfileWeight (std::vector <ProfilePair>* profile,
        std::string highway_type)
{
    float wt = -1.0; // negative return if highway_type not in profile

    for (std::vector <ProfilePair>::iterator itr = (*profile).begin();
            itr != (*profile).end(); itr++)
    {
        if ((*itr).first == highway_type)
        {
            wt = (*itr).second;
            break;
        }
    }

    return wt;
}
