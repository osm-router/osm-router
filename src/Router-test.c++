/***************************************************************************
 *  Project:    osm-router
 *  File:       Router.c++
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

#include "Router-test.h"

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                                MAIN                                **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int main(int argc, char *argv[]) {
    bool compact;
    float xfrom, yfrom, xto, yto;
    Bbox bbox;
    std::string xml_file, profile_file;

    try {
        boost::program_options::options_description generic("Generic options");
        generic.add_options()
            ("version,v", "print version std::string")
            ("help", "produce help message")    
            ;

        boost::program_options::options_description config("Configuration");
        config.add_options()
            ("graph compact, g", boost::program_options::value <bool>
                (&compact)->default_value (true),
                "use compact graph")
            ("xml_file,f", boost::program_options::value <std::string> 
                (&xml_file)->default_value ("xmldat.xml"), 
                "xml_file name (.xml will be appended)")
            ("profile_file,p", boost::program_options::value <std::string> 
                (&profile_file)->default_value ("../profile.cfg"), 
                "profile file name")
            ("xfrom,a", boost::program_options::value <float> 
                (&xfrom)->default_value (-0.12), "xfrom")
            ("yfrom,b", boost::program_options::value <float> 
                (&yfrom)->default_value (51.515), "yfrom")
            ("xto,c", boost::program_options::value <float> 
                (&xto)->default_value (-0.115), "xto")
            ("yto,d", boost::program_options::value <float> 
                (&yto)->default_value (51.52), "yto")
            ;

        boost::program_options::options_description cmdline_options;
        cmdline_options.add(generic).add(config);

        boost::program_options::options_description visible("Allowed options");
        visible.add(generic).add(config);

        boost::program_options::variables_map vm;
        store(boost::program_options::command_line_parser(argc, argv).
                options(cmdline_options).run(), vm);

        notify(vm);

        if (vm.count("help")) {
            std::cout << visible << std::endl;
            return 0;
        }

        if (vm.count("version")) {
            std::cout << "osm-router, version 1.0" << std::endl;
            return 0;
        }

    }
    catch(std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }    

    bbox = get_bbox (xfrom, yfrom, xto, yto);

    Router router (xml_file, profile_file, 
            bbox.lonmin, bbox.latmin, bbox.lonmax, bbox.latmax, compact); 
};


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                               GETBBOX                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

Bbox get_bbox (float xfrom, float yfrom, float xto, float yto, float expand)
{
    // expand is in both directions, so actually twice that value
    float xmin, ymin, xmax, ymax, xrange, yrange;
    Bbox bbox;

    if (xfrom < xto)
    {
        xmin = xfrom;
        xmax = xto;
    } else
    {
        xmin = xto;
        xmax = xfrom;
    }
    if (yfrom < yto)
    {
        ymin = yfrom;
        ymax = yto;
    } else
    {
        ymin = yto;
        ymax = yfrom;
    }

    xrange = xmax - xmin;
    yrange = ymax - ymin;

    bbox.lonmin = xmin - xrange * expand;
    bbox.lonmax = xmax + xrange * expand;
    bbox.latmin = ymin - yrange * expand;
    bbox.latmax = ymax + yrange * expand;

    return bbox;
}

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         READROUTINGPOINTS                          **
 **                                                                    **
 ************************************************************************
 ************************************************************************/


int Router::readRoutingPoints ()
{
    int station_id, ipos = 0;
    std::string linetxt, txt, fname, stationNodeFileName;

    /*
     * If the routing points are matched to the OSM nodes the first time, 
     * the matching results are written into the file stationNodeFileName.
     * If the file already exists, the matching results are taken from there
     * instead of newly calculating them.
     */

    if (access (stationNodeFileName.c_str (), F_OK) != -1)
    {
        std::ifstream stationNodeFile;
        stationNodeFile.open (stationNodeFileName, std::ifstream::in);
        assert (!stationNodeFile.fail ());

        stationNodeFile.clear ();
        stationNodeFile.seekg (0);
        getline (stationNodeFile, linetxt, '\n'); //header
        while (getline (stationNodeFile, linetxt, '\n'))
            ipos++;

        stationNodeFile.clear ();
        stationNodeFile.seekg (0);
        getline (stationNodeFile, linetxt, '\n'); //header

        std::cout << "Matching " << ipos << " routing points to nearest OSM nodes, using previously calculated nodes from file " << stationNodeFileName << "...";

        RoutingPoint rpoint;

        while (getline (stationNodeFile, linetxt, '\n'))
        {
            ipos = linetxt.find (",");
            station_id = atof (linetxt.substr (0, ipos).c_str());
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            ipos = linetxt.find (",");
            rpoint.node = atof (linetxt.substr (0, ipos).c_str());
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            ipos = linetxt.find (",");
            rpoint.nodeIndex = atof (linetxt.substr (0, ipos).c_str());
            RoutingPointsList.push_back (rpoint);

        }
        stationNodeFile.close ();

        std::cout << " done." << std::endl;
    }
    else
    {
        std::ifstream in_file;

        RoutingPoint rpoint;

        in_file.open (fname.c_str (), std::ifstream::in);
        assert (!in_file.fail ());

        in_file.clear ();
        in_file.seekg (0); 
        getline (in_file, linetxt, '\n'); // header
        while (getline (in_file, linetxt, '\n'))
            ipos++;

        in_file.clear ();
        in_file.seekg (0); 
        getline (in_file, linetxt, '\n'); 

        RoutingPointsList.resize (0);
        std::cout << "Matching " << ipos << " routing points to nearest OSM nodes ...";
        std::cout.flush ();


        /*
         * the main task of this routine is to assign routing points to their
         * nearest highway nodes. These nodes must be within the largest connected
         * component of the OSM graph, which boost *seems* always seems to number
         * with 0.
         *
         * Note that nearestNode is applied to gFull, not gCompact, because the
         * nearest nodes are included as terminalNodes which are used to build
         * gComact.  remapRoutingPoints then realigns routing point indexes into
         * gCompact.
         */

        std::ofstream stationNodeFile;
        stationNodeFile.open (stationNodeFileName);
        stationNodeFile << "id, node, nodeIndex" << '\n';

        while (getline (in_file, linetxt, '\n'))
        {
            ipos = linetxt.find (",");
            station_id = atof (linetxt.substr (0, ipos).c_str());
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            ipos = linetxt.find (",");
            rpoint.lat = atof (linetxt.substr (0, ipos).c_str());
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            ipos = linetxt.find (",");
            rpoint.lon = atof (linetxt.substr (0, ipos).c_str());
            rpoint.node = nearestNode (rpoint.lon, rpoint.lat);
            //assert (nodeNames.find (rpoint.node) != nodeNames.end());
            //rpoint.nodeIndex = nodeNames.find (rpoint.node)->second;
            RoutingPointsList.push_back (rpoint);

            stationNodeFile << station_id << ", " << rpoint.node << ", " << rpoint.nodeIndex << '\n';

        } // end while getline
        stationNodeFile.close ();
        in_file.close ();
        std::cout << " done." << std::endl;
    }

    return RoutingPointsList.size ();
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                            NEARESTNODE                             **
 **                                                                    **
 ************************************************************************
 ************************************************************************/


long long Router::nearestNode (float lon0, float lat0)
{
    int id;
    long long node = -INT_MAX;
    float d, lon, lat, xd, yd, dmin = FLOAT_MAX;

    /*
     * This is the only place in the whole program that requires explicit
     * looping over vertices, and takes longer than all other bits. Calculating
     * distances to nearest vertices as sums of the 2 abs distances in lats &
     * lons saves only about 2.5s.
     *
     * nodeNames lists the highway vertices, which must be matched to the full
     * list of allNodes to find corresponding lats&lons.
     */

    boost::property_map< Graph_t, long long bundled_vertex_type::* >::type 
        vertex_id = boost::get(&bundled_vertex_type::id, gFull);
    boost::property_map< Graph_t, float bundled_vertex_type::* >::type 
        vertex_lat = boost::get(&bundled_vertex_type::lat, gFull);
    boost::property_map< Graph_t, float bundled_vertex_type::* >::type 
        vertex_lon = boost::get(&bundled_vertex_type::lon, gFull);
    boost::property_map< Graph_t, int bundled_vertex_type::* >::type 
        vertex_component = boost::get(&bundled_vertex_type::component, gFull);

    auto vs = boost::vertices (gFull);
    for (auto vit = vs.first; vit != vs.second; ++vit)
    {
        if (vertex_component [*vit] == 0)
        {
            lat = vertex_lat [*vit];
            lon = vertex_lon [*vit];

            //d = std::abs (lon - lon0) + std::abs (lat - lat0);

            xd = (lon - lon0) * PI / 180.0;
            yd = (lat - lat0) * PI / 180.0;
            d = sin (yd / 2.0) * sin (yd / 2.0) + cos (lat * PI / 180.0) *
                cos (lat0 * PI / 180.0) * sin (xd / 2.0) * sin (xd / 2.0);
            d = 2.0 * atan2 (sqrt (d), sqrt (1.0 - d));
            d = d * 6371.0;

            if (d < dmin)
            {
                dmin = d;
                node = vertex_id [*vit];
            }
        }
    }

    return node;
} // end function calcDist



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                              DIJKSTRA                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Router::dijkstra (long long fromNode)
{
    std::vector<Vertex> predecessors (boost::num_vertices(gCompact)); 
    std::vector<Weight> distances (boost::num_vertices(gCompact)); 

    boost::property_map< Graph_t, long long bundled_vertex_type::* >::type 
        vertex_id = boost::get(&bundled_vertex_type::id, gCompact);
    boost::property_map< Graph_t, float bundled_vertex_type::* >::type 
        vertex_lat = boost::get(&bundled_vertex_type::lat, gCompact);
    boost::property_map< Graph_t, float bundled_vertex_type::* >::type 
        vertex_lon = boost::get(&bundled_vertex_type::lon, gCompact);

    auto p_map = boost::make_iterator_property_map
        (&predecessors[0], boost::get(boost::vertex_index, gCompact));
    auto d_map = boost::make_iterator_property_map
        (&distances[0], boost::get(boost::vertex_index, gCompact));
    auto w_map = boost::get(&bundled_edge_type::weight, gCompact); 

    int start = vertex (fromNode, gCompact);
    boost::dijkstra_shortest_paths(gCompact, start,
            weight_map(w_map). 
            predecessor_map(p_map).
            distance_map(d_map));

    // Check that all routing points have been reached
    int nvalid = 0;
    for (std::vector <RoutingPoint>::iterator itr = RoutingPointsList.begin();
            itr != RoutingPointsList.end (); itr++)
        if (distances [(*itr).nodeIndex] < FLOAT_MAX)
            nvalid++;
    std::cout << "nvalid = " << nvalid << " / " << RoutingPointsList.size ();
    assert (nvalid == RoutingPointsList.size ());

    // Trace back from each routing point 
    Vertex v0, v;
    dists.resize (0);
    float dist;

    for (std::vector <RoutingPoint>::iterator itr = RoutingPointsList.begin();
            itr != RoutingPointsList.end (); itr++)
    {
        v0 = itr->nodeIndex;
        v = v0;
        dist = 0.0;
        for (Vertex u = p_map[v]; u != v; v = u, u = p_map[v]) 
        {
            std::pair<Graph_t::edge_descriptor, bool> edgePair = 
                boost::edge(u, v, gCompact);
            Graph_t::edge_descriptor edge = edgePair.first;
            dist += boost::get (&bundled_edge_type::dist, gCompact, edge);
        }
        dists.push_back (dist);
    }

    assert (dists.size () == RoutingPointsList.size ());
    /*
     * First have to match long long fromNode to index# within
     * RoutingPointsList.  There are cases where two routing points match to
     * same OSM node, so distmats have to be filled for all such multiples
     */
    std::vector <int> id;
    id.resize (0);
    for (int i=0; i<RoutingPointsList.size(); i++)
        if (RoutingPointsList [i].nodeIndex == fromNode)
            id.push_back (i);
    assert (id.size () > 0);

    return 0;
}
