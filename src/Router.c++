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
 *  NeutralClusters.  If not, see <http://www.gnu.org/licenses/>.
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

#include "Router.h"

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                                MAIN                                **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int main(int argc, char *argv[]) {
    std::string city;
    float stdDev;

    try {
        boost::program_options::options_description generic("Generic options");
        generic.add_options()
            ("version,v", "print version std::string")
            ("help", "produce help message")    
            ;

        boost::program_options::options_description config("Configuration");
        config.add_options()
            ("city,c", boost::program_options::value <std::string> 
                (&city)->default_value ("london"), "city")
            ("stdDev,s", boost::program_options::value <float>
                (&stdDev)->default_value (0), "standard deviation of edge weights")
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

    std::transform (city.begin(), city.end(), city.begin(), ::tolower);
    if (city.substr (0, 2) == "lo")
        city = "london";
    else if (city.substr (0, 2) == "ny")
        city = "nyc";

    Ways ways (city, stdDev);
};


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                               GETBBOX                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Ways::getBBox ()
{
    const float buffer = 0.01;
    float lat, lon;
    int nskips, ipos = 0;

    lonmin = latmin = FLOAT_MAX;
    lonmax = latmax = -FLOAT_MAX;

    std::string linetxt, txt, fname;
    fname = "../data/routing_points_" + getCity () + ".txt";
    std::ifstream in_file;
    
    in_file.open (fname.c_str (), std::ifstream::in);
    assert (!in_file.fail ());

    in_file.clear ();
    in_file.seekg (0); 
    getline (in_file, linetxt, '\n'); // header

    while (getline (in_file, linetxt, '\n'))
    {
        ipos = linetxt.find (",");
        linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
        ipos = linetxt.find (",");

        lat = atof (linetxt.substr (0, ipos).c_str());
        if (lat < latmin)
            latmin = lat;
        else if (lat > latmax)
            latmax = lat;

        linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
        ipos = linetxt.find (",");
        lon = atof (linetxt.substr (0, ipos).c_str());
        if (lon < lonmin)
            lonmin = lon;
        else if (lon > lonmax)
            lonmax = lon;
    } // end while getline
    in_file.close ();

    lonmin -= buffer;
    latmin -= buffer;
    lonmax += buffer;
    latmax += buffer;

    return 0;
}; // end function getBBox


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                              READNODES                             **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Ways::readNodes ()
{
    int ipos;
    long long node;
    float lat, lon;
    std::string linetxt, txt;
    std::ifstream in_file;

    /*
     * Note that boost::iostreams::bzip2_decompressor could be used to unzip the
     * .bz2, but they don't seem to compile, and it's not for lack of the
     * appropriate libraries. The compiler suggests an internal error in the
     * decompressor.
     *
     * TODO: Check this out further
     */

    in_file.open (osmFile.c_str (), std::ifstream::in);
    assert (!in_file.fail ());
    in_file.clear ();
    in_file.seekg (0); 

    std::cout << "Reading nodes ...";
    std::cout.flush ();

    float nodeLonmin = FLOAT_MAX, nodeLatmin = FLOAT_MAX,
        nodeLonmax = -FLOAT_MAX, nodeLatmax = -FLOAT_MAX;

    while (getline (in_file, linetxt, '\n'))
    {
        if (linetxt.find ("<node") != std::string::npos)
        {
            ipos = linetxt.find ("id=\"");
            linetxt = linetxt.substr (ipos + 4, 
                    linetxt.length () - ipos - 4);
            ipos = linetxt.find ("\"");
            node = atoll (linetxt.substr (0, ipos).c_str());

            // Note that this presumes that lat always comes before lon!
            assert (linetxt.find ("lat=\"") < linetxt.find ("lon=\""));
            ipos = linetxt.find ("lat=\"");
            linetxt = linetxt.substr (ipos + 5, linetxt.length () - ipos - 5);
            ipos = linetxt.find ("\"");
            lat = atof (linetxt.substr (0, ipos).c_str());

            if (lat < nodeLatmin)
                nodeLatmin = lat;
            else if (lat > nodeLatmax)
                nodeLatmax = lat;

            if (lat >= latmin && lat <= latmax)
            {
                ipos = linetxt.find ("lon=\"");
                linetxt = linetxt.substr (ipos + 5, linetxt.length () - ipos - 5);
                ipos = linetxt.find ("\"");
                lon = atof (linetxt.substr (0, ipos).c_str());

                // lon ranges are only calculated within BBOX of RoutingPoints
                if (lon < nodeLonmin)
                    nodeLonmin = lon;
                else if (lon > nodeLonmax)
                    nodeLonmax = lon;

                if (lon >= lonmin && lon <= lonmax)
                    allNodes [node] = std::make_pair (lat, lon);
            }
        }
    } 
    assert (nodeLonmin < lonmin);
    assert (nodeLonmax > lonmax);
    assert (nodeLatmin < latmin);
    assert (nodeLatmax > latmax);

    /*
     * Then rewind file and read any extra nodes that are part of ways with
     * those in allNodes.
     */
    in_file.clear ();
    in_file.seekg (0);

    bool inway = false, highway = false, nodeFound;
    int id0, id1, nways = 0;
    umapPair_Itr umapitr;
    std::vector <long long> waynodes;
    boost::unordered_set <long long> extraNodes;
    
    while (getline (in_file, linetxt, '\n'))
    {
        if (linetxt.find ("<way") != std::string::npos)
        {
            inway = true;
            highway = false;
            nodeFound = false;
            waynodes.resize (0);
        }
        else if (linetxt.find ("</way>") != std::string::npos)
        {
            if (highway & nodeFound)
            {
                for (std::vector <long long>::iterator itr=waynodes.begin ();
                        itr != waynodes.end (); itr++)
                    if (allNodes.find (*itr) == allNodes.end () &&
                        extraNodes.find (*itr) == extraNodes.end ())
                        extraNodes.insert (*itr);
            } // end if highway
            inway = false;
            nodeFound = false;
        } // end else if end way
        else if (inway)
        {
            if (linetxt.find ("<nd") != std::string::npos)
            {
                ipos = linetxt.find ("<nd ref=");
                linetxt = linetxt.substr (ipos + 9, 
                        linetxt.length () - ipos - 9);
                node = atoll (linetxt.c_str ());
                waynodes.push_back (node);
                if (allNodes.find (node) != allNodes.end ())
                    nodeFound = true;
            }
            else if (linetxt.find ("k=\"highway\"") != std::string::npos)
            {
                // highway is only true if it has one of the values listed in
                // the profile
                for (std::vector<ProfilePair>::iterator itr = profile.begin();
                        itr != profile.end(); itr++)
                {
                    tempstr = "v=\"" + (*itr).first;
                    if (linetxt.find (tempstr) != std::string::npos)
                    {
                        highway = true;
                        break;
                    }
                }
            }
        } // end else if inway
    } // end while getline

    /*
     * And then finally store the extraNodes in allNodes
     */
    in_file.clear ();
    in_file.seekg (0);

    while (getline (in_file, linetxt, '\n'))
    {
        if (linetxt.find ("<node") != std::string::npos)
        {
            ipos = linetxt.find ("id=\"");
            linetxt = linetxt.substr (ipos + 4, 
                    linetxt.length () - ipos - 4);
            ipos = linetxt.find ("\"");
            node = atoll (linetxt.substr (0, ipos).c_str());
            if (extraNodes.find (node) != extraNodes.end ())
            {
                ipos = linetxt.find ("lat=\"");
                linetxt = linetxt.substr (ipos + 5, linetxt.length () - ipos - 5);
                ipos = linetxt.find ("\"");
                lat = atof (linetxt.substr (0, ipos).c_str());

                ipos = linetxt.find ("lon=\"");
                linetxt = linetxt.substr (ipos + 5, linetxt.length () - ipos - 5);
                ipos = linetxt.find ("\"");
                lon = atof (linetxt.substr (0, ipos).c_str());

                allNodes [node] = std::make_pair (lat, lon);
            }
        }
    } 
    in_file.close ();

    std::cout << "\rRead coordinates of " << allNodes.size () << " nodes." <<
        std::endl;

    return allNodes.size ();
} // end Way::readNodes


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                            READALLWAYS                             **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Ways::readAllWays ()
{
    /*
     * The sole point of this routine is to make the boost::graph "gFull", which
     * is in turn used for the three purposes of (1) determining the largest
     * connected component, (2) enabling routing points to be allocated to
     * the nearest highway vertex that is part of this connected component, and
     * (3) making a list of terminal nodes which are any that appear in multiple
     * ways plus the routing point nodes.
     *
     * The lat-lons of vertices are needed in gFull for matching routing points,
     * but edge distances are not relevant, so default to 1.0.
     */
    bool inBBox, inway = false, highway = false, oneway;
    int ipos, id0, id1, nodeCount = 0, nways = 0;
    long long node;
    float lat0, lon0, lat1, lon1;
    umapPair_Itr umapitr;
    boost::unordered_set <long long> nodeList;
    std::vector <long long> waynodes;
    std::string linetxt, tempstr;
    std::ifstream in_file;
    
    in_file.open (osmFile.c_str (), std::ifstream::in);
    assert (!in_file.fail ());
    in_file.clear ();
    in_file.seekg (0); 

    // oneways should only be v="yes", but wiki allows the other two as well
    typedef std::vector <std::string> strvec;
    strvec oneWayList;
    oneWayList.push_back ("k=\"oneway\" v=\"yes");
    oneWayList.push_back ("k=\"oneway\" v=\"1");
    oneWayList.push_back ("k=\"oneway\" v=\"true");

    /*
     * The boost::graph is only given size through the following vertex
     * additions. Vertex numbers are implicit and sequential, enumeratred by
     * nodeCount, and edges must reference these numbers.
     */
    bundled_vertex_type oneVert;
    bundled_edge_type oneEdge;

    std::cout << "Reading ways ...";
    std::cout.flush ();

    // Can be used to plot a map of the network:
    //std::ofstream out_file;
    //out_file.open ("map.csv", std::ofstream::out);

    while (getline (in_file, linetxt, '\n'))
    {
        if (linetxt.find ("<way") != std::string::npos)
        {
            inway = true;
            highway = false;
            oneway = false;
            inBBox = false;
            waynodes.resize (0);
        }
        else if (linetxt.find ("</way>") != std::string::npos)
        {
            /*
             * Only highways that have every waynode in allNodes are inBBox, so
             * this first has to be double checked.
             */
            if (highway && inBBox)
            {
                for (std::vector <long long>::iterator itr=waynodes.begin();
                        itr != waynodes.end(); itr++)
                    if (allNodes.find (*itr) == allNodes.end ())
                        inBBox = false;
            }
            if (highway && inBBox)
            {
                node = waynodes.front ();
                assert ((umapitr = allNodes.find (node)) != allNodes.end());
                lat0 = ((*umapitr).second).first;
                lon0 = ((*umapitr).second).second;

                // Presume all first nodes are terminal nodes
                if (terminalNodes.find (node) == terminalNodes.end())
                    terminalNodes.insert (node);

                if (nodeNames.find (node) == nodeNames.end())
                {
                    id0 = nodeCount++;
                    nodeNames [node] = id0;
                    oneVert.id = node;
                    oneVert.lat = lat0;
                    oneVert.lon = lon0;
                    boost::add_vertex (oneVert, gFull);
                }
                else
                    id0 = (*nodeNames.find (node)).second;

                // Note std::next (c++11) in loop
                for (std::vector <long long>::iterator 
                        itr=std::next (waynodes.begin());
                        itr != waynodes.end(); itr++)
                {
                    node = (*itr);
                    assert ((umapitr = allNodes.find (node)) != allNodes.end());
                    lat1 = ((*umapitr).second).first;
                    lon1 = ((*umapitr).second).second;
                    //out_file << lon0 << "," << lat0 << "," <<
                    //    lon1 << "," << lat1 << std::endl;

                    if (nodeList.find (node) == nodeList.end())
                        nodeList.insert (node);
                    else if (terminalNodes.find (node) == terminalNodes.end())
                        terminalNodes.insert (node);

                    if (nodeNames.find (node) == nodeNames.end())
                    {
                        id1 = nodeCount++;
                        nodeNames [node] = id1;
                        oneVert.id = node;
                        oneVert.lat = lat1;
                        oneVert.lon = lon1;
                        boost::add_vertex (oneVert, gFull);
                    }
                    else
                        id1 = (*nodeNames.find (node)).second;

                    oneEdge.weight = oneEdge.dist = 1.0;
                    boost::add_edge(id0, id1, oneEdge, gFull);
                    if (!oneway)
                        boost::add_edge(id1, id0, oneEdge, gFull);

                    nways++;
                    id0 = id1;
                    lat0 = lat1;
                    lon0 = lon1;
                }
            } // end if highway
            inway = false;
            inBBox = false;
        } // end else if end way
        else if (inway)
        {
            /*
             * Nodes have to first be stored, because they are only subsequently
             * analysed if they are part of a highway.
             */
            if (linetxt.find ("<nd") != std::string::npos)
            {
                ipos = linetxt.find ("<nd ref=");
                linetxt = linetxt.substr (ipos + 9, 
                        linetxt.length () - ipos - 9);
                node = atoll (linetxt.c_str ());
                waynodes.push_back (atoll (linetxt.c_str ()));
                if (allNodes.find (node) != allNodes.end ())
                    inBBox = true;
            }
            else if (linetxt.find ("k=\"highway\"") != std::string::npos)
            {
                /* 
                 * highway is only true if it has one of the values listed in
                 * the profile
                 */
                for (std::vector<ProfilePair>::iterator itr = profile.begin();
                        itr != profile.end(); itr++)
                {
                    tempstr = "v=\"" + (*itr).first;
                    if (linetxt.find (tempstr) != std::string::npos)
                    {
                        highway = true;
                        break;
                    }
                    for (strvec::iterator oit = oneWayList.begin();
                            oit != oneWayList.end(); oit++)
                        if (linetxt.find (*oit) != std::string::npos)
                            oneway = true;
                }
            }
        } // end else if inway
    } // end while getline
    in_file.close ();
    oneWayList.resize (0);
    //out_file.close ();

    /*
     * nodes and corresponding indices can then be respectively obtained with
     * node = nodeNames.find (node)->first; // redundant
     * index = nodeNames.find (node)->second;
     */

    std::vector<int> compvec(num_vertices(gFull));
    ipos = boost::connected_components(gFull, &compvec[0]);
    std::cout << "\rRead " << nways << " ways with " << 
        num_vertices (gFull) << " vertices and " << ipos << 
        " components" << std::endl;

    return 0;
} // end function readAllWays



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          READCOMPACTWAYS                           **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Ways::readCompactWays (std::mt19937& mTwister, std::normal_distribution<>& norm_dist)
{
    bool inBBox, inway = false, highway = false, oneway;
    int ipos, id0, id1, nodeCount = 0, nways = 0;
    long long node;
    float d, weight, tempWeight;
    umapPair_Itr umapitr;
    boost::unordered_set <long long> nodeList;
    std::vector <long long> waynodes;
    std::vector <float> lats, lons;
    std::string linetxt, tempstr, highwayType;
    std::ifstream in_file;
    
    in_file.open (osmFile.c_str (), std::ifstream::in);
    assert (!in_file.fail ());
    in_file.clear ();
    in_file.seekg (0); 

    // oneways should only be "yes", but wiki allows the other two as well
    typedef std::vector <std::string> strvec;
    strvec oneWayList;
    oneWayList.push_back ("k=\"oneway\" v=\"yes");
    oneWayList.push_back ("k=\"oneway\" v=\"0");
    oneWayList.push_back ("k=\"oneway\" v=\"true");

    bundled_vertex_type oneVert;
    bundled_edge_type oneEdge;

    nodeNames.clear ();

    std::cout << "Reading compact ways ...";
    std::cout.flush ();

    while (getline (in_file, linetxt, '\n'))
    {
        if (linetxt.find ("<way") != std::string::npos)
        {
            inway = true;
            highway = false;
            oneway = false;
            inBBox = false;
            weight = -FLOAT_MAX;
            waynodes.resize (0);
        }
        else if (linetxt.find ("</way>") != std::string::npos)
        {
            if (highway && inBBox)
            {
                for (std::vector <long long>::iterator itr=waynodes.begin();
                        itr != waynodes.end(); itr++)
                    if (allNodes.find (*itr) == allNodes.end ())
                        inBBox = false;
            }
            if (highway && inBBox)
            {
                lats.resize (0);
                lons.resize (0);

                node = waynodes.front ();

                if (nodeNames.find (node) == nodeNames.end())
                {
                    id0 = nodeCount++;
                    nodeNames [node] = id0;
                    umapitr = allNodes.find (node);
                    oneVert.lat = ((*umapitr).second).first;
                    oneVert.lon = ((*umapitr).second).second;
                    oneVert.id = node;
                    boost::add_vertex (oneVert, gCompact);
                }
                else
                    id0 = (*nodeNames.find (node)).second;

                for (std::vector <long long>::iterator itr=waynodes.begin();
                        itr != waynodes.end(); itr++)
                {
                    assert ((umapitr = allNodes.find (*itr)) != allNodes.end());
                    lats.push_back (((*umapitr).second).first);
                    lons.push_back (((*umapitr).second).second);
                    /*
                     * If a single way crosses a terminal node, then break it
                     * into two separate ways either side thereof.
                     */
                    if (itr > waynodes.begin() &&
                            (terminalNodes.find (*itr)) != terminalNodes.end())
                    {
                        node = (*itr);
                        if (nodeNames.find (node) == nodeNames.end())
                        {
                            id1 = nodeCount++;
                            nodeNames [node] = id1;
                            oneVert.id = node;
                            oneVert.lat = ((*umapitr).second).first;
                            oneVert.lon = ((*umapitr).second).second;
                            boost::add_vertex (oneVert, gCompact);
                        }
                        else
                            id1 = (*nodeNames.find (node)).second;

                        oneEdge.dist = calcDist (lons, lats);
                        if (weight == 0.0)
                            oneEdge.weight = FLOAT_MAX;
                        else
                            oneEdge.weight = oneEdge.dist / weight;
                            /*
                             * A normally distributed random factor with mean 0
                             * and a user defined standard deviation is added to
                             * the weight. If this results in a negative edge
                             * weight, a new factor is calculated.
                             */
                            tempWeight = oneEdge.weight + norm_dist (mTwister);
                            while (tempWeight < 0)
                            {
                                tempWeight = oneEdge.weight + norm_dist (mTwister);
                            }
                            oneEdge.weight = tempWeight;

                        boost::add_edge (id0, id1, oneEdge, gCompact);
                        if (!oneway)
                            boost::add_edge(id1, id0, oneEdge, gCompact);

                        lats.resize (0);
                        lons.resize (0);
                        lats.push_back (((*umapitr).second).first);
                        lons.push_back (((*umapitr).second).second);

                        nways++;
                        id0 = id1;
                    }
                }
            } // end if highway
            inway = false;
        } // end else if highway
        else if (inway)
        {
            if (linetxt.find ("<nd") != std::string::npos)
            {
                ipos = linetxt.find ("<nd ref=");
                linetxt = linetxt.substr (ipos + 9, 
                        linetxt.length () - ipos - 9);
                node = atoll (linetxt.c_str ());
                waynodes.push_back (atoll (linetxt.c_str ()));
                if (allNodes.find (node) != allNodes.end ())
                    inBBox = true;
            }
            else if (linetxt.find ("k=\"highway\"") != std::string::npos)
            {
                for (std::vector<ProfilePair>::iterator itr = profile.begin();
                        itr != profile.end(); itr++)
                {
                    tempstr = "v=\"" + (*itr).first;
                    if (linetxt.find (tempstr) != std::string::npos)
                    {
                        highway = true;
                        highwayType = (*itr).first;
                        weight = (*itr).second;
                        break;
                    }
                    for (strvec::iterator oit = oneWayList.begin();
                            oit != oneWayList.end(); oit++)
                        if (linetxt.find (*oit) != std::string::npos)
                            oneway = true;
                }
            }
        } // end else if inway
    } // end while getline
    in_file.close ();
    oneWayList.resize (0);

    std::vector<int> compvec(num_vertices(gCompact));
    ipos = boost::connected_components(gCompact, &compvec[0]);
    std::cout << "\rRead " << nways << " compact ways with " << 
        num_vertices (gCompact) << " vertices and " << ipos << 
        " components" << std::endl;

    return 0;
} // end function readCompactWays


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                              CALCDIST                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/


float Ways::calcDist (std::vector <float> x, std::vector <float> y)
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
} // end function calcDist


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                            GETCONNECTED                            **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Ways::getConnected ()
{
    std::vector<int> compvec(num_vertices(gFull));
    int num = boost::connected_components(gFull, &compvec[0]);

    // Then store component info in vertices
    /*
    typedef boost::graph_traits <Graph_t>::vertex_iterator viter;
    std::pair <viter, viter> vp;
    for (vp = vertices(gFull); vp.first != vp.second; ++vp.first)
        vertex_component [*vp.first] = compvec [*vp.first];
     */

    // Alternative:
    boost::property_map< Graph_t, int bundled_vertex_type::* >::type 
        vertex_component = boost::get(&bundled_vertex_type::component, gFull);
    auto vs = boost::vertices (gFull);
    for (auto vit = vs.first; vit != vs.second; ++vit)
        vertex_component [*vit] = compvec [*vit];

    // Optional filtering of component = 0:
    /*
    in_component_0 <VertMap> filter (boost::get 
        (&bundled_vertex_type::component, gFull));
    boost::filtered_graph <Graph, in_component_0 <VertMap> > fg (g, filter);
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
 **                         READROUTINGPOINTS                          **
 **                                                                    **
 ************************************************************************
 ************************************************************************/


int Ways::readRoutingPoints ()
{
    int station_id, ipos = 0;
    std::string linetxt, txt, fname, stationNodeFileName;

    stationNodeFileName = "routing_point_nodes_" + getCity () + ".txt";
    fname = "../data/routing_points_" + getCity () + ".txt";

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

            if (terminalNodes.find (rpoint.node) == terminalNodes.end())
                terminalNodes.insert (rpoint.node);
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
            assert (nodeNames.find (rpoint.node) != nodeNames.end());
            rpoint.nodeIndex = nodeNames.find (rpoint.node)->second;
            RoutingPointsList.push_back (rpoint);

            stationNodeFile << station_id << ", " << rpoint.node << ", " << rpoint.nodeIndex << '\n';

            // Also add rpoint nodes to terminalNodes
            if (terminalNodes.find (rpoint.node) == terminalNodes.end())
                terminalNodes.insert (rpoint.node);
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


long long Ways::nearestNode (float lon0, float lat0)
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
 **                        REMAPROUTINGPOINTS                          **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Ways::remapRoutingPoints ()
{
    /*
     * Dijkstra needs the vertex index number, which now has to changed from the
     * index into gFull to the index into gCompact
     */

    boost::property_map< Graph_t, long long bundled_vertex_type::* >::type 
        vertex_id = boost::get(&bundled_vertex_type::id, gCompact);

    auto vs = boost::vertices (gCompact);

    for (std::vector <RoutingPoint>::iterator sitr = RoutingPointsList.begin();
            sitr != RoutingPointsList.end(); sitr++)
    {
        (*sitr).nodeIndex = -INT_MAX;
        for (auto vit = vs.first; vit != vs.second; ++vit)
            if (vertex_id [*vit] == (*sitr).node)
            {
                (*sitr).nodeIndex = *vit;
                break;
            }
        assert ((*sitr).nodeIndex >= 0);
    }

    return 0;
}

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                              DIJKSTRA                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int Ways::dijkstra (long long fromNode)
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
    for (int i=0; i<id.size (); i++)
        idDone [id[i]] = true;

    for (int i=0; i<RoutingPointsList.size (); i++)
        for (std::vector <int>::iterator itr=id.begin(); itr != id.end(); itr++)
            distMat (*itr, i) = dists [i];

    return 0;
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                             WRITEDMAT                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void Ways::writeDMat ()
{
    std::string fname = "RoutingPointDistsMat-" + getCityProfile () + ".csv";

    std::ofstream out_file;

    out_file.open (fname, std::ios::out);

    for (int i=0; i<distMat.size1(); i++)
    {
        for (int j=0; j<(distMat.size2() - 1); j++)
            out_file << distMat (i, j) << ", ";
        out_file << distMat (i, distMat.size2() - 1) << std::endl;
    }

    out_file.close ();
};
