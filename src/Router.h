/***************************************************************************
 *  Project:    osm-router
 *  File:       Router.h
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


/*
 * Stores a planet-osm file as a boost::graph. The latter rely on property maps,
 * which are statically-typed entities, and so data are initially read and
 * stored as unordered_maps to be subsequently read into the graph. This initial
 * storage as unordered_maps also speeds up mapping node IDs in ways onto
 * lat-lon coordinates.
 *
 * The boost::graphs are constructed in two phases:
 * (1) A graph "gFull" is constructed from all OSM nodes and highways, and used
 * to identify both the largest connected component and all terminal or junction
 * nodes (simply as those that appear in multiple different ways). The routing
 * points are then mapped onto the nearest vertices within this largest
 * component, and these vertices are also added to terminalNodes.  (2) A second
 * reading of the data is used to make "gCompact" which only has terminalNodes
 * as vertices, and edge distances as traced along all intermediate nodes.
 */

#include "Utils.h"

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/graph/connected_components.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>


typedef std::pair <double, double> ddPair;

typedef boost::unordered_map <long long, int> umapInt;
typedef boost::unordered_map <long long, int>::iterator umapInt_Itr;
typedef boost::unordered_map <long long, ddPair> umapPair;
typedef boost::unordered_map <long long, ddPair>::iterator umapPair_Itr;

typedef std::pair <std::string, float> ProfilePair;
typedef float Weight;

struct RoutingPoint
{
    float lon, lat;
    long long node;
    int nodeIndex;
};

// See http://theboostcpplibraries.com/boost.unordered
std::size_t hash_value(const int &i)
{
    std::size_t seed = 0;
    boost::hash_combine(seed, i);
    return seed;
}
std::size_t hash_value(const ddPair &d)
{
    std::size_t seed = 0;
    boost::hash_combine(seed, d.first);
    boost::hash_combine(seed, d.second);
    return seed;
}

// edge and vertex bundles for the boost::graph
struct bundled_edge_type
{ 
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

/*
 * VertMap can be used to filter the largest connected component, but a filtered
 * graph retains all vertices regardless, so there's no real advantage here to
 * filtering.
 */
template <typename VertMap>
struct in_component_0 {
    in_component_0 () { }
    in_component_0 (VertMap comp0) : m_comp0(comp0) { }
    template <typename Vertex>
        bool operator()(const Vertex& v) const {
            return boost::get (m_comp0, v) == 0; // largest connected has id=0
        }
    VertMap m_comp0;
};

class Ways 
{
    using Graph_t = boost::adjacency_list< boost::vecS, boost::vecS, 
          boost::directedS, bundled_vertex_type, bundled_edge_type >;

    using Vertex = boost::graph_traits<Graph_t>::vertex_descriptor;

    private:
    Graph_t gFull, gCompact;
    std::string profileName;
    unsigned numWeightingProfiles=0, countWeightingProfiles;
    protected:
    float latmin, lonmin, latmax, lonmax;
    std::string _dirName;
    const std::string _city;
    const std::string osmDir = "../data/";
    const std::string profileDir = "../data/weighting_profiles/";
    std::string osmFile;
    std::vector <ProfilePair> profile;
    boost::unordered_set <long long> terminalNodes;
    dmat distMat;

    public:
    int err, count;
    long long node;
    float d;
    std::string tempstr;


    /*
     * The storage of nodeNames is done in readNodes, while they are then
     * only subsequently used to replace the long long OSM numbers with
     * corresponding ints when storing the boost::graph in the readWays
     * routine.
     */
    umapPair allNodes;
    umapInt nodeNames;
    bool firstRun = true;
    std::vector <RoutingPoint> RoutingPointsList;
    std::vector <float> dists;
    Ways (std::string str, float stdDev)
        : _city (str)
    {
        tempstr = _city;
        std::transform (tempstr.begin(), tempstr.end(), 
                tempstr.begin(), ::toupper);
        std::cout << "---" << tempstr << "---" << std::endl;
        osmFile = osmDir + "planet-" + _city + ".osm";

        boost::filesystem::path profiles (profileDir);
        boost::filesystem::directory_iterator it (profiles), eod;

        //count weighting profiles

        for (boost::filesystem::directory_iterator itCount (profiles);
        itCount != boost::filesystem::directory_iterator (); itCount++)
        {
            numWeightingProfiles++;
        }

        // initialize number generator for randomized edge weighting
        std::random_device rd;
        std::mt19937 mTwister (rd ());
        std::normal_distribution<> norm_dist (0, stdDev);

        BOOST_FOREACH (boost::filesystem::path const &p, std::make_pair (it, eod))
        {
            if (is_regular_file(p))
            {
                countWeightingProfiles++;
                profileName = p.stem ().c_str ();
                profileName = profileName.substr (profileName.find ("_") + 1);
                std::cout << "Routing weight profile: " << profileName << " ("
                    << countWeightingProfiles << "/" << numWeightingProfiles
                    << ")" << std::endl;

                setProfile (profileName.c_str ());

                // These operations are only called once
                if (firstRun)
                {
                    firstRun = false;
                    err = getBBox ();
                    err = readNodes();
                    err = readAllWays ();
                    err = getConnected ();
                    err = readRoutingPoints ();
                    distMat.resize (RoutingPointsList.size (), RoutingPointsList.size ());
                }

                // gFull is no longer needed, so can be destroyed
                gFull.clear ();

                err = readCompactWays (&mTwister, &norm_dist);
                err = remapRoutingPoints ();

                std::cout << "Getting distances between routing points";
                std::cout.flush ();
                count = 0;
                for (std::vector<RoutingPoint>::iterator itr=RoutingPointsList.begin();
                        itr != RoutingPointsList.end(); itr++)
                {
                    err = dijkstra (itr->nodeIndex);
                    assert (dists.size () == RoutingPointsList.size ());
                    std::cout << "\rGetting distances between routing points " <<
                        count << "/" << RoutingPointsList.size () << " ";
                    std::cout.flush ();
                    count++;
                }
                std::cout << "\rGetting distances between routing points " <<
                    RoutingPointsList.size () << "/" << RoutingPointsList.size
                    () << std::endl;
                std::cout << "done." << std::endl;

                writeDMat ();
                gCompact.clear ();
            }
        }
    }
    ~Ways ()
    {
        RoutingPointsList.resize (0);
        dists.resize (0);
    }

    std::string returnDirName () { return _dirName; }
    std::string returnCity () { return _city; }

    void setProfile (const std::string& profileName);
    int getBBox ();
    int readNodes ();
    int readAllWays ();
    int getConnected ();
    int readRoutingPoints ();
    int remapRoutingPoints ();
    int readCompactWays (std::mt19937* mTwister, std::normal_distribution<>* norm_dist);
    int dijkstra (long long fromNode);
    void writeDMat ();

    float calcDist (std::vector <float> x, std::vector <float> y);
    long long nearestNode (float lon0, float lat0);
    std::string getCity () { return _city; }
    std::string getCityProfile () { return _city + "_" + profileName;  }
};

void Ways::setProfile (const std::string& profileName)
{
    const std::string configfile = "../data/weighting_profiles/profile_" + \
                                    profileName + ".cfg"; 
    int ipos;
    float value;
    std::string line, field;
    std::ifstream in_file;

    Ways::profile.resize (0);

    in_file.open (configfile.c_str (), std::ifstream::in);
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

            Ways::profile.push_back (std::make_pair (field, value));
        }
    }
    in_file.close ();
};
