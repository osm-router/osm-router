/***************************************************************************
 *  Project:    osm-router
 *  File:       xml-parser.h
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
 *  Description:    test of boost xml parser
 *
 *  Limitations:
 *
 *  Dependencies:       libboost
 *
 *  Compiler Options:   -std=c++11 -lboost_program_options 
 ***************************************************************************/



#include <curl/curl.h>

#include <vector>
#include <string>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>

struct Node
{
    long long id;
    float lat, lon;
};

typedef std::vector <Node> Nodes;

struct Way
{
    long long id;
    std::string type, name; // type is highway type (value for highway key)
    std::vector <std::string> key, value;
    std::vector <long long> nodes;
};

typedef std::vector <Way> Ways;


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         CLASS::CURLPLUPLUS                         **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

class CURLplusplus
{
    private:
        CURL* curl;
        std::stringstream ss;
        long http_code;
    public:
        CURLplusplus()
            : curl(curl_easy_init())
              , http_code(0)
        {
        }
        ~CURLplusplus()
        {
            if (curl) curl_easy_cleanup(curl);
        }
        std::string Get(const std::string& url)
        {
            CURLcode res;
            curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
            curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
            curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
            curl_easy_setopt(curl, CURLOPT_WRITEDATA, this);

            ss.str("");
            http_code = 0;
            res = curl_easy_perform(curl);
            if (res != CURLE_OK)
            {
                throw std::runtime_error(curl_easy_strerror(res));
            }
            curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &http_code);
            return ss.str();
        }
        long GetHttpCode()
        {
            return http_code;
        }
    private:
        static size_t write_data(void *buffer, size_t size, size_t nmemb, void *userp)
        {
            return static_cast<CURLplusplus*>(userp)->Write(buffer,size,nmemb);
        }
        size_t Write(void *buffer, size_t size, size_t nmemb)
        {
            ss.write((const char*)buffer,size*nmemb);
            return size*nmemb;
        }
};


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                             CLASS::XML                             **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

class Xml 
{
    /*
     * Downloads OSM data from the overpass API and parses the XML structure to
     * extract all nodes and ways.
     */
    private:
        const std::string _file;

    protected:
        bool file_exists;
        float latmin, lonmin, latmax, lonmax;

    public:
        std::string tempstr;
        Nodes nodes;
        Ways ways;

    Xml (std::string file)
        : _file (file)
    {
        nodes.resize (0);
        ways.resize (0);

        boost::filesystem::path p (_file);
        if (boost::filesystem::exists (p))
        {
            std::ifstream in_file;
            in_file.open (_file, std::ifstream::in);
            assert (!in_file.fail ());
            std::stringstream ss;
            ss << in_file.rdbuf ();
            tempstr = ss.str ();
        } else
        {
            getBBox ();
            std::cout << "Downloading overpass query ... ";
            std::cout.flush ();
            tempstr = readOverpass ();
            // Write raw xml data to _file:
            std::ofstream out_file;
            out_file.open (_file, std::ofstream::out);
            out_file << tempstr;
            out_file.flush ();
            out_file.close ();
            std::cout << " stored in " << _file << std::endl;
        }
        
        parseXML (tempstr);
    }
    ~Xml ()
    {
        nodes.resize (0);
        ways.resize (0);
    }

    std::string get_file () { return _file; }

    void getBBox ();
    std::string readOverpass ();
    void parseXML ( std::string & is );
    void traverseXML (const boost::property_tree::ptree& pt);
    Way traverseWay (const boost::property_tree::ptree& pt, Way way);
    Node traverseNode (const boost::property_tree::ptree& pt, Node node);
};


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         FUNCTION::GETBBOX                          **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void Xml::getBBox ()
{
    // TODO: Copy old getBBox from Graph.c++
    lonmin = -0.12;
    lonmax = -0.115;
    latmin = 51.515;
    latmax = 51.52;
}; // end function getBBox


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                       FUNCTION::READOVERPASS                       **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

std::string Xml::readOverpass ()
{
    const std::string key = "['highway']",
                url_base = "http://overpass-api.de/api/interpreter?data=";
    std::stringstream bbox, query, url;

    bbox << "";
    bbox << "(" << latmin << "," << lonmin << "," << 
        latmax << "," << lonmax << ")";

    query << "";
    query << "(node" << key << bbox.str() << ";way" << key << bbox.str() <<
                    ";rel" << key << bbox.str() << ";";
    url << "";
    url << url_base << query.str() << ");(._;>;);out;";

    CURLplusplus client;
    std::string x = client.Get (url.str().c_str ());

    return x;
}; // end function readOverpass


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         FUNCTION::PARSEXML                         **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void Xml::parseXML ( std::string & is )
{
    // populate tree structure pt
    using boost::property_tree::ptree;
    ptree pt;
    std::stringstream istream (is, std::stringstream::in);
    read_xml (istream, pt);

    // std::cout << is << std::endl; // The overpass XML data
    traverseXML (pt);
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                        FUNCTION::TRAVERSEXML                       **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void Xml::traverseXML (const boost::property_tree::ptree& pt)
{
    int tempi;
    Way way;
    Node node;

    for (boost::property_tree::ptree::const_iterator it = pt.begin ();
            it != pt.end (); ++it)
    {
        if (it->first == "node")
        {
            node = traverseNode (it->second, node);
            nodes.push_back (node);

            // ------------ just text output guff ---------------
            /*
            std::cout << "-----Node ID = " << node.id << " (" <<
                node.lat << ", " << node.lon << ")" << std::endl;
            */
        }
        if (it->first == "way")
        {
            way.key.resize (0);
            way.value.resize (0);
            way.nodes.resize (0);

            way = traverseWay (it->second, way);
            assert (way.key.size () == way.value.size ());
            // get name from key-value pairs
            for (int i=0; i<way.key.size (); i++)
                if (way.key [i] == "name")
                {
                    way.name = way.value [i];
                    tempi = i;
                    break;
                }
            way.key.erase (way.key.begin() + tempi);
            way.value.erase (way.value.begin() + tempi);
            ways.push_back (way);

            // ------------ just text output guff ---------------
            /*
            std::cout << "-----Way ID = " << way.id << std::endl;
            std::cout << "nodes = (";
            for (int i=0; i<way.nodes.size (); i++)
                std::cout << way.nodes [i] << ", ";
            std::cout << ")" << std::endl;
            for (int i=0; i<way.key.size (); i++)
                if (way.key [i] == "name")
                    std::cout << "NAME = " << way.value [i] << std::endl;
                else if (way.key [i] == "highway")
                    std::cout << "highway: " << way.value [i] << std::endl;
            */
        } else
            traverseXML (it->second);
    }
}

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                        FUNCTION::TRAVERSEWAY                       **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

Way Xml::traverseWay (const boost::property_tree::ptree& pt, Way way)
{
    for (boost::property_tree::ptree::const_iterator it = pt.begin ();
            it != pt.end (); ++it)
    {
        if (it->first == "k")
            way.key.push_back (it->second.get_value <std::string> ());
        else if (it->first == "v")
            way.value.push_back (it->second.get_value <std::string> ());
        else if (it->first == "id")
            way.id = it->second.get_value <long long> ();
        else if (it->first == "ref")
            way.nodes.push_back (it->second.get_value <long long> ());
        way = traverseWay (it->second, way);
    }

    return way;
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                       FUNCTION::TRAVERSENODE                       **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

Node Xml::traverseNode (const boost::property_tree::ptree& pt, Node node)
{
    for (boost::property_tree::ptree::const_iterator it = pt.begin ();
            it != pt.end (); ++it)
    {
        if (it->first == "id")
            node.id = it->second.get_value <long long> ();
        else if (it->first == "lat")
            node.lat = it->second.get_value <float> ();
        else if (it->first == "lon")
            node.lon = it->second.get_value <float> ();
        // No other key-value pairs currently extracted for nodes
        node = traverseNode (it->second, node);
    }

    return node;
}
