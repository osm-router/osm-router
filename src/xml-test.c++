/***************************************************************************
 *  Project:    osm-router
 *  File:       xml-test.c++
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

#include "xml-test.h"

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                                MAIN                                **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int main(int argc, char *argv[]) {
    std::string city;

    Xml xml (true); // returns all nodes and ways
    std::cout << "Parsed " << xml.nodes.size () << " nodes and " << 
        xml.ways.size () << " ways" << std::endl;

    return 0;
};

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                               GETBBOX                              **
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
 **                            READOVERPASS                            **
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
 **                              PARSEXML                              **
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
 **                            TRAVERSEXML                             **
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
 **                            TRAVERSEWAY                             **
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
 **                            TRAVERSENODE                            **
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
