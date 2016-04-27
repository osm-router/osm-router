/***************************************************************************
 *  Project:    osm-router
 *  File:       xml-test.h
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

class Xml 
{
    private:

    protected:
        float latmin, lonmin, latmax, lonmax;

    public:
        std::string tempstr;
        Nodes nodes;
        Ways ways;

    Xml ()
    {
        nodes.resize (0);
        ways.resize (0);
        getBBox ();
        std::cout << "Downloading overpass query ... ";
        std::cout.flush ();
        tempstr = readOverpass ();
        parseXML (tempstr);
        std::cout << "Parsed " << nodes.size () << " nodes and " << 
            ways.size () << " ways" << std::endl;
    }
    ~Xml ()
    {
        nodes.resize (0);
        ways.resize (0);
    }

    void getBBox ();
    std::string readOverpass ();
    void parseXML ( std::string & is );
    void traverseXML (const boost::property_tree::ptree& pt);
    Way traverseWay (const boost::property_tree::ptree& pt, Way way);
    Node traverseNode (const boost::property_tree::ptree& pt, Node node);
};

