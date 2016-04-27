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



#include "Utils.h"
#include <curl/curl.h>

#include <boost/filesystem.hpp>

#include <vector>
#include <string>
#include <fstream>
#include <locale>
#include <boost/foreach.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>

struct ptree
{
    std::string data;                          // data associated with the node
    std::list< std::pair<std::string, ptree> > children; // ordered list of named children
};

struct Way
{
    long long waynodes;
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
    const std::string _city;
    float latmin, lonmin, latmax, lonmax;
    std::string _dirName;

    public:
    int err, count;
    long long node;
    float d;
    std::string tempstr;
    Ways ways;


    /*
     * The storage of nodeNames is done in readNodes, while they are then
     * only subsequently used to replace the long long OSM numbers with
     * corresponding ints when storing the boost::graph in the readWays
     * routine.
     */
    Xml (std::string str)
        : _city (str)
    {
        tempstr = _city;
        std::transform (tempstr.begin(), tempstr.end(), 
                tempstr.begin(), ::toupper);

        err = getBBox ();
        std::cout << "downloading overpass query ... ";
        std::cout.flush ();
        tempstr = readOverpass ();
        std::cout << "done; processing ..." << std::endl;
        ways = readNodes (tempstr);
    }
    ~Xml ()
    {
    }

    std::string returnDirName () { return _dirName; }
    std::string returnCity () { return _city;   }

    int getBBox ();
    std::string readOverpass ();
    Ways readNodes ( std::string & is );
    void traverse (const boost::property_tree::ptree& pt);
    void traverseWay (const boost::property_tree::ptree& pt);
    void traverseNode (const boost::property_tree::ptree& pt);
};

