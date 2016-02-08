/***************************************************************************
 *  Project:    osm-router
 *  File:       StationData.c++
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
 *  Limitations:
 *
 *  Dependencies:       libboost
 *
 *  Compiler Options:   -std=c++11 -lboost_program_options 
 ***************************************************************************/

#include "StationData.h"

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                           GETDIRNAME                               **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

std::string Stations::GetDirName ()
{
    std::ifstream in_file;
    std::string dirtxt;
    std::string configfile = "../bikes.cfg"; // Contains name of data directory
    struct dirent *ent;

    in_file.open (configfile.c_str (), std::ifstream::in);
    assert (!in_file.fail ());

    while (!in_file.eof ())
    {
        getline (in_file, dirtxt, '\n');
        if (dirtxt.find (':') < std::string::npos) 
        {
            dirtxt = dirtxt.substr (0, dirtxt.find (':'));
            std::transform (dirtxt.begin(), dirtxt.end(), 
                    dirtxt.begin(), ::tolower);
            if (dirtxt.substr (0, 3) == _city.substr (0, 3)) 
            {
                getline (in_file, dirtxt, '\n');
                break;
            }
        }
    }
    in_file.close ();
    return dirtxt;
} // end Stations::GetDirName

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                           GETDIRLIST                               **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void StationData::GetDirList ()
{
    std::string fname;
    DIR *dir;
    struct dirent *ent;

    filelist.resize (0);
    if ((dir = opendir (_dirName.c_str())) != NULL) 
    {
        while ((ent = readdir (dir)) != NULL) 
        {
            fname = ent->d_name;
            if (fname != "." && fname != "..") 
            {
                fname = fname;
                filelist.push_back (fname);
            }
        }
        closedir (dir);
        std::sort (filelist.begin(), filelist.end());
    } else {
        std::string outstr = "ERROR: Directory for city = " +\
                              _city + " at " + _dirName + " does not exist";
        perror ("");
        std::cout << outstr << std::endl;
        //return EXIT_FAILURE;
    }
} // end StationData::GetDirList


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          GETSTATIONS                               **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int StationData::GetStations ()
{
    /*
     * Reads from station_latlons which is constructed with getLatLons.py and
     * *MUST* be ordered numerically. Returns _numStations
     *
     * TODO: Write handlers for cases where there are trips to/from stations not
     * in StationList.
     * TODO: Update StationList for nyc, because there seem to be trips to/from
     * one station that is not in list.
     * 
     */
    const std::string dir = "data/"; 
    int ipos, tempi, count;
    bool tube;
    OneStation oneStation;
    std::string fname;
    std::ifstream in_file;
    std::string linetxt;

    StationList.resize (0);
    count = 0;
    oneStation.name = "";

    if (_city == "london" || _city == "nyc" || _city == "washingtondc")
    {
        fname = dir + "station_latlons_" + _city + ".txt";
        in_file.open (fname.c_str (), std::ifstream::in);
        assert (!in_file.fail ());
        in_file.clear ();
        in_file.seekg (0); 
        getline (in_file, linetxt, '\n'); // header
        while (!in_file.eof ()) 
        {
            getline (in_file, linetxt,'\n');
            if (linetxt.length () > 1) 
            {
                ipos = linetxt.find(',',0);
                tempi = atoi (linetxt.substr (0, ipos).c_str());
                if (tempi > count) 
                    count = tempi;
                oneStation.ID = tempi;
                linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
                ipos = linetxt.find (',', 0);
                oneStation.lat = atof (linetxt.substr (0, ipos).c_str());
                linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
                ipos = linetxt.find (',', 0);
                oneStation.lon = atof (linetxt.substr (0, ipos).c_str());
                StationList.push_back (oneStation);
            }
        }
        in_file.close();
    } 
    else if (_city == "boston")
    {
        fname = dir + "hubway_stations.csv";
        in_file.open (fname.c_str (), std::ifstream::in);
        assert (!in_file.fail ());
        in_file.clear ();
        in_file.seekg (0); 
        getline (in_file, linetxt, '\n'); // header
        while (getline (in_file, linetxt,'\n'))
        {
            ipos = linetxt.find(',',0);
            tempi = atoi (linetxt.substr (0, ipos).c_str());
            if (tempi > count) 
                count = tempi;
            oneStation.ID = tempi;
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            for (int i=0; i<3; i++)
            {
                ipos = linetxt.find (',', 0);
                linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            }
            ipos = linetxt.find (',', 0);
            oneStation.lat = atof (linetxt.substr (0, ipos).c_str());
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            ipos = linetxt.find (',', 0);
            oneStation.lon = atof (linetxt.substr (0, ipos).c_str());
            StationList.push_back (oneStation);
        }
        in_file.close();
    }
    else if (_city == "chicago")
    {
        fname = dir + "Divvy_Stations_2014-Q3Q4.csv";
        in_file.open (fname.c_str (), std::ifstream::in);
        assert (!in_file.fail ());
        in_file.clear ();
        in_file.seekg (0); 
        getline (in_file, linetxt, '\n'); // header

        while (getline (in_file, linetxt,'\n'))
        {
            ipos = linetxt.find(',',0);
            tempi = atoi (linetxt.substr (0, ipos).c_str());
            if (tempi > count) 
                count = tempi;
            oneStation.ID = tempi;
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            ipos = linetxt.find (',', 0);
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            ipos = linetxt.find (',', 0);
            oneStation.lat = atof (linetxt.substr (0, ipos).c_str());
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            ipos = linetxt.find (',', 0);
            oneStation.lon = atof (linetxt.substr (0, ipos).c_str());
            StationList.push_back (oneStation);
        }
        in_file.close();
    }
    else if (_city == "oysterRail" || _city == "oysterTube")
    {
        oneStation.ID = INT_MIN;
        if (_city == "oysterRail")
        {
            tube = false;
            fname = dir + "London-rail-stations.txt";
        }
        else
        {
            tube = true;
            fname = dir + "London-tube-stations.txt";
        }
        in_file.open (fname.c_str (), std::ifstream::in);
        assert (!in_file.fail ());
        in_file.clear ();
        in_file.seekg (0); 
        getline (in_file, linetxt,'\n');
        while (!in_file.eof ()) 
        {
            getline (in_file, linetxt,'\n');
            if (linetxt.length () > 1) 
            {
                ipos = linetxt.find(',',0);
                oneStation.name = standardise (linetxt.substr (0, ipos).c_str());
                oneStation.name = substituteNames (tube, oneStation.name);
                linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
                ipos = linetxt.find (',', 0);
                oneStation.lat = atof (linetxt.substr (0, ipos).c_str());
                linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
                oneStation.lon = atof (linetxt.c_str());
                StationList.push_back (oneStation);
            }
        }
        in_file.close();
    }

    return count;
} // end StationData::GetStations


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          MAKESTATIONINDEX                          **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void StationData::MakeStationIndex ()
{
    // First station is #1 and last is _maxStation, so _StationIndex has 
    // len (_maxStns + 1), with _StationIndex [sti.ID=1] = 0 and
    // _StationIndex [sti.ID=_maxStation] = _numStations.
    OneStation sti;

    _StationIndex.resize (_maxStation + 1);
    for (std::vector <int>::iterator pos=_StationIndex.begin();
            pos != _StationIndex.end(); pos++)
        *pos = INT_MIN;
    for (int i=0; i<StationList.size (); i++) 
    {
        sti = StationList [i];
        _StationIndex [sti.ID] = i;
    }
} // end StationData::MakeStationIndex


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                             COUNTTRIPS                             **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

double StationData::CountTrips ()
{
    double count = 0.0;

    for (int i=0; i<ntrips.size1(); i++)
        for (int j=0; j<ntrips.size2(); j++)
            count += ntrips (i, j);

    return (count);
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                              READDMAT                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int StationData::readDMat ()
{
    std::string distFile;
    if (_city == "london")
        distFile = "data/station_dists_london.txt";
    else if (_city == "nyc")
        distFile = "data/station_dists_nyc.txt";
    else if (_city == "boston" || _city == "chicago" || _city == "washingtondc")
        distFile = "results/stationDistsMat-" + _city + ".csv";
    else if (_city == "oysterTube")
        distFile = "data/London-tube-station-dists.txt";
    else if (_city == "oysterRail")
        distFile = "data/London-rail-station-dists.txt";

    int count, ipos, tempi [2];
    double d;
    std::ifstream in_file;
    std::string linetxt;
    in_file.open (distFile.c_str (), std::ifstream::in);
    if (in_file.fail ()) {
        // TODO: INSERT ERROR HANDLER
        return -1;
    } 
    in_file.clear ();
    in_file.seekg (0); 
    count = 0;
    
    if (_city == "london" || _city == "nyc")
    {
        while (getline (in_file, linetxt, '\n')) {
            ipos = linetxt.find(',',0);
            tempi [0] = atoi (linetxt.substr (0, ipos).c_str()); // Start Station ID
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            ipos = linetxt.find(',',0);
            tempi [1] = atoi (linetxt.substr (0, ipos).c_str()); // End Station ID
            tempi [0] = _StationIndex [tempi[0]];
            tempi [1] = _StationIndex [tempi[1]];
            if (tempi[0] < 0 | tempi[0] > returnNumStations ())
                std::cout << "ERROR: stn[0]#" << tempi[0] << std::endl;
            if (tempi[1] < 0 | tempi[1] > returnNumStations ())
                std::cout << "ERROR: stn[1]#" << tempi[1] << std::endl;
            linetxt = linetxt.substr (ipos + 1, linetxt.length () - ipos - 1);
            d = atof (linetxt.c_str());
            dists (tempi [0], tempi [1]) = dists (tempi [1], tempi [0]) = d;
            count++;
        }
    }
    else 
    {
        /*
         * Distances between bike stations in boston, chicago, and washington,
         * and between train stations are read direct from matrix
         */
        while (getline (in_file, linetxt, '\n'))
            count++;
        if (count != _numStations)
        {
            std::cout << "ERROR: " << distFile << " does not have " <<
                _numStations << " columns and rows!" << std::endl;
            return -1;
        }
        else
        {
            in_file.clear ();
            in_file.seekg (0); 
            for (int i=0; i<_numStations; i++)
            {
                getline (in_file, linetxt, '\n');
                for (int j=0; j<(_numStations - 1); j++)
                {
                    ipos = linetxt.find (',',0);
                    dists (i, j) = dists (j, i) = 
                        atof (linetxt.substr (0, ipos).c_str());
                    linetxt = linetxt.substr (ipos + 1, 
                            linetxt.length () - ipos - 1);
                } // end for j over columns
                dists (i, _numStations - 1) = dists (_numStations - 1, i) =
                    atof (linetxt.c_str ());
            } // end for i over rows
        } // end else count == _numStations
    }

    in_file.close ();

    return 0;
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                             WRITEDMAT                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int StationData::writeDMat ()
{
    int numStations = StationList.size ();
    std::string distFile, nameFile;
    std::ofstream dists_out, names_out;

    if (_city == "london")
        distFile = "stationDistsMat_london.csv";
    else if (_city == "nyc")
        distFile = "stationDistsMat_nyc.csv";
    else if (_city == "oysterTube")
    {
        distFile = "DistMatTube.csv";
        nameFile = "DistMatTube_StationList.csv";
    }
    else if (_city == "oysterRail")
    {
        distFile = "DistMatRail.csv";
        nameFile = "DistMatRail_StationList.csv";
    }
    
    dists_out.open (distFile.c_str (), std::ofstream::out);
    names_out.open (nameFile.c_str (), std::ofstream::out);

    for (int i=0; i<numStations; i++)
    {
        names_out << StationList[i].name << ", " << StationList [i].lat <<
            ", " << StationList [i].lon << std::endl;
        for (int j=0; j<(numStations - 1); j++)
            dists_out << dists (i, j) << ", ";
        dists_out << dists (i, numStations - 1) << std::endl;
    }
    dists_out.close ();
    names_out.close ();

    std::cout << numStations << " inter-station distances written to " << 
        distFile.c_str () << std::endl;
    if (_city.substr (0, 6) == "oyster")
        std::cout << "\t... and corresponding station names written to " <<
            nameFile.c_str () << std::endl;

    return 0;
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                           WRITENUMTRIPS                            **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int StationData::writeNumTrips (std::string fname)
{
    int numStations = ntrips.size1 (); // < _numStations for rail data

    std::ofstream out_file;
    out_file.open (fname.c_str (), std::ofstream::out);
    for (int i=0; i<numStations; i++)
    {
        for (int j=0; j<numStations; j++)
        {
            out_file << ntrips (i, j);
            if (j == (numStations - 1))
                out_file << std::endl;
            else
                out_file << ", ";
        }
    }
    out_file.close ();
    std::cout << "Numbers of trips written to " << fname.c_str () << std::endl;

    return 0;
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                              CALCR2                                **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

// Also calculates covariance

int StationData::calcR2 (bool from)
{
    int tempi, numStations = StationList.size (), decLimsi [2];
    double tempd, decLims [2];
    std::vector <double> x0, y0, x1, x2, y2, d;
    RegrResults regrResults;

    d.resize (0);

    /* Double loop over numStations requires multiple (x,y) vectors
     * x0 is held in outer loop for i over (numStations - 1), potentially in
     * standardised (0-1) form.
     * y0 is held in inner loop for j over (i+1):numStations, also potentially
     * in standardised (0-1) form.
     *
     * For each of these inner loops, both vectors may be modified through
     * reduction to near/far stations only and/or removal of zeros. In each
     * case/both cases, the modified vectors are stored as (x2, y2), and then
     * copied after modification to (x1, y0), from which the correlations are
     * evaluated.
     *
     * x1 is necessary to allow it to revert to x0 in each inner loop.
     * The code ensures that all vectors have the same lengths at all times, and
     * so explicit loops rather than iterators are used for clarity.
     */

    for (int i=0; i<(numStations-1); i++)
    {
        x0.resize (0);
        for (int j=0; j<numStations; j++)
            if (from)
                x0.push_back (ntrips (i, j));
            else
                x0.push_back (ntrips (j, i));
        if (_standardise)
        {
            tempd = 0.0;
            for (int j=0; j<numStations; j++)
                tempd += x0 [j];
            for (int j=0; j<numStations; j++)
                x0 [j] = x0 [j] / tempd;
        }
        for (int j=(i+1); j<numStations; j++)
        {
            x1.resize (0);
            for (int k=0; k<numStations; k++)
                x1.push_back (x0 [k]);

            y0.resize (0);
            for (int k=0; k<numStations; k++)
                if (from)
                    y0.push_back (ntrips (j, k));
                else
                    y0.push_back (ntrips (k, j));

            if (_standardise)
            {
                tempd = 0.0;
                for (int k=0; k<numStations; k++)
                    tempd += y0 [k];
                for (int k=0; k<numStations; k++)
                    y0 [k] = y0 [k] / tempd;
            }

            if (nearfar > 0 && nearfar < 3) // Remove half of stations from lists
            {
                d.resize (0);
                for (int k=0; k<numStations; k++)
                    if (dists (i, k) > 0.0 && dists (j, k) > 0.0)
                        d.push_back (dists (i, k) + dists (j, k));
                if (d.size () > 10)
                {
                    std::sort (d.begin(), d.end());
                    tempi = floor (d.size () / 2);
                    tempd = (d [tempi] + d [tempi + 1]) / 2.0;
                    d.resize (0);
                    for (int k=0; k<numStations; k++)
                        d.push_back (dists (i, k) + dists (j, k));
                    x2.resize (0);
                    y2.resize (0);
                    for (int k=0; k<numStations; k++)
                    {
                        if (nearfar == 1 && d[k] < tempd)
                        {
                            x2.push_back (x1[k]);
                            y2.push_back (y0[k]);
                        } else if (nearfar == 2 && d[k] > tempd)
                        {
                            x2.push_back (x1[k]);
                            y2.push_back (y0[k]);
                        }
                    }
                    x1.resize (0);
                    y0.resize (0);
                    for (int k=0; k<x2.size(); k++)
                    {
                        x1.push_back (x2 [k]);
                        y0.push_back (y2 [k]);
                    }
                }
            } else if (nearfar >= 10) { // Analyse deciles
                d.resize (0);
                for (int k=0; k<numStations; k++)
                    if (dists (i, k) > 0.0 && dists (j, k) > 0.0)
                        d.push_back (dists (i, k) + dists (j, k));
                if (d.size () > 10)
                {
                    std::sort (d.begin(), d.end());
                    tempi = nearfar - 10;
                    decLimsi [0] = floor ((double) tempi * d.size () / 10.0);
                    decLimsi [1] = floor (((double) tempi + 1.0) * 
                                        d.size () / 10.0) - 1;
                    if (decLimsi [0] > 0)
                        decLims [0] = (d [decLimsi [0]] + d [decLimsi [0] + 1]) / 2.0;
                    if (tempi < 9)
                        decLims [1] = (d [decLimsi [1]] + d [decLimsi [1] + 1]) / 2.0;
                    
                    x2.resize (0);
                    y2.resize (0);
                    for (int k=0; k<numStations; k++)
                        if (d [k] > decLims [0] && d [k] < decLims [1])
                        {
                            x2.push_back (x1[k]);
                            y2.push_back (y0[k]);
                        }
                    x1.resize (0);
                    y0.resize (0);
                    for (int k=0; k<x2.size(); k++)
                    {
                        x1.push_back (x2 [k]);
                        y0.push_back (y2 [k]);
                    }
                }
            } // end if nearfar
            regrResults = regression (x1, y0);
            r2 (i, j) = r2 (j, i) = regrResults.r2;
            cov (i, j) = cov (j, i) = regrResults.cov;
            MI (i, j) = MI (j, i) = calcMI (x1, y0);
        } // end for j over (i+1):numStations
    } // end for i over (numStations - 1)
    x0.resize (0);
    x1.resize (0);
    x2.resize (0);
    y0.resize (0);
    y2.resize (0);

    return 0;
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                           GETDISTDECILES                           **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void StationData::writeDistDeciles ()
{
    int numStations = StationList.size ();
    double tempd, decLo, decHi;
    std::vector <int> counts1 (10), counts2 (10);
    std::vector <double> d, dtemp (10), decileDists (10);

    d.resize (0);

    for (int i=0; i<10; i++)
    {
        counts1 [i] = 0;
        decileDists [i] = 0.0;
    }

    // Some stations are in station lists, yet don't have lat-lons, and so can't
    // be mapped onto distance matrices. Thus distance matrices sometimes have
    // entire rows/columns of INT_MIN.
    for (int i=0; i<(numStations-1); i++)
    {
        for (int j=(i+1); j<numStations; j++)
        {
            d.resize (0);
            for (int k=0; k<numStations; k++)
                if (dists (i, k) > 0.0 && dists (j, k) > 0.0)
                    d.push_back (dists (i, k) + dists (j, k));
            if (d.size () > 10) 
            {
                std::sort (d.begin(), d.end());
                for (int k=0; k<10; k++)
                {
                    decLo = d [floor ((double) k * d.size () / 10.0)];
                    decHi = d [floor (((double) k + 1.0) * d.size () / 10.0) - 1];
                    counts2 [k] = 0;
                    dtemp [k] = 0.0;
                    for (int m=0; m<d.size (); m++)
                        if (d [m] > decLo && d [m] < decHi && d [m] > 0.0)
                        {
                            dtemp [k] += d [m] / 2.0;
                            counts2 [k] += 1;
                        }
                    decileDists [k] += dtemp [k] / (double) counts2 [k];
                    counts1 [k]++;
                }
            }
        } 
    }

    std::string fname = "DistDeciles_" + returnCity () + ".csv";
    std::ofstream out_file;
    out_file.open (fname.c_str (), std::ofstream::out);

    for (int i=0; i<10; i++)
        out_file << decileDists [i] / counts1 [i] << std::endl;
    out_file.close ();
    std::cout << "Distance deciles written to " << fname.c_str () << std::endl;

    counts1.resize (0);
    counts2.resize (0);
    dtemp.resize (0);
    d.resize (0);
    decileDists.resize (0);
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                             WRITER2MAT                             **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int StationData::writeR2Mat (std::string fname)
{
    int numStations = StationList.size ();

    std::ofstream out_file;
    out_file.open (fname.c_str (), std::ofstream::out);
    for (int i=0; i<numStations; i++) {
        for (int j=0; j<numStations; j++) {
            out_file << r2 (i, j);
            if (j == (numStations - 1))
                out_file << std::endl;
            else
                out_file << ", ";
        }
    }
    out_file.close ();
    std::cout << "Correlations written to " << fname.c_str () << std::endl;

    return 0;
}

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                            WRITECOVMAT                             **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int StationData::writeCovMat (std::string fname)
{
    int numStations = StationList.size ();

    std::ofstream out_file;
    out_file.open (fname.c_str (), std::ofstream::out);
    for (int i=0; i<numStations; i++) {
        for (int j=0; j<numStations; j++) {
            out_file << cov (i, j);
            if (j == (numStations - 1))
                out_file << std::endl;
            else
                out_file << ", ";
        }
    }
    out_file.close ();
    std::cout << "Covariances written to " << fname.c_str () << std::endl;

    return 0;
}



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                            WRITEMIMAT                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int StationData::writeMIMat (std::string fname)
{
    int numStations = StationList.size ();

    std::ofstream out_file;
    out_file.open (fname.c_str (), std::ofstream::out);
    for (int i=0; i<numStations; i++) {
        for (int j=0; j<numStations; j++) {
            out_file << MI (i, j);
            if (j == (numStations - 1))
                out_file << std::endl;
            else
                out_file << ", ";
        }
    }
    out_file.close ();
    std::cout << "Mutual Information written to " << fname.c_str () << std::endl;

    return 0;
}

