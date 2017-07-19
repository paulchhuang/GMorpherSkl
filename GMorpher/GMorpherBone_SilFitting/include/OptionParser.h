/* *************************************************
 * Copyright (2011) : Cedric Cagniart
 * *************************************************/
#ifndef OPTIONPARSER_H_DEFINED
#define OPTIONPARSER_H_DEFINED

#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <cstdio>
#include <stdexcept>
#include <iostream>




inline std::string buildFilename(const std::string &baseName, int id)
{
	char cfilename[1024]; //buffer for the filename
	sprintf_s(cfilename, baseName.c_str(), id);
	return std::string(cfilename);
}

inline std::string buildFilename(const std::string &baseName, int id1, int id2)
{
	char cfilename[1024]; //buffer for the filename
	sprintf_s(cfilename, baseName.c_str(), id1, id2);
	return std::string(cfilename);
}

inline std::string buildFilename(const std::string &baseName, const std::string& s, int id)
{
	char cfilename[1024]; //buffer for the filename
	sprintf_s(cfilename, baseName.c_str(), s.c_str(), id);
	return std::string(cfilename);
}





/// current restriction, flags have to be separated from the argument
class OptionParser
{
    public :
    inline OptionParser(int argc, char**argv);

    template <class T>
    const T getOption(const char *optName, bool verbose = true);

    protected :
    std::vector<std::string > mArgV;
};



inline OptionParser::OptionParser(int argc, char**argv)
{
    for(int i=0;i<argc;++i) mArgV.push_back(std::string(argv[i]));
}


template <class T>
const T OptionParser::getOption(const char *optName, bool verbose )
{
    T r;
    for(unsigned int i=0;i<mArgV.size()-1;++i)
    {
        if(mArgV[i].compare(optName) == 0)
        {
            std::stringstream ss(mArgV[i+1]);
            ss >> r;
            if(verbose) std::cout<<"OPT "<<optName<<" : "<<r<<std::endl;
            return r;
        }
    }
    throw std::runtime_error( std::string("could not find option in the command line ") + std::string(optName) );
}



#endif
