//
//  config_file.h
//  DSC_compare
//
//  Created by Tuan Nguyen Trung on 01/02/2018.
//  Copyright © 2018 Asger Nyman Christiansen. All rights reserved.
//

#ifndef config_file_h
#define config_file_h

#include <string>
#include <map>
#include <iomanip>
#include <iostream>
#include <fstream>

// Load config file
// Format:
// key = value
// Any line without the '=' token will be considered a comment
// No space allowed in the key as well as value
// support value type: int, string, double, bool (yes/no)

class config_file
{
public:
    config_file(std::string path){
        std::ifstream infile(path);
        
        for (std::string line; std::getline(infile, line); )
        {
            line.erase(std::remove(line.begin(),line.end(),' '),line.end());
            if (line.empty())
            {continue;}
            
            size_t pos = line.find_first_of("=");
            if(pos > line.size() -1)
                continue; // This is a comment line
            
            std::string key = line.substr(0, pos);
            std::string val = line.substr(pos+1, line.size()-1);
            
            options[key] = val;
        }
    };
    ~config_file(){};
    
    int get_int(const std::string key, int default_=0)
    {
        if (options.find(key) == options.end())
        {
            return default_;
        }
        else
        {
            return atoi(options[key].c_str());
        }
    };
    
    float get_double(const std::string key, float default_=0)
    {
        if (options.find(key) == options.end())
        {
            return default_;
        }
        else
        {
            return atof(options[key].c_str());
        }
    };
    
    std::string get_string(const std::string key, std::string default_="")
    {
        if (options.find(key) == options.end())
        {
            return default_;
        }
        else
        {
            return options[key];
        }
    }
    
private:
    std::map<std::string, std::string> options;
};

#endif /* config_file_h */
