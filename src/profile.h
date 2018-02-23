//
//  profile.hpp
//  DSC
//
//  Created by Tuan Nguyen Trung on 2/4/16.
//  Copyright © 2016 Asger Nyman Christiansen. All rights reserved.
//

#ifndef profile_hpp
#define profile_hpp

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <chrono>
#include <map>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>

#define P_TIME_NOW (std::chrono::system_clock::now())
typedef std::chrono::duration<double> p_duration_t;
typedef std::chrono::system_clock::time_point p_time_point;

struct profile_att
{
    int count = 0;
    double total_time = 0;
    p_time_point m_start;
};

class profile_temp
{
    std::string m_message;
    p_time_point m_start;
public:
    profile_temp(std::string m = ""){
        m_message = m;
        m_start = P_TIME_NOW;
    }
    ~profile_temp(){
        p_duration_t t = P_TIME_NOW - m_start;
        std::cout << "Time " << m_message << ": " << t.count() << std::endl;
    }
};

class profile
{
public:
    static void init();
    static void close();
    
    profile(std::string name);
    
    void change(std::string name);
    void done();
    
    profile();
    ~profile();
private:

    static std::map<std::string, profile_att> m_objects;
    static profile_att * get_object(const std::string &  name );
private:
    
    std::string m_name;
};

class profile_progress
{
public:
    profile_progress(std::string name){
        _file_name = name;
        _last_point = P_TIME_NOW;
    }
    ~profile_progress(){
        std::ofstream f(std::string("LOG/") + _file_name + std::string(".txt"));
        if (f.is_open())
        {
            for (auto t : _time_iter)
            {
                f << t << std::endl;
            }
            
            f.close();
        }
    }
    
    void add()
    {
        p_duration_t t = P_TIME_NOW - _last_point;
        _time_iter.push_back(t.count());
        _last_point = P_TIME_NOW;
    }
    
private:
    std::vector<double> _time_iter;// By iteration
    std::string _file_name;
    
    p_time_point _last_point;
};

#endif /* profile_hpp */
