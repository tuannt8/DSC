//
//  debugger.h
//  DSC_compare
//
//  Created by Tuan Nguyen Trung on 12/02/2018.
//  Copyright © 2018 Asger Nyman Christiansen. All rights reserved.
//

#ifndef debugger_h
#define debugger_h
#include <map>
#include <iostream>
#include <string>

template <typename T = int>
class debugger
{
    std::map<std::string, T> m_debug_value;
public:
    debugger(){};
    ~debugger(){
    };
    
    void set(std::string key, T value){
        m_debug_value[key] = value;
    };
    T& get(std::string key, T default_v)
    {
        if (m_debug_value.find(key) == m_debug_value.end())
        {
            m_debug_value[key] = default_v;
        }
        
        return m_debug_value[key];
    }
    
    static int & get_int(std::string key, int default_v = 0)
    {
        static debugger<int> s_debug_i;
        return s_debug_i.get(key, default_v);
    }
};


#endif /* debugger_h */
