//
//  cache.hpp
//  interface_tracking
//
//  Created by Tuan Nguyen Trung on 11/7/16.
//  Copyright © 2016 Tuan Nguyen Trung. All rights reserved.
//

#ifndef cache_hpp
#define cache_hpp

#include <stdio.h>
#include <vector>
#include <map>
#include <iostream>
#include <queue>



class cache_item_base
{
public:
    cache_item_base(){};
    ~cache_item_base(){};
    
    virtual void release(){};
    virtual void * get(int idx){return nullptr;};
    virtual void set(void * newdata, int idx){};
    virtual void mark_dirty(int idx){};
};

template<typename type>
class dynamic_vector : public std::vector<type>
{
public:
    type & operator [](const size_t index)
    {
        if(index > this->size())
            this->resize(index+1);
            
        return this->at(index);
    }
};

template<typename type>
class cache_item : public cache_item_base
{
    dynamic_vector<type *> _data;
    
public:
    cache_item(){

    }
    
    ~cache_item(){
        
    }
    
    virtual void release(){
        for(auto & d : _data){
            if(d){
                delete d;
                d = nullptr;
            }
        }
    };
    
    virtual void * get(int idx){
        return _data[idx];
    }
    
    virtual void set(void * newdata, int idx){
        if(_data[idx])
        {
            delete _data[idx];
            _data[idx] = nullptr;
        }
        _data[idx] = (type*)newdata;
    };
    
    virtual void mark_dirty(int idx){
        if(_data[idx]){
            delete _data[idx];
            _data[idx] = nullptr;
        }
    }
};

class cache
{
public:
    template<typename type>
    type* get_cache(int cache_id, int item_id, std::function<type*(int)>compute)
    {
        // create cache buffer if not exist.
        // Does it take time???
        if (_data.find(cache_id) == _data.end())
        {
            add_cache_type<type>(cache_id, compute);
        }
        
        // If cache does not exist, build
        if(!(_data[cache_id])->get(item_id))
        {
            _data[cache_id]->set((void*)compute(item_id), item_id);
        }
        
        // Return result
        return (type*)(_data[cache_id])->get(item_id);
    }
    
    void mark_dirty(int cache_id, int item_id)
    {
        _data[cache_id]->mark_dirty(item_id);
    }
    
    static cache * get_instance(){
        static cache instance;
        return &instance;
    }
    
private:
    template<typename type>
    bool add_cache_type(int id, std::function<type*(int)>compute){
        cache_item<type> * new_array = new cache_item<type>;
        _data.insert(std::make_pair(id, (cache_item_base*)new_array));
        
        //Sample time and memory
        
        return true;
    }
    
private:
    std::map<int, cache_item_base*> _data;
    std::queue<int> _cache_track; // For first in first out cache
    
    cache(){}
    ~cache(){
        for(auto it = _data.begin(); it != _data.end(); it++)
        {
            it->second->release();
            delete it->second;
        }
        _data.clear();
    }
};


#endif /* cache_hpp */
