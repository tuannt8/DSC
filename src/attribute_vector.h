//
//  attibute_vector.h
//  3D_DSC
//
//  Created by Morten Nobel-Jørgensen.
//  Copyright (c) 2014 Morten Nobel-Jørgensen. All rights reserved.
//

#include <vector>
#include <cassert>

namespace DSC {

template<typename ITEMID, typename ITEM>
class AttributeVector
{
    std::vector<ITEM> items;
public:
    /// Construct from optional size and item
    AttributeVector(size_t size = 0, ITEM item = ITEM())
    {
        items = std::vector<ITEM>(size, item);
    }

    /// const reference to item given by ID
    const ITEM& get(ITEMID ID) const
    {
        assert(static_cast<unsigned int>(ID) < items.size());
        return items[ID];
    }

    /// reference to item given by ID
    ITEM& get(ITEMID ID)
    {
        if(ID >= items.size())
        {
            items.resize(ID + 1);
        }
        return items[ID];
    }

    /// const reference to item given by ID
    const ITEM& operator [](ITEMID ID) const
    {
        assert(static_cast<unsigned int>(ID) < items.size());
        return items[ID];
    }

    /// reference to item given by ID
    ITEM& operator [](ITEMID ID)
    {
        if(ID >= items.size())
        {
            items.resize(ID + 1);
        }
        return items[ID];
    }

    size_t size() const
    {
        return items.size();
    }

    void erase(std::vector<ITEMID> ids, ITEM item = ITEM()){
        for (auto id : ids){
            items[id] = item;
        }
    }
};
}