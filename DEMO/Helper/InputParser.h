//
//  InputParser.h
//  DEMO
//
//  Created by Tuan Nguyen Trung on 01/02/2018.
//

// Parser input arguments of the main function

#ifndef InputParser_h
#define InputParser_h

#include <vector>
#include <string>

class InputParser{
public:
    InputParser (int &argc, char **argv){
        for (int i=1; i < argc; ++i)
            this->tokens.push_back(std::string(argv[i]));
    }
    /// @author iain
    const std::string& getCmdOption(const std::string &option) const{
        std::vector<std::string>::const_iterator itr;
        itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
        if (itr != this->tokens.end() && ++itr != this->tokens.end()){
            return *itr;
        }
        static const std::string empty_string("");
        return empty_string;
    }
    
    std::string getCmdOption(const std::string &option, const std::string &default_){
        std::vector<std::string>::const_iterator itr;
        itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
        if (itr != this->tokens.end() && ++itr != this->tokens.end()){
            return *itr;
        }

        return default_;
    }
    /// @author iain
    bool cmdOptionExists(const std::string &option) const{
        return std::find(this->tokens.begin(), this->tokens.end(), option)
        != this->tokens.end();
    }
private:
    std::vector <std::string> tokens;
};

#endif /* InputParser_h */
