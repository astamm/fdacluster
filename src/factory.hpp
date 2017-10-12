#ifndef FACTORY_HPP_
#define FACTORY_HPP_

#include<RcppArmadillo.h>
#include<memory>

using namespace arma;


// Factory class
template<typename D>
class SharedFactory
{

public:
    typedef std::unordered_map< std::string, std::function< std::shared_ptr<D>() > > registry_map;

    registry_map map;

    // use this to instantiate the proper Derived class
    std::shared_ptr<D> instantiate(const std::string& name)
    {
        auto it = map.find(name);
        return it == map.end() ? nullptr : (it->second)();
    }

    template<typename T>
    void FactoryRegister(std::string name)
    {
        map[name] = []()
        {
            return std::make_shared<T>();
        };
        //std::cout << "Registering class '" << name << "'\n";
    }

};

#endif
