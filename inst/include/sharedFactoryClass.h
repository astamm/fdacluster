#ifndef SHAREDFACTORYCLASS_H
#define SHAREDFACTORYCLASS_H

#include <functional> // for std::function
#include <memory> // for std::shared_ptr...
#include <string>
#include <unordered_map>

/// Factory class
template <class BaseObjectType>
class SharedFactory
{
public:
  using SharedPointerType = std::shared_ptr<BaseObjectType>;
  using RegistryMapType = std::unordered_map<std::string, std::function<SharedPointerType()> >;

  // Use this to instantiate the proper Derived class
  SharedPointerType Instantiate(const std::string &name);

  template <class DerivedObjectType> void Register(const std::string &name);

private:
  RegistryMapType m_RegistryMap;
};

#include "sharedFactoryClass.hpp"

#endif /* SHAREDFACTORYCLASS_H */
