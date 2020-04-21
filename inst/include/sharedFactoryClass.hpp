#include "sharedFactoryClass.h"

template <class BaseObjectType>
typename SharedFactory<BaseObjectType>::SharedPointerType
SharedFactory<BaseObjectType>::Instantiate(const std::string &name)
{
    auto it = m_RegistryMap.find(name);
    return it == m_RegistryMap.end() ? nullptr : (it->second)();
}

template <class BaseObjectType>
template <class DerivedObjectType>
void SharedFactory<BaseObjectType>::Register(const std::string &name)
{
  m_RegistryMap[name] = []()
  {
    return std::make_shared<DerivedObjectType>();
  };
}
