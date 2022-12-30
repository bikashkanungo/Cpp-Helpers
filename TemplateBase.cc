template <typename T>
Base<T>::Base()
{

  std::cout << "Base Constructor" << std::endl;

}

template <typename T>
Base<T>::~Base()
{

}

template <typename T>
void
Base<T>::setValue(T val)
{

  d_val = val;

}

template <typename T>
T
Base<T>::getValue()
{

  return d_val;

}

