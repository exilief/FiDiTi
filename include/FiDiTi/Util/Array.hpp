#ifndef FIDITI_UTIL_ARRAY_HPP
#define FIDITI_UTIL_ARRAY_HPP


namespace fidi
{

// a[i] = a[0]
template <class T>
class RepeatArray
{
    T x = {};

 public:
    const T& operator[](int i) const { return x; }
    T&       operator[](int i) { return x; }
};


} // namespace fidi


#endif // FIDITI_UTIL_ARRAY_HPP
