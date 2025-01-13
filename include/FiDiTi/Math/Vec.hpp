#ifndef FIDITI_MATH_VEC_HPP
#define FIDITI_MATH_VEC_HPP

#include "../config.hpp"

#include <cassert>
#include <utility>


namespace fidi
{

template <int N, class T>
struct Vec
{
    T v[N] = {};

    constexpr Vec() = default;

    // Fill Vec with x
    constexpr explicit Vec(T x);

    //Vec(initializer_list<T>)

    template <class U>
    constexpr explicit Vec(const Vec<N,U>& vec);

    constexpr T& operator[](int i) { return v[i]; }

    constexpr const T& operator[](int i) const { return v[i]; }

    //static Vec filled(T x) { return Vec(x); }
    //static Vec basis(int i, T x=1, T fill={}) { Vec v = filled(fill);  v[i] = x;  return v; }
};

template <int N, class T>
Vec<N, T> basisVec(int i, T x = 1, T fill = {});

template <int N, int N0, class T>
Vec<N, T> resize(const Vec<N0, T>& v);


// -1, +1 (includes 0)
template <int D, class T>
Vec<D, int> sign(const Vec<D, T>& v);

// -1, 0, +1
template <int D, class T>
Vec<D, int> signum(const Vec<D, T>& v);


template <int D, class T>
T dot(const Vec<D, T>&, const Vec<D, T>&);

// Project onto coordinate plane with axis i (set component i to zero)
template <int D, class T>
Vec<D, T> project(Vec<D, T> v, int i);


/**
// =========== //
// [Operators] //
// =========== //

// Unary op: + -

Vec<N,T> operator +(const Vec<N,T>& right);


// Binary op: + - * /

// Vector op= Vector
Vec<N,T>& operator +=(Vec<N,T>& left, const Vec<N,T>& right);

// Vector op= Scalar
Vec<N,T>& operator +=(Vec<N,T>& left, T right);

// Vector op Vector
Vec<N,T> operator +(Vec<N,T> left, const Vec<N,T>& right);

// Vector op Scalar
Vec<N,T> operator +(Vec<N,T> left, T right);

// Scalar op Vector
Vec<N,T> operator +(T left, Vec<N,T> right);

**/

template <int N, class T>
constexpr bool operator ==(const Vec<N,T>& left, const Vec<N,T>& right);

template <int N, class T>
constexpr bool operator !=(const Vec<N,T>& left, const Vec<N,T>& right);

#include "Vec.inl"





// Specializations to allow v.x, v.y, v.z

template <class T>
struct Vec<1, T>
{
    T x = {};

    constexpr Vec() = default;

    constexpr Vec(T X) : x(X) { }

    template <class U>
    constexpr explicit Vec(const Vec<1,U>& vec) : x(T(vec.x)) { }

    constexpr const T& operator[](int i) const
    {
        assert(i == 0);
        return x;
    }

    constexpr T& operator[](int i)
    {
        assert(i == 0);
        return x;
    }
};


template <class T>
struct Vec<2, T>
{
    T x = {}, y = {};

    constexpr Vec() = default;

    // Fill Vec with X
    constexpr explicit Vec(T X) : x(X), y(X) { }

    constexpr Vec(T X, T Y) : x(X), y(Y) { }

    template <class U>
    constexpr explicit Vec(const Vec<2,U>& vec) : x(T(vec.x)), y(T(vec.y)) { }

    constexpr const T& operator[](int i) const
    {
        assert(i < 2);

        switch(i)
        {
          case 0:  return x;
          default: return y;
        }
    }

    constexpr T& operator[](int i)
    {
        return const_cast<T&>(std::as_const(*this)[i]);
    }
};


template <class T>
struct Vec<3, T>
{
    T x = {}, y = {}, z = {};

    constexpr Vec() = default;

    // Fill Vec with X
    constexpr explicit Vec(T X) : x(X), y(X), z(X) { }

    constexpr Vec(T X, T Y, T Z) : x(X), y(Y), z(Z) { }

    template <class U>
    constexpr explicit Vec(const Vec<3,U>& vec) : x(T(vec.x)), y(T(vec.y)), z(T(vec.z)) { }

    constexpr const T& operator[](int i) const
    {
        assert(i < 3);

        switch(i)
        {
          case 0:  return x;
          case 1:  return y;
          default: return z;
        }
    }

    constexpr T& operator[](int i)
    {
        return const_cast<T&>(std::as_const(*this)[i]);
    }
};


template <class T> using Vec2 = Vec<2, T>;
template <class T> using Vec3 = Vec<3, T>;

using Vec2i = Vec2<int>;
using Vec2u = Vec2<unsigned>;
using Vec2f = Vec2<float>;
using Vec2d = Vec2<double>;

using Vec3i = Vec3<int>;
using Vec3u = Vec3<unsigned>;
using Vec3f = Vec3<float>;
using Vec3d = Vec3<double>;

template <int N> using VecNi = Vec<N, int>;
template <int N> using VecNu = Vec<N, unsigned>;
template <int N> using VecNf = Vec<N, float>;
template <int N> using VecNd = Vec<N, double>;

} // namespace fidi


#endif // FIDITI_MATH_VEC_HPP
