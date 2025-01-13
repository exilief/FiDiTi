#ifndef FIDITI_MATH_RECT_HPP
#define FIDITI_MATH_RECT_HPP

#include "Vec.hpp"


namespace fidi
{

template <int D, class T>
struct Rect
{
    Vec<D, T> min, max;

    Vec<D, T> size() const { return max - min; }

    Vec<D, T> corner(int i) const
    {
        assert(i < (1 << D));
        auto p = min, d = size();

        for (int j = 0; j < D; ++j)
        {
            // x: 01010101, y: 00110011, z: 00001111
            p[j] += d[j] * ((i / (1 << j)) % 2);
        }
        return p;
    }

    // Hypervolume (length, area, volume)
    T volume() const
    {
        T vol = 1;
        for (int i = 0; i < D; ++i)
            vol *= max[i] - min[i];
        return vol;
    }
};


template <int D, class T>
Vec<D, T> clamp(Vec<D, T> v, const Rect<D, T>& rect)
{
    for (int i = 0; i < D; ++i)
        v[i] = v[i] < rect.min[i] ? rect.min[i] :
               v[i] > rect.max[i] ? rect.max[i] : v[i];
    return v;
}

// Clamp low = clamp(v, Rect{vmin, +inf})
template <int D, class T>
Vec<D, T> clamp_l(Vec<D, T> v, const Vec<D, T>& vmin)
{
    for (int i = 0; i < D; ++i)
        if (v[i] < vmin[i])
            v[i] = vmin[i];
    return v;
}

// Clamp high = clamp(v, Rect{-inf, vmax})
template <int D, class T>
Vec<D, T> clamp_h(Vec<D, T> v, const Vec<D, T>& vmax)
{
    for (int i = 0; i < D; ++i)
        if (v[i] > vmax[i])
            v[i] = vmax[i];
    return v;
}

template <int D, class T>
Rect<D, T> clamp(const Rect<D, T>& r, const Rect<D, T>& bounds)
{
    return {clamp(r.min, bounds), clamp(r.max, bounds)};
}


template <int D, class T>
bool contains(const Rect<D, T>& r, const Vec<D, T>& v)
{
    for (int i = 0; i < D; ++i)
        if (v[i] < r.min[i] || v[i] > r.max[i])
            return false;
    return true;
}


template <int D, class F>
void forEachPoint(const Rect<D, int>& grid, F f)
{
    static_assert(D > 0);

    if constexpr(D == 1)
    {
        for (int i = grid.min.x; i <= grid.max.x; ++i)
            f(i);
    }
    else if constexpr(D == 2)
    {
        for (int j = grid.min.y; j <= grid.max.y; ++j)
            for (int i = grid.min.x; i <= grid.max.x; ++i)
                f(Vec2i{i, j});
    }
    else if constexpr(D == 3)
    {
        for (int k = grid.min.z; k <= grid.max.z; ++k)
            for (int j = grid.min.y; j <= grid.max.y; ++j)
                for (int i = grid.min.x; i <= grid.max.x; ++i)
                    f(Vec3i{i, j, k});
    }
    else
    {
        for (auto i = grid.min; i[D-1] <= grid.max[D-1]; )
        {
            f(i);

            int n;
            for (n = 0; n < D-1; ++n)
            {
                if (i[n] < grid.max[n])
                    break;
                i[n] = grid.min[n];
            }
            ++i[n];
        }
    }
}

template <int D, class F>
void forEachCell(const Rect<D, int>& grid, F f)  // forEachSubRect
{
    forEachPoint(Rect<D, int>{grid.min, grid.max - 1}, std::move(f));
}


template <class T> using Rect2 = Rect<2, T>;
template <class T> using Rect3 = Rect<3, T>;

using Rect2i = Rect2<int>;
using Rect2u = Rect2<unsigned>;
using Rect2f = Rect2<float>;
using Rect2d = Rect2<double>;

using Rect3i = Rect3<int>;
using Rect3u = Rect3<unsigned>;
using Rect3f = Rect3<float>;
using Rect3d = Rect3<double>;

template <int N> using RectNi = Rect<N, int>;
template <int N> using RectNu = Rect<N, unsigned>;
template <int N> using RectNf = Rect<N, float>;
template <int N> using RectNd = Rect<N, double>;

} // namespace fidi


#endif // FIDITI_MATH_RECT_HPP
