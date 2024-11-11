#ifndef FIDITI_MATH_RECT_HPP
#define FIDITI_MATH_RECT_HPP

#include "Vec.hpp"


namespace fidi
{

template <size_t D, class T>
struct Rect
{
    Vec<D, T> min, max;

    Vec<D, T> size() const { return max - min; }

    Vec<D, T> corner(size_t i) const
    {
        assert(i < (1 << D));
        auto p = min, d = size();

        for (size_t j = 0; j < D; ++j)
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
        for (size_t i = 0; i < D; ++i)
            vol *= max[i] - min[i];
        return vol;
    }
};


template <size_t D, class T>
Vec<D, T> clamp(Vec<D, T> v, const Rect<D, T>& rect)
{
    for (size_t i = 0; i < D; ++i)
        v[i] = v[i] < rect.min[i] ? rect.min[i] :
               v[i] > rect.max[i] ? rect.max[i] : v[i];
    return v;
}

// Clamp low = clamp(v, Rect{vmin, +inf})
template <size_t D, class T>
Vec<D, T> clamp_l(Vec<D, T> v, const Vec<D, T>& vmin)
{
    for (size_t i = 0; i < D; ++i)
        if (v[i] < vmin[i])
            v[i] = vmin[i];
    return v;
}

// Clamp high = clamp(v, Rect{-inf, vmax})
template <size_t D, class T>
Vec<D, T> clamp_h(Vec<D, T> v, const Vec<D, T>& vmax)
{
    for (size_t i = 0; i < D; ++i)
        if (v[i] > vmax[i])
            v[i] = vmax[i];
    return v;
}

template <size_t D, class T>
Rect<D, T> clamp(const Rect<D, T>& r, const Rect<D, T>& bounds)
{
    return {clamp(r.min, bounds), clamp(r.max, bounds)};
}


template <size_t D, class F>
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

template <size_t D, class F>
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

template <size_t N> using RectNi = Rect<N, int>;
template <size_t N> using RectNu = Rect<N, unsigned>;
template <size_t N> using RectNf = Rect<N, float>;
template <size_t N> using RectNd = Rect<N, double>;

} // namespace fidi


#endif // FIDITI_MATH_RECT_HPP
