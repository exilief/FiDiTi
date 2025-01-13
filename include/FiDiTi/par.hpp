#ifndef FIDITI_PAR_HPP
#define FIDITI_PAR_HPP

#include "Math/Rect.hpp"

#ifdef _OPENMP
  #include <omp.h>
#endif


namespace fidi::par
{

int get_thread_num()
{
  #ifdef _OPENMP
    return omp_get_thread_num();
  #endif
    return 0;
}

int get_num_threads()
{
  #ifdef _OPENMP
    return omp_get_num_threads();
  #endif
    return 1;
}

int get_max_threads()
{
  #ifdef _OPENMP
    return omp_get_max_threads();
  #endif
    return 1;
}

void set_num_threads(int N)
{
  #ifdef _OPENMP
    omp_set_num_threads(N);
  #endif
}

void set_num_threads_max()
{
    set_num_threads(get_max_threads());
}


template <int D, class F>
void forEachPoint(const Rect<D, int>& grid, F f)
{
    static_assert(D > 0);

    if constexpr(D == 1)
    {
        #pragma omp parallel for
        for (int i = grid.min.x; i <= grid.max.x; ++i)
            f(i);
    }
    else if constexpr(D == 2)
    {
        #pragma omp parallel for
        for (int j = grid.min.y; j <= grid.max.y; ++j)
            for (int i = grid.min.x; i <= grid.max.x; ++i)
                f(Vec2i{i, j});
    }
    else if constexpr(D == 3)
    {
        #pragma omp parallel for
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
void forEachCell(const Rect<D, int>& grid, F f)
{
    par::forEachPoint(Rect<D, int>{grid.min, grid.max - 1}, std::move(f));
}


} // namespace fidi::par


#endif // FIDITI_PAR_HPP
