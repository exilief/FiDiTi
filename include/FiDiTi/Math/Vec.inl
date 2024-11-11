
template <size_t N, class T>
constexpr Vec<N,T>::Vec(T x)
{ 
    for (auto& xi : v)
        xi = x;
}


template <size_t N, class T>
template <class U>
constexpr Vec<N,T>::Vec(const Vec<N,U>& vec)
{
    for (size_t i = 0; i < N; ++i)
        v[i] = T(vec[i]);
}


template <size_t N, class T>
Vec<N,T> basisVec(size_t i, T x, T fill)
{
    Vec<N,T> v(fill);
    v[i] = x;
    return v;
}

template <size_t N, size_t N0, class T>
Vec<N, T> resize(const Vec<N0, T>& v)
{
    if constexpr(N == N0)
        return v;

    Vec<N, T> u;
    for (size_t i = 0; i < (N < N0 ? N : N0); ++i)
        u[i] = v[i];
    return u;
}



template <size_t D, class T>
Vec<D, int> sign(const Vec<D, T>& v)
{
    Vec<D, int> s;
    for (int i = 0; i < D; ++i)
        s[i] = v[i] < 0 ? -1 : 1;
    return s;
}

template <size_t D, class T>
Vec<D, int> signum(const Vec<D, T>& v)
{
    Vec<D, int> s;
    for (int i = 0; i < D; ++i)
        s[i] = v[i] < 0 ? -1 : (v[i] > 0 ? 1 : 0);
    return s;
}


template <size_t D, class T>
T dot(const Vec<D, T>& v1, const Vec<D, T>& v2)
{
    T sum = 0;
    for (int i = 0; i < D; ++i)
        sum += v1[i] * v2[i];
    return sum;
}

template <size_t D, class T>
Vec<D, T> project(Vec<D, T> v, size_t i)
{
    assert(i < D);
    v[i] = 0;
    return v;
}





#define VEC_DEF_OP_UNARY(op)                    \
    template <size_t N, class T>           \
    constexpr Vec<N,T> operator op(Vec<N,T> v)  \
    {                                           \
        for (size_t i=0; i < N; ++i)        \
            v[i] = op v[i];                     \
        return v;                               \
    }

#define VEC_DEF_OP_BINARY(op, op_eq)                                            \
    template <size_t N, class T>                                           \
    constexpr Vec<N,T>& operator op_eq(Vec<N,T>& left, const Vec<N,T>& right)   \
    {                                                                           \
        for (size_t i=0; i < N; ++i)                                        \
            left[i] op_eq right[i];                                             \
        return left;                                                            \
    }                                                                           \
    template <size_t N, class T>                                           \
    constexpr Vec<N,T>& operator op_eq(Vec<N,T>& left, T right)                 \
    {                                                                           \
        for (size_t i=0; i < N; ++i)                                        \
            left[i] op_eq right;                                                \
        return left;                                                            \
    }                                                                           \
\
    template <size_t N, class T>                                       \
    constexpr Vec<N,T> operator op(Vec<N,T> left, const Vec<N,T>& right)    \
    {                                                                       \
        return left op_eq right;                                            \
    }                                                                       \
    template <size_t N, class T>                                       \
    constexpr Vec<N,T> operator op(Vec<N,T> left, T right)                  \
    {                                                                       \
        return left op_eq right;                                            \
    }                                                                       \
    template <size_t N, class T>                                       \
    constexpr Vec<N,T> operator op(T left, Vec<N,T> right)                  \
    {                                                                       \
        for (size_t i=0; i < N; ++i)                                    \
            right[i] = left op right[i];                                    \
        return right;                                                       \
    }

VEC_DEF_OP_UNARY(+)
VEC_DEF_OP_UNARY(-)

VEC_DEF_OP_BINARY(+, +=)
VEC_DEF_OP_BINARY(-, -=)
VEC_DEF_OP_BINARY(*, *=)
VEC_DEF_OP_BINARY(/, /=)

#undef VEC_DEF_OP_UNARY
#undef VEC_DEF_OP_BINARY




template <size_t N, class T>
constexpr bool operator ==(const Vec<N,T>& left, const Vec<N,T>& right)
{
    for (size_t i=0; i < N; ++i)
        if(left[i] != right[i])
            return false;
    return true;
}


template <size_t N, class T>
constexpr bool operator !=(const Vec<N,T>& left, const Vec<N,T>& right)
{
    return !(left == right);
}
