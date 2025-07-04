
template <int N, class T>
constexpr Vec<N,T>::Vec(T x)
{ 
    for (auto& xi : v)
        xi = x;
}


template <int N, class T>
template <class U>
constexpr Vec<N,T>::Vec(const Vec<N,U>& vec)
{
    for (int i = 0; i < N; ++i)
        v[i] = T(vec[i]);
}


template <int N, class T>
Vec<N,T> basisVec(int i, T x, T fill)
{
    assert(i < N);

    Vec<N,T> v(fill);
    v[i] = x;
    return v;
}

template <int N, class T>
Vec<N,T> basisVec0(int i, T x, T fill)
{
    Vec<N,T> v(fill);
    if (i < N)
        v[i] = x;
    return v;
}

template <int N, int N0, class T>
Vec<N, T> resize(const Vec<N0, T>& v)
{
    if constexpr(N == N0)
        return v;

    Vec<N, T> u;
    for (int i = 0; i < (N < N0 ? N : N0); ++i)
        u[i] = v[i];
    return u;
}

template <int N, class T>
Vec<N, T> toVec(const T* data)
{
    Vec<N, T> v;

    for (int i = 0; i < N; ++i)
        v[i] = data[i];

    return v;
}



template <int D, class T>
Vec<D, int> sign(const Vec<D, T>& v)
{
    Vec<D, int> s;
    for (int i = 0; i < D; ++i)
        s[i] = v[i] < 0 ? -1 : 1;
    return s;
}

template <int D, class T>
Vec<D, int> signum(const Vec<D, T>& v)
{
    Vec<D, int> s;
    for (int i = 0; i < D; ++i)
        s[i] = v[i] < 0 ? -1 : (v[i] > 0 ? 1 : 0);
    return s;
}


template <int D, class T>
T len_sq(const Vec<D, T>& v)
{
    return dot(v, v);
}

template <int D, class T>
Vec<D, T> project(Vec<D, T> v, int i)
{
    assert(i < D);
    v[i] = 0;
    return v;
}

template <int D, class T>
T dot(const Vec<D, T>& v1, const Vec<D, T>& v2)
{
    T sum = 0;
    for (int i = 0; i < D; ++i)
        sum += v1[i] * v2[i];
    return sum;
}

template <class T>
Vec<3, T> cross(const Vec<3, T>& u, const Vec<3, T>& v)
{
    Vec<3, T> w;
    for (int i = 0; i < 3; ++i)
    {
        int j = (i+1) % 3;
        int k = (i+2) % 3;
        w[i] += u[j] * v[k] - u[k] * v[j];
    }
    return w;
}





#define VEC_DEF_OP_UNARY(op)                    \
    template <int N, class T>                   \
    constexpr Vec<N,T> operator op(Vec<N,T> v)  \
    {                                           \
        for (int i=0; i < N; ++i)               \
            v[i] = op v[i];                     \
        return v;                               \
    }

#define VEC_DEF_OP_BINARY(op, op_eq)                                            \
    template <int N, class T>                                                   \
    constexpr Vec<N,T>& operator op_eq(Vec<N,T>& left, const Vec<N,T>& right)   \
    {                                                                           \
        for (int i=0; i < N; ++i)                                               \
            left[i] op_eq right[i];                                             \
        return left;                                                            \
    }                                                                           \
    template <int N, class T, class U>                                          \
    constexpr Vec<N,T>& operator op_eq(Vec<N,T>& left, U right)                 \
    {                                                                           \
        for (int i=0; i < N; ++i)                                               \
            left[i] op_eq right;                                                \
        return left;                                                            \
    }                                                                           \
\
    template <int N, class T>                                               \
    constexpr Vec<N,T> operator op(Vec<N,T> left, const Vec<N,T>& right)    \
    {                                                                       \
        return left op_eq right;                                            \
    }                                                                       \
    template <int N, class T, class U>                                      \
    constexpr Vec<N,T> operator op(Vec<N,T> left, U right)                  \
    {                                                                       \
        return left op_eq right;                                            \
    }                                                                       \
    template <int N, class T, class U>                                      \
    constexpr Vec<N,T> operator op(U left, Vec<N,T> right)                  \
    {                                                                       \
        for (int i=0; i < N; ++i)                                           \
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




template <int N, class T>
constexpr bool operator ==(const Vec<N,T>& left, const Vec<N,T>& right)
{
    for (int i=0; i < N; ++i)
        if(left[i] != right[i])
            return false;
    return true;
}


template <int N, class T>
constexpr bool operator !=(const Vec<N,T>& left, const Vec<N,T>& right)
{
    return !(left == right);
}
