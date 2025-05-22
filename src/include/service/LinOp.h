/*<<???????? ???????? ?????? (????????? ? ??????)>>
 @author: ??????? ????
 2021-2024*/
#pragma once
#include<vector>
#include<iomanip>
#include<iostream>
using namespace std;

template <typename T> inline constexpr
int signum(T x, std::false_type is_signed) {
    return T(0) < x;
}

template <typename T> inline constexpr
int signum(T x, std::true_type is_signed) {
    return (T(0) < x) - (x < T(0));
}

template <typename T> inline constexpr
int signum(T x) {
    return signum(x, std::is_signed<T>());
}

// ???????? ??????
template<typename LT>
vector<vector<LT>> mat_sum(const vector<vector<LT>>& matr1, const vector<vector<LT>>& matr2)
{
    vector<vector <LT>> product(matr1);
    for (int i = 0; i < matr1.size(); ++i)
        for (int j = 0; j < matr2.size(); ++j)
            product[i][j] += matr2[i][j];
    return product;
}
template<typename LT>
vector<vector<LT>> operator + (const vector<vector<LT>> matr1, const vector<vector<LT>> matr2)
{
    return mat_sum(matr1, matr2);
}
//????????? ??????? ?? ???????
template<typename LT>
vector<vector<LT>> mat_mult(const vector<vector<LT>> matr1, const vector<vector<LT>> matr2)
{
    vector<vector<LT>> product;
    int n = matr1[0].size();
    product.push_back({});
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            LT temp_sum = 0;
            for (int k = 0; k < n; ++k)
                temp_sum += matr1[i][k] * matr2[k][j];
            product[i].push_back(temp_sum);
        }
        product.push_back({});
    }
    return product;
}

template<typename LT>
vector<vector<LT>> operator * (const vector<vector<LT>> matr1, const vector<vector<LT>> matr2)
{
    return mat_mult(matr1, matr2);
}

//????????? ??????? ?? ?????
template<typename LT, typename LT2>
vector<vector<LT>> operator * (const vector<vector<LT>>& matr, const LT2 lambda)
{
    vector<vector<LT>> result(matr);
    for (int i = 0; i < matr.size(); ++i)
        for (int j = 0; j < matr[0].size(); ++j)
            result[i][j] *= lambda;
    return result;
}
template<typename LT, typename LT2>
vector<vector<LT>> operator * (const LT2 lambda, const vector<vector<LT>>& matr)
{
    return matr * lambda;
}
//????????? ??????? ?? ??????
template<typename LT>
vector<LT> operator * (const vector<vector<LT>>& matr, const vector<LT> vec)
{
    vector<LT> result;
    for (int i = 0; i < vec.size(); ++i)
    {
        LT tmp_val = 0;
        for (int j = 0; j < vec.size(); ++j)
        {
            tmp_val += matr[i][j] * vec[j];
        }
        result.push_back(tmp_val);
    }
    return result;
}
//????????? ??????? ?? ?????
template<typename LT, typename LT2>
vector<LT> operator * (const vector<LT>& vec, const LT2& lambda)
{
    vector<LT> result;
    for (auto el : vec)
    {
        result.push_back(el * lambda);
    }
    return result;
}
template<typename LT, typename LT2>
vector<LT> operator * (const LT2& lambda, const vector<LT>& vec)
{
    return vec * lambda;
}

//????????? ????????????
template <typename LT>
LT scmult(vector<LT> v1, vector<LT> v2)
{
    LT result = 0;
    for (auto v = v1.begin(), p = v2.begin(); v != v1.end() && p != v2.end(); v++, p++)
    {
        result += *v * *p;
    }
    return result;
}
template<typename LT>
LT operator * (const vector<LT>& v1, const vector<LT>& v2)
{
    return scmult(v1, v2);
}
//?????
template <typename LT>
LT vec_norm(vector<LT> v1)
{
    return sqrt(scmult<LT>(v1, v1));
}
//c????
template <typename LT>
vector<LT> vsum(vector<LT>& v1, vector<LT>& v2)
{
    vector<LT> result;
    auto v = v1.begin();
    auto p = v2.begin();
    for (; v != v1.end() && p != v2.end(); v++, p++)
    {
        result.push_back(*v + *p);
    }
    for (; v != v1.end(); v++)
    {
        result.push_back(*v);
    }
    for (; p != v2.end(); p++)
    {
        result.push_back(*p);
    }
    return result;
}
//????????
template <typename LT>
vector<LT> vdiff(vector<LT>& v1, vector<LT>& v2)
{
    vector<LT> result;
    auto v = v1.begin();
    auto p = v2.begin();
    for (; v != v1.end() && p != v2.end(); v++, p++)
    {
        result.push_back(*v - *p);
    }
    for (; v != v1.end(); v++)
    {
        result.push_back(*v);
    }
    for (; p != v2.end(); p++)
    {
        result.push_back(-*p);
    }
    return result;
}

//???????? ????????
template<typename LT>
vector<LT> operator + (vector<LT> r_vec, vector<LT> l_vec)
{
    return vsum(r_vec, l_vec);
}
//????????? ????????
template<typename LT>
vector<LT> operator - (vector<LT> r_vec, vector<LT> l_vec)
{
    return vdiff(r_vec, l_vec);
}

//??????????? ?????????*lambda  ? ???????
template<typename LT, typename LT2>
vector<vector<LT>> operator + (vector<vector<LT>> matr, LT2 lambda)
{
    vector<vector<LT>> result(matr);
    for (int i = 0; i < matr[0].size(); ++i)
        result[i][i] += lambda;
    return result;
}
template<typename LT, typename LT2>
vector<vector<LT>> operator + (LT2 lambda, vector<vector<LT>> matr)
{
    return matr + lambda;
}
//????? ??????
template<typename LT>
LT norm1(vector<vector<LT>> A)
{
    int n = A[0].size();
    LT max_val = 0;
    for (int j = 0; j < n; ++j)
    {
        LT temp_sum = 0;
        for (int i = 0; i < n; ++i)
        {
            temp_sum += fabs(A[i][j]);
        }
        if (temp_sum > max_val)
            max_val = temp_sum;
    }
    return max_val;
}

template<typename LT>
LT norm_inf(vector<vector<LT>> A)
{
    int n = A[0].size();
    LT max_val = 0;
    for (int i = 0; i < n; ++i)
    {
        LT temp_sum = 0;
        for (int j = 0; j < n; ++j)
        {
            temp_sum += fabs(A[i][j]);
        }
        if (temp_sum > max_val)
            max_val = temp_sum;
    }
    return max_val;
}
template<typename LT>
LT vec_norm_inf(vector<LT> vec)
{
    LT max_el = 0;
    for (auto el : vec)
        if (max_el < fabs(el))
            max_el = fabs(el);
    return max_el;
}

template<typename LT>
LT vec_norm1(vector<LT> vec)
{
    LT res = 0;
    for (auto el : vec)
        res += fabs(el);
    return res;
}

//?????? ???????
template<typename LT>
void out(vector<LT> vec)
{
    int n = vec.size();
    for (int i = 0; i < n; ++i)
    {
        cout << fixed << setprecision(6) << setw(8) << setfill(' ') << vec[i] << "  ";
    }
    cout << endl;

}

template<typename LT>
LT det(vector<vector<LT>> matrix)
{
    int dim_size = matrix.size();
    if (dim_size == matrix[0].size())
    {
        if (dim_size == 2)
        {
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        }
        else if (dim_size == 3)
        {
            return	matrix[0][0] * matrix[1][1] * matrix[2][2] +
                      matrix[2][0] * matrix[0][1] * matrix[1][2] +
                      matrix[0][2] * matrix[1][0] * matrix[2][1] -
                      matrix[2][0] * matrix[1][1] * matrix[0][2] -
                      matrix[1][0] * matrix[0][1] * matrix[2][2] -
                      matrix[0][0] * matrix[2][1] * matrix[1][2];
        }
    }
    return;
}

//???????????????? ???????
template<typename LT>
vector<vector<LT>> transpose(vector<vector<LT>> matr)
{
    int n = matr[0].size();
    vector<vector<LT>> product(matr);
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (i != j)
            {
                product[i][j] = matr[j][i];
                product[j][i] = matr[i][j];
            }
    return product;
}