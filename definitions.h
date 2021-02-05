#ifndef EVOLBIO_DEFINITIONS_H
#define EVOLBIO_DEFINITIONS_H

#pragma GCC target ("avx2")
#pragma GCC optimize ("O3")
#pragma GCC optimize ("unroll-loops")
#include <bits/stdc++.h>
#include <random>
#include <math.h>
#include <chrono>

void __print(bool x) {std::cerr << (x ? "true" : "false");}void __print(int x) {std::cerr << x;} void __print(long x) {std::cerr << x;} void __print(long long x) {std::cerr << x;} void __print(unsigned x) {std::cerr << x;} void __print(unsigned long x) {std::cerr << x;} void __print(unsigned long long x) {std::cerr << x;} void __print(float x) {std::cerr << x;} void __print(double x) {std::cerr << x;} void __print(long double x) {std::cerr << x;} void __print(char x) {std::cerr << '\'' << x << '\'';} void __print(const char *x) {std::cerr << '"' << x << '"';} void __print(const std::string &x) {std::cerr << '"' << x << '"';} template<typename T, typename V> void __print(const std::pair<T, V> &x) {std::cerr << '{'; __print(x.first); std::cerr << ','; __print(x.second); std::cerr << '}';} template<typename T> void __print(const T &x) {int f = 0; std::cerr << '{'; for (auto &i: x) std::cerr << (f++ ? "," : ""), __print(i); std::cerr << "}";} void _print() {std::cerr << "]\n";} template <typename T, typename... V> void _print(T t, V... v) {__print(t); if (sizeof...(v)) std::cerr << ", "; _print(v...);}
#define D(x...) std::cerr << "[" << #x << "] = ["; _print(x);

#define P(x...) std::cout << "[" << #x << "] = [" << x << "]" << std::endl;

typedef long long int ll;
typedef long double ld;
typedef std::pair<ll, ll> pl;
#define F(i, end) for(ll i = 0; i < (end); ++i)
#define A(x) (x).begin(), (x).end()
#define X first
#define Y second
#define tr(x, v) for (auto x : v)

#define LOCI 4

#define SEX 0
// MALE = 1
#define MALE(x) (x.gene[SEX][0] + x.gene[SEX][1] > 0)
#define FEMALE (!MALE(x))

#define TYPE_1 0
#define TYPE_2 1

std::mt19937 generator(std::chrono::steady_clock::now().time_since_epoch().count());
std::uniform_real_distribution<ld> uniform_distribution_0_1(0.0,1.0);
std::uniform_int_distribution<ll> coin_flip(0,1);

#endif //EVOLBIO_DEFINITIONS_H
