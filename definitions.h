#ifndef EVOLBIO_DEFINITIONS_H
#define EVOLBIO_DEFINITIONS_H

#pragma GCC target ("avx2")
#pragma GCC optimize ("O3")
#pragma GCC optimize ("unroll-loops")
#include <iostream>
#include <random>
#include <cmath>
#include <chrono>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <array>
#include <cassert>

void __print(bool x) {std::cerr << (x ? "true" : "false");}void __print(int x) {std::cerr << x;} void __print(long x) {std::cerr << x;} void __print(long long x) {std::cerr << x;} void __print(unsigned x) {std::cerr << x;} void __print(unsigned long x) {std::cerr << x;} void __print(unsigned long long x) {std::cerr << x;} void __print(float x) {std::cerr << x;} void __print(double x) {std::cerr << x;} void __print(long double x) {std::cerr << x;} void __print(char x) {std::cerr << '\'' << x << '\'';} void __print(const char *x) {std::cerr << '"' << x << '"';} void __print(const std::string &x) {std::cerr << '"' << x << '"';} template<typename T, typename V> void __print(const std::pair<T, V> &x) {std::cerr << '{'; __print(x.first); std::cerr << ','; __print(x.second); std::cerr << '}';} template<typename T> void __print(const T &x) {int f = 0; std::cerr << '{'; for (auto &i: x) std::cerr << (f++ ? "," : ""), __print(i); std::cerr << "}";} void _print() {std::cerr << "]\n";} template <typename T, typename... V> void _print(T t, V... v) {__print(t); if (sizeof...(v)) std::cerr << ", "; _print(v...);}
#define D(x...) std::cerr << "[" << #x << "] = ["; _print(x);
#define P(x...) std::cout << "[" << #x << "] = [" << x << "]" << std::endl;

typedef long long ll;
typedef unsigned long long ull;
typedef long double ld;
typedef std::pair<ll, ll> pl;
typedef std::vector<ll> vl;
typedef std::set<ll> sl;
typedef std::vector<pl> vpl;
typedef std::vector<vl> vvl;
typedef std::vector<ld> vld;
typedef std::vector<vld> vvld;
typedef std::vector<bool> vb;
typedef std::vector<vb> vvb;
typedef std::vector<std::string> vs;
#define F(i, end) for(ll i = 0; i < (end); ++i)
#define FS(i, start, end) for(ll i = (start); i < (end); ++i)
#define A(x) (x).begin(), (x).end()
#define X first
#define Y second
#define TR(x, v) for (const auto &(x) : v)
template <typename T> std::vector<T> set_to_vector(std::set<T> s) { return {A(s)}; }
#define DS(x) D(set_to_vector(x))

#define LOCI 4

#define SEX 0
// MALE = 1
#define MALE(x) (x.gene[SEX][0] + x.gene[SEX][1] > 0)
#define FEMALE (!MALE(x))

#define TYPE_1 0
#define TYPE_2 1
#define FIXED_RAND true

std::mt19937 generator(FIXED_RAND ? 0 : std::chrono::steady_clock::now().time_since_epoch().count());
std::uniform_real_distribution<ld> uniform_distribution_0_1(0.0,1.0);
std::uniform_int_distribution<ll> coin_flip(0,1);

ull codon_to_mask[4] = {0b0001, 0b0010, 0b0100, 0b1000};


#endif //EVOLBIO_DEFINITIONS_H
