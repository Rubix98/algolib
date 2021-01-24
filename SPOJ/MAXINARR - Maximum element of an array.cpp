// Task: MAXINARR - Maximum element of an array (SPOJ)
// URL: https://www.spoj.com/problems/MAXINARR/
// Author: Marcin Wojdat
// Date: 24.01.2021
// Algorithms from: https://github.com/Rubix98/algorithms

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Biblioteczki
#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <queue>
#include <stack>
#include <numeric>
using namespace std;
// Typy proste
typedef unsigned int        UI;
typedef long long           LL;
typedef unsigned long long  ULL;
// Pary
typedef pair<int, int>  PI;
#define MP              make_pair
#define ST              first
#define ND              second
// Vectory
typedef vector<int>  VI;
typedef vector<bool> VB;
typedef vector<double> VD;
typedef vector<string> VS;
typedef vector<VI>   VVI;
typedef vector<PI>   VPI;
#define VT           vector<T>
#define VVT          vector<VT >
#define PB           push_back
#define SIZE(v)      (int(v.size()))
#define ALL(v)       v.begin(), v.end()
// Pętle
#define LOOP(i, a, b)  for (int i = (a); i < (int)(b); ++i)
#define FOREACH(x, v)  for (auto x: v)
#define REP(i, n)      LOOP(i, 0, n)
// Stałe
const int INF = 1e9+9;
const double EPS = 1e-9;
// Grafy
#define ET   Edge<T>
#define GT   Graph<T>
#define VE   vector<ET*>
#define VVE  vector<VE > 
#define PVTVE pair<VT, VE >
// Geometria
#define VP vector<Point>
#define ISZERO(x) abs(x) < EPS
// Inne
#define MOD(a, m) (a % m + m) % m
#define TEMPL template <typename T>
#define NP nullptr
#define BETWEEN(x, a, b) (((a) <= x && x <= (b)) || ((b) <= x && x <= (a)))

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Biblioteczki
#include <string>
#include <map>
#include <set>
#include <list>
#include <bitset>
#include <stack>
#include <queue>

// Vectory
typedef vector<UI>      VUI;
typedef vector<LL>      VLL;
typedef vector<ULL>     VULL;

// Pary
typedef pair<double, double>  PD;
typedef pair<string, string>  PS;
typedef pair<UI, UI>          PUI;
typedef pair<LL, LL>          PLL;
typedef pair<ULL, ULL>        PULL;

// Pętle
#define DLOOP(i, a, b)  for (int i = (a); i <= (b); ++i)
#define RLOOP(i, a, b)  for (int i = (a); i >= (b); --i)
#define TESTS(t)        int t; cin >> t; REP(i, t) 

// Stałe
const LL LLINF = 1e18+9;

// Konwersja char na int i odwrotnie
#define CTOI(c) (int(c))
#define ITOC(x) (char(x))

// I/O
#define FASTIO()  ios_base::sync_with_stdio(0); cin.tie(NULL);
TEMPL void printTab(VT &v) {
    cout << "[";
    REP(i, v.size())
        cout << (i ? ", " : "") << v[i];
    cout << "]\n";
}
TEMPL void readTab(VT &v) {
    int n;
    cin >> n;
    v.resize(n);
    REP(i, n) cin >> v[i];
}

// Debugowanie
#ifdef DEBUGGING
#define DEBUG(x) cout << #x << ": " << x << endl;
#else
#define DEBUG(x) 
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEMPL struct MinMax {
    VVT minTab, maxTab;
    VI lg;
    int n;
    T inf;

    MinMax(VT &tab = VT(), T inf = INF): n(tab.size()), inf(inf) {
        if (n) {
            minTab.resize(log2(n)+1);
            maxTab.resize(log2(n)+1);
            lg.resize(n+1);
            minTab[0] = maxTab[0] = tab;
            LOOP(i, 1, minTab.size()) {
                minTab[i].resize(n + 1 - (1<<i));
                maxTab[i].resize(n + 1 - (1<<i));
                REP(j, minTab[i].size()) {
                    int x = j+(1<<(i-1));
                    minTab[i][j] = min(minTab[i-1][j], minTab[i-1][x]);
                    maxTab[i][j] = max(maxTab[i-1][j], maxTab[i-1][x]);
                }
            }
            LOOP(i, 1, n+1)lg[i] = log2(i);
        }
    }
    
    T minimum(int a, int b) {
        if (a >= b) return inf;
        int x = lg[b-a];
        return min(minTab[x][a], minTab[x][b-(1<<x)]);
    }
    
    T maximum(int a, int b) {
        if (a >= b) return -inf;
        int x = lg[b-a];
        return max(maxTab[x][a], maxTab[x][b-(1<<x)]);
    }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main() {
    FASTIO();
    while (true) {
        VI tab;
        readTab(tab);
        if (tab.empty()) break;
        MinMax<int> mm(tab);
        cout << mm.maximum(0, tab.size()) << endl;
    } 
}



