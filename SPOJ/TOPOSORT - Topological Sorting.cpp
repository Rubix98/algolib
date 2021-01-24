// Task: TOPOSORT - Topological Sorting (SPOJ)
// URL: https://www.spoj.com/problems/TOPOSORT/
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

TEMPL struct Edge {
    int a, b;
    T w;
    ET(int a, int b, T w = 1): a(a), b(b), w(w) {}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEMPL struct Graph {
    VVE adj;
    VE E;
    int n;
    T inf;
    Graph(int n = 0, T inf = INF): n(n), inf(inf) {
        adj.resize(n);
    }
    void addEdgeD(int a, int b, int w = 1) {
        ET *e = new ET(a, b, w);
        E.PB(e);
        adj[a].PB(e);
    }
    void addEdgeD(ET e) {
        addEdgeD(e.a, e.b, e.w);
    }
    void addEdge(int a, int b, int w = 1) {
        addEdgeD(a, b, w);
        if (a != b) adj[b].PB(new ET(b, a, w));
    }
    void addEdge(ET e) {
        addEdge(e.a, e.b, e.w);
    }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEMPL void readGraph(GT &g, bool isD, bool isW, int id = 0) {
    int n, m, a, b;
    T w = 1;
    cin >> n >> m;
    g = GT(n);
    REP(i, m) {
        cin >> a >> b;
        a -= id;
        b -= id;
        if (isW) cin >> w;
        if (isD) g.addEdgeD(a, b, w);
        else g.addEdge(a, b, w);
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEMPL VI topoSort(GT &g) {
    VI res;
    priority_queue<int> q;
    VI rank = VI(g.n, 0);
    FOREACH(e, g.E)
        rank[e->b]++;
    REP(i, g.n) 
        if (rank[i] == 0)
            q.push(-i);
    while(!q.empty()) {
        int x = -q.top(); q.pop();
        res.PB(x);
        FOREACH(e, g.adj[x])
            if (--rank[e->b] == 0) 
                q.push(-e->b);
    }
    return res;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main() {
    FASTIO();
    Graph<LL> g;
    readGraph(g, true, false, 1);
    VI result = topoSort(g);
    if ((int)result.size() == g.n) {
        FOREACH(x, result) 
            cout << x + 1 << " ";
        cout << endl;
    }
    else cout << "Sandro fails." << endl;
}



