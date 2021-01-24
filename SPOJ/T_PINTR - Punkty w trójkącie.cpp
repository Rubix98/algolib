// Task: T_PINTR - Punkty w trójkącie (SPOJ)
// URL: https://pl.spoj.com/problems/T_PINTR/
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

struct Point {
    int x, y;
    
    Point(): x(0), y(0) {}
    
    Point(double x, double y): x(x), y(y) {}
    
    double dist(Point p = Point()) {
        p = p - *this;
        return sqrt(p.x*p.x + p.y*p.y);
    }
    
    Point rotate(double a, Point c = Point()) {
        Point p = *this - c;
        return c + Point(p.x*cos(a) - p.y*sin(a), p.x*sin(a) + p.y*cos(a));
    }
    
    Point operator +(Point p) {
        return Point(x + p.x, y + p.y);
    }
    
    Point operator -(Point p) {
        return Point(x - p.x, y - p.y);
    }
    
    Point operator *(double a) {
        return Point(a*x, a*y);
    }
    
    bool operator <(Point p) {
        return MP(x, y) < MP(p.x, p.y);
    }
    
    bool operator ==(Point p) {
        p = p - *this;
        return ISZERO(p.x) && ISZERO(p.y);
    }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Vector {
    int a, b;
    Point s;
    
    Vector(): a(0), b(0) {}
    
    Vector(double a, double b): a(a), b(b) {}
    
    Vector(Point p1, Point p2): s(p1) {
        p1 = p2 - p1;
        a = p1.x;
        b = p1.y;
    }
    
    double length() {
        return Point(a, b).dist();
    }
    
    double scalarProduct(Vector v) {
        return a*v.a + b*v.b;
    }
    
    double det(Vector v) {
        return a*v.b - b*v.a;
    }
    
    double alpha(Vector v) {
        double res = acos(scalarProduct(v) / (length() * v.length()));
        if (det(v) >= 0) {
            return res;
        }
        else {
            return 2*M_PI - res;
        }
    }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main() {
    FASTIO();
    int x1, y1, x2, y2, x3, y3, x4, y4;
    cin >> x1 >> y1 >> x2 >> y2 >> x3 >> y3 >> x4 >> y4;
    Point p1(x1, y1), p2(x2, y2), p3(x3, y3), p4(x4, y4);
    Vector v1(p1, p2), v2(p3, p4);
    if (v1.a == v2.a && v1.b == v2.b) printf("tak [%d;%d]\n", v1.a, v2.b);
    else printf("nie\n");
}
