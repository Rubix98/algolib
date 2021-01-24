// Task: CALCAREA - Calculate the Area (SPOJ)
// URL: https://www.spoj.com/problems/CALCAREA/
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

VD rqe(double a, double b, double c) {
    if (ISZERO(a)) return ISZERO(b) ? VD() : VD(1, -c/b);
    else {
        double dlt = b*b - 4*a*c;
        if (ISZERO(dlt)) return VD(1, -b / (2*a));
        else if (dlt < 0) return VD();
        else {
            VD res;
            res.PB((-b-sqrt(dlt)) / (2*a));
            res.PB((-b+sqrt(dlt)) / (2*a));
            return res;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Point {
    double x, y;
    
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
    double a, b;
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

struct Line {
    double a, b, c;
    int pts;
    VP P;
    
    Line(): a(0), b(0), c(0), pts(0) {}
    
    Line(double a, double b, double c): a(a), b(b), c(c), pts(0) {}
    
    Line(Point p1, Point p2, int pts = 0): pts(pts), P({p1, p2}){
        double dx = p1.x - p2.x;
        double dy = p1.y - p2.y;
        a = -dy;
        b = dx;
        c = -(a*p1.x + b*p1.y);
        if (pts == 1) {
            P[1].x = p1.x < p2.x ? INF : -INF;
            P[1].y = p1.y < p2.y ? INF : -INF;
        }
    }
    
    double dist(Point p) {
        if (pts) {
            REP(i, 2) {
                Vector v1(P[i], P[1-i]), v2(P[i], p);
                if (min(v1.alpha(v2), v2.alpha(v1)) > M_PI_2) {
                    return min(p.dist(P[0]), p.dist(P[1]));
                }
            }
        }
        return abs(a*p.x + b*p.y + c) / sqrt(a*a + b*b);
    }
    
    double length() {
        return (pts == 2) ? P[0].dist(P[1]) : INF;
    }
    
    Line perpendicular(Point p) {
        return Line(b, -a, -b*p.x + a*p.y);
    }
    
    Line parallel(Point p) {
        return Line(a, b, -a*p.x - b*p.y);
    }
    
    bool hasPoint(Point p) {
        return ISZERO(dist(p));
    }
    
    VP getPointsWithX(double x) {
        if (ISZERO(b)) return {};
        Point p = Point(x, -(a*x+c)/b);
        return hasPoint(p) ? VP(1, p) : VP();
    }
    
    VP getPointsWithY(double y) {
        if (ISZERO(a)) return {};
        Point p = Point(-(b*y+c)/a, y);
        return hasPoint(p) ? VP(1, p) : VP();
    }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Circle {
    double a, b, c;
    double r;
    Point s;
    
    Circle(): a(0), b(0), c(0), r(0) {}
    
    Circle(double a, double b, double c): a(a), b(b), c(c) {
        s = Point(a, b) * (-0.5);
        r = sqrt(s.x*s.x + s.y*s.y - c);
    }
    
    Circle(Point s, double r): r(r), s(s){
        a = -2*s.x;
        b = -2*s.y;
        c = s.x*s.x + s.y*s.y - r*r;
    }
    
    Circle(Point s, Point p): s(s) {
        *this = Circle(s, s.dist(p));
    }
    
    double area() {
        return M_PI * r * r;
    }
    
    double circuit() {
        return 2 * M_PI * r;
    }
    
    VP getPointsWithX(double x) {
        VP res;
        FOREACH(y, rqe(1, b, x*x + a*x + c)) {
            res.PB(Point(x, y));
        }
        return res;
    }
    
    VP getPointsWithY(double y) {
        VP res;
        FOREACH(x, rqe(1, a, y*y + b*y + c)) {
            res.PB(Point(x, y));
        }
        return res;
    }
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Polygon {
    VP P;
    vector <Line> E;
    int n;
    
    Polygon(VP P): P(P), n(P.size()) {
        REP(i, n) {
            E.PB(Line(P[i], P[(i+1)%n], 2));
        }
    }
    
    double area() {
        double res = 0;
        REP(i, n-2) {
            Vector v1 = Vector(P[n-1], P[i]);
            Vector v2 = Vector(P[n-1], P[i+1]);
            res += v1.det(v2);
        }
        return abs(res)/2;
    }
    
    double circuit() {
        double res = 0;
        FOREACH(e, E) {
            res += e.length();
        }
        return res;
    }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main() {
    FASTIO();
    int n;
    cin >> n;
    VP points;
    while (n--) {
        int x, y;
        cin >> x >> y;
        points.PB(Point(x, y));
    }
    Polygon polygon(points);
    cout << round(polygon.area()) << endl;
}


