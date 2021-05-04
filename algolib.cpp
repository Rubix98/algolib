// Autor: Marcin Wojdat
// Data: 24.01.2021
// Opis: Zbiór wszystkich implementacji z pracy "Biblioteka algorytmów i struktur danych używanych w zawodach programistycznych"

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

ULL gcd(LL a, LL b) {
    a = abs(a);
    b = abs(b);
    if (a < b) swap(a, b);
    return b ? gcd(b, a % b) : a;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ULL lcm(int a, int b) {
    a = abs(a);
    b = abs(b);
    return (a || b) ? (ULL) a / gcd(a, b) * b : 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

UI invMod(LL a, UI m) {
    LL b = m, x = 1 % m, y = 0;
    a = MOD(a, m); 
    while (a > 1) {
        if (!b) return 0;
        x -= a / b * y;
        a %= b;
        swap(x, y); 
        swap(a, b); 
    }
    return MOD(x, m); 
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

UI powMod(LL a, ULL b, UI m) {
    if (!b) return 1 % m;
    a = MOD(a, m);
    UI res = powMod((ULL) a * a % m, b >> 1, m);
    if (b & 1) res = res * a % m;
    return res;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

UI newton(int n, int k, UI p = INF) {
    if (n < 0 || k < 0 || n < k) return 0;
    ULL res = 1;
    LOOP(i, 0, min(k, n-k)) {
        res = res * (n-i) % p;
        res = res * invMod(i+1, p) % p;
    }
    return res;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

VB eratosthenes(int n) {
    VB res = VB(n, 1);
    res[0] = res[1] = 0;
    int sq = sqrt(n);
    LOOP(i, 2, sq+1)
        if (res[i]) 
            for (int j = i+i; j < n; j += i)
                res[j] = 0;
    return res;
}

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

template <typename T, typename T2 = T> struct PrefixSums {
    vector<T2> sumTab;
    PrefixSums(VT &tab = VT()) {
        sumTab.PB(0);
        FOREACH(x, tab) sumTab.PB(sumTab.back() + x);
    }
    
    T2 sum(int a, int b) {
        return a < b ? sumTab[b] - sumTab[a] : 0;
    }
};


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

struct FU {
    VI P, R;
    int n;

    FU(int n = 0): n(n) {
        P.resize(n);
        R.resize(n, 0);
        REP(i, n) P[i] = i;
    }

    int find(int a) {
        return P[a] == a ? a : (P[a] = find(P[a]));
    }

    bool Union(int a, int b) {
        a = find(a);
        b = find(b);
        if (a == b) return false;
        if (R[a] < R[b]) P[a] = b;
        else P[b] = a;
        if (R[a] == R[b]) R[a]++;
        return true;
    }

    bool isConnected(int a, int b) {
        return find(a) == find(b);
    }
};

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

TEMPL void makeGraph(GT &g, VE E, bool isD, bool isW) {
    g = GT(g.n);
    FOREACH(e, E) {
        ET edge(e->a, e->b, isW ? e->w : 1);
        if (isD) g.addEdgeD(edge);
        else g.addEdge(edge);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEMPL T sumEdges(VE &E) {
    T res = 0;
    FOREACH(e, E)
        if (e != NP)
            res += e->w;
    return res;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEMPL void dfsR(VE &res, GT &g, int v, int r, ET *e = NP) {
    if (res[v] == NP) {
        res[v] = e;
        FOREACH(e, g.adj[v]) {
            if (e->b != r) {
                dfsR(res, g, e->b, r, e);
            }
        }
    }
}

TEMPL VE dfs(GT &g, int v = -1) {
    VE res(g.n, NP);
    if (v != -1)
        dfsR(res, g, v, v);
    else 
        REP(i, g.n) 
            dfsR(res, g, i, i);
    return res;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEMPL PVTVE bfs(GT &g, int r) {
    VT dist(g.n, g.inf);
    VE E(g.n, NP);
    queue<int> q;
    dist[r] = 0;
    q.push(r);
    while (!q.empty()) {
        int v = q.front();
        q.pop();
        FOREACH(e, g.adj[v]) {
            T d = dist[v] + 1;
            if (d < dist[e->b]) {
                dist[e->b] = d;
                E[e->b] = e;
                q.push(e->b);
            }
        }
    }
    return MP(dist, E);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEMPL PVTVE dijkstra(GT &g, int r) {
    VT dist(g.n, g.inf);
    VE E(g.n, NP);
    priority_queue<pair<T, int> > q;
    dist[r] = 0;
    q.push(MP(0, r));
    while (!q.empty()) {
        auto x = q.top();
        q.pop();
        if (x.ST > dist[x.ND]) continue;
        FOREACH(e, g.adj[x.ND]) {
            T d = -x.ST + e->w;
            if (d < dist[e->b]) {
                dist[e->b] = d;
                E[e->b] = e;
                q.push(MP(-d, e->b));
            }
        }
    }
    return MP(dist, E);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEMPL VI topoSort(GT &g) {
    VI res;
    VI rank = VI(g.n, 0);
    FOREACH(e, g.E)
        rank[e->b]++;
    REP(i, g.n) 
        if (rank[i] == 0)
            res.PB(i);
    REP(i, res.size())
        FOREACH(e, g.adj[res[i]])
            if (--rank[e->b] == 0) 
                res.PB(e->b);
    return res;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEMPL VE mst(GT &g) {
    VE res;
    VE E = g.E;
    sort(ALL(E), [](ET *e1, ET *e2){return e1->w < e2->w;});
    FU fu(g.n);
    FOREACH(e, E)
        if (fu.Union(e->a, e->b))
            res.PB(e);
    return res;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

VI kmp(string s) {
    VI res(s.size()+1, 0);
    int x = 0;
    LOOP(i, 2, res.size()) {
        while (x && s[i-1] != s[x]) x = res[x];
        if (s[i-1] == s[x]) x++;
        res[i] = x;
    }
    return res;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

VI findPattern(string p, string t) {
    VI res;
    VI tab = kmp(p+"#"+t);
    int ps = p.size();
    REP(i, tab.size()) 
        if (tab[i] == ps) 
            res.PB(i-2*ps-1);
    return res;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int minPeriod(string &s) {
    return s.size() - kmp(s).back();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct TNode {
    vector<TNode*> cld;
    TNode* fail = NP;
    TNode* suf = NP;
    VI P;
    
    TNode(int ch){
        cld.resize(ch, NP);
    }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct Trie {
    TNode *root;
    VS P;
    int ch;
    
    Trie(VS P, int ch = 62): ch(ch) {
        root = new TNode(ch);
        FOREACH(p, P) insert(p);
        setFails();
    }
    
    int id(char c) {
        if ('a' <= c && c <= 'z') return c-'a';
        else if ('A' <= c && c <= 'Z') return c-'A'+26;
        else return c-'0'+52;
    }
    
    void insert(string s) {
        TNode *node = root;
        FOREACH(c, s) {
            if (node->cld[id(c)] == NP)
                node->cld[id(c)] = new TNode(ch);
            node = node->cld[id(c)];
        }
        node->P.PB(P.size());
        P.PB(s);
    }
    
    void setFails() {
        queue <TNode*> q;
        root->fail = root;
        q.push(root);
        while (!q.empty()) {
            TNode* node = q.front();
            q.pop();
            REP(i, ch) {
                TNode *cld = node->cld[i];
                if (cld != NP) {
                    if (node != root) {
                        TNode *tmp = node->fail;
                        while (tmp != root && tmp->cld[i] == NP)
                            tmp = tmp->fail;
                        cld->fail = (tmp->cld[i] == NP) ? 
                            root : tmp->cld[i];
                        cld->suf = cld->fail->P.empty() ?
                            cld->fail->suf : cld->fail;
                        q.push(cld);
                    }
                    else {
                        cld->fail = root;
                        q.push(cld);
                    }
                }
            }
        }
    }
    
    VVI search(string s) {
        VVI res(P.size());
        TNode *node = root;
        REP(i, s.size()) {
            char c = s[i];
            while (node != root && node->cld[id(c)] == NP)
                node = node->fail;
            if (node->cld[id(c)] != NP)
                node = node->cld[id(c)];
            TNode *tmp = node;
            do {
                FOREACH(p, tmp->P)
                    res[p].PB(i-P[p].size()+1);
                tmp = tmp->suf;
            } while (tmp != NP);
        }
        return res;
    }
};

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
