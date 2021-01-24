// Task: FR_09_19 - Słówka, słóweczka (SPOJ)
// URL: https://pl.spoj.com/problems/FR_09_19/
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

char tab[200][200];

int main() {
    FASTIO();
    VS patterns;
    readTab(patterns);
    Trie trie(patterns);
    int n;
    cin >> n;
    REP(i, n) {
        REP(j, n) {
            cin >> tab[i][j];
        }
    }
    
    string s1, s2;
    REP(i, n) {
        REP(j, n) {
            s1 += tab[i][j];
            s2 += tab[j][i];
        }
        s1 += "x";
        s2 += "x";
    }
    
    auto res1 = trie.search(s1), res2 = trie.search(s2);
    int answer = 0;
    REP(i, patterns.size()) {
        answer += res1[i].size() + res2[i].size();
    }
    cout << answer << endl;
    
}

