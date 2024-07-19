#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <queue>
#include <deque>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <list>
#include <bitset>
#include <complex>

using namespace std;

typedef long long ll;
typedef long double ld;
typedef pair<int, int> pii;
typedef pair<ll, ll> pLL;

#define FOR(i, a, b) for (int i(a); i < int(b); ++i)
#define ROF(i, a, b) for (int i(a); i >= int(b); --i)
#define FOREACH(i, c) for (__typeof((c).begin()) i = (c).begin(); i != (c).end(); i++)
#define ALL(c) (c).begin(), (c).end()
#define UNIQUE(c) (c).resize(unique(ALL(c)) - (c).begin())
#define pb push_back
#define mp make_pair
#define X first
#define Y second

const int MAX_N = 100000 + 20;
pii a[MAX_N];
int n;

int main()
{
	ios::sync_with_stdio(false);
	cin >> n;
	FOR(i, 0, n)
	cin >> a[i].X >> a[i].Y;
	sort(a, a + n);
	FOR(i, 1, n)
	if (a[i].Y < a[i - 1].Y)
	{
		cout << "Happy Alex" << endl;
		return 0;
	}
	cout << "Poor Alex" << endl;
	return 0;
}
