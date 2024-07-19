#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef vector<int> vi;
ll n, m, k, q, l, r, x, y, z;
const ll template_array_size = 1e6 + 585;
ll a[template_array_size];
ll b[template_array_size];
ll c[template_array_size];
string s, t;
ll ans = 0;

void solve(int tc = 0)
{
	int n, k;
	cin >> n >> k;
	string s;
	cin >> s;

	vector<pair<char, int>> c(n);
	for (int i = 0; i < n; ++i)
		c[i] = make_pair(s[i], i);
	sort(c.begin(), c.end());
	sort(c.begin() + k, c.end(), [&](const pair<char, int> &a, const pair<char, int> &b)
		 { return a.second < b.second; });

	for (int i = k; i < n; ++i)
		cout << c[i].first;
	cout << '\n';
}

int main()
{

	int tc = 1;
	// cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}