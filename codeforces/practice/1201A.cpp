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
	int n, m;
	cin >> n >> m;
	vector<map<int, int>> ans(m);
	for (int i = 0; i < n; i++)
	{
		string s;
		cin >> s;
		for (int j = 0; j < m; j++)
		{
			ans[j][s[j]]++;
		}
	}
	vector<int> points(m);
	int total = 0;
	for (int i = 0; i < m; i++)
	{
		cin >> points[i];
		int most = 0;
		for (auto &it : ans[i])
		{
			most = max(most, it.second);
		}
		total += points[i] * most;
	}
	cout << total << '\n';
}

int main()
{

	int tc = 1;
	// cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}