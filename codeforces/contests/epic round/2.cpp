#include <bits/stdc++.h>
using namespace std;
#include <iostream>

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
	int n;
	cin >> n;
	vector<int> a;
	int flag=0;
	for (int i = 0; i < n; i++)
	{
		int u;
		cin >> u;
		a.push_back(u);
		
	}
	ll ans = 0;
	// vector<int> b;
	
	for (int i = 0; i < n; i++)
	{
		if (a[i + 1] < a[i])
		{
			ans += abs(a[i + 1] - a[i]);
			a[i + 1] = a[i];
		}
	}
	cout << ans + 1 << '\n';
}

int main()
{

	int tc = 1;
	cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}