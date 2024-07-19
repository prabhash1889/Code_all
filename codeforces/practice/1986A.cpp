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
	vi a;
	for (int i = 0; i < 3; i++)
	{
		int x;
		cin >> x;
		a.push_back(x);
	}
	sort(a.begin(), a.end());
	cout << (a[2] - a[0]) << '\n';
}

int main()
{

	int tc = 1;
	cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}