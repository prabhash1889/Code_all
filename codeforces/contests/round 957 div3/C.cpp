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
	ll n, m, k;
	cin >> n >> m >> k;
	ll tt = 0;
	for (int i = 0; (n - i) >= k; i++)
	{
		if ((n - i) >= k)
		{
			cout << n - i << " ";
		}
		tt++;
	}
	for (int i = tt; i < n - m; i++)
	{
		cout << n - tt << " ";
		tt++;
	}
	for (int i = 0; i < (n - tt); i++)
	{
		cout << i + 1 << " ";
	}
	cout << '\n';
}

int main()
{

	int tc = 1;
	cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}