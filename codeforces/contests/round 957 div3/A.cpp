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
	vector<ll> a(3);
	for (ll i = 0; i < 3; i++)
	{
		cin >> a[i];
	}
	sort(a.begin(), a.end());
	a[0]++;
	for (int i = 0; i < 4; i++)
	{
		if (a[0] > a[1])
		{
			a[1]++;
		}
		else if (a[1] > a[2])
		{
			a[2]++;
		}
		else if (a[0] <= a[1])
		{
			a[0]++;
		}
		else if (a[1] <= a[2])
		{
			a[1]++;
		}
	}

	cout << a[0] * a[1] * a[2] << '\n';
}

int main()
{

	int tc = 1;
	cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}