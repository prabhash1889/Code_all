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
	ll n, k;
	cin >> n >> k;
	vector<ll> a(k);
	ll ones = 0;
	for (ll i = 0; i < k; i++)
	{
		cin >> a[i];
		if (a[i] == 1)
			ones++;
	}
	sort(a.begin(), a.end());
	ll ans = 0;

	ans += ones;
	if (ans == n)
	{
		cout << n - 1 << '\n';
		return;
	}
	for (ll i = ones; i < k - 1; i++)
	{
		ans = ans + (a[i] + a[i] - 1);
	}
	cout << ans << '\n';
}

int main()
{

	int tc = 1;
	cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}