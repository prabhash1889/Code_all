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
	ll n;
	cin >> n;
	vector<ll> A(n), B(n + 1);
	for (int i = 0; i < n; i++)
	{
		cin >> A[i];
	}
	for (int i = 0; i < n + 1; i++)
	{
		cin >> B[i];
	}
	// vector<ll> C;
	ll ans = 0;
	ll mn = INT_MAX;
	for (ll i = 0; i < n; i++)
	{
		ll a = abs(A[i] - B[i]);
		ans += a;
		mn = min(mn, abs(A[i] - B[n]));
		mn = min(mn, abs(B[i] - B[n]));
		if (min(A[i], B[i]) <= B[n] && B[n] <= max(A[i], B[i]))
		{
			mn = 0;
		}
	}

	cout << ans + mn + 1 << '\n';
}

int main()
{

	int tc = 1;
	cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}