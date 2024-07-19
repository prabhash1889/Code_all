#include <bits/stdc++.h>
using namespace std;
#define ll long long

void solve(int tc = 1)
{
	int n;
	ll x;
	cin >> n >> x;
	vector<ll> a(n + 1);
	for (int i = 1; i <= n; ++i)
		cin >> a[i];
	partial_sum(a.begin() + 1, a.end(), a.begin() + 1);
	vector<int> dp(n + 2);
	for (int i = n - 1; i >= 0; --i)
	{
		int q = upper_bound(a.begin(), a.end(), a[i] + x) - a.begin();
		dp[i] = dp[q] + q - i - 1;
	}
	cout << accumulate(dp.begin(), dp.end(), 0ll) << '\n';
}
int main()
{
	ios::sync_with_stdio(false);
	cin.tie(0);

	int tc = 1;
	cin >> tc;
	for (int t = 0; t < tc; t++)
	{
		solve(t);
	}
}
