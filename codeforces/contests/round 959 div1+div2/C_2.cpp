#include <bits/stdc++.h>
using namespace std;

void solve(int tc = 1)
{
	int n, x;
	cin >> n >> x;
	vector<int> a(n);
	for (int i = 0; i < n; i++)
	{
		cin >> a[i];
	}
	vector<int64_t> pref(n + 1);
	for (int i = 0; i < n; i++)
	{
		pref[i + 1] = pref[i] + a[i];
	}
	int64_t ans = int64_t(n) * (n + 1) / 2;
	vector<int> dp(n + 1);
	for (int i = n - 1; i >= 0; i--)
	{
		auto it = lower_bound(pref.begin(), pref.end(), pref[i] + x + 1);
		int id = int(it - pref.begin());
		dp[i] = (id == n + 1 ? 0 : dp[id] + 1);
		ans -= dp[i];
	}
	cout << ans << '\n';
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
