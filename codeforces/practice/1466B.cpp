#include <bits/stdc++.h>
using namespace std;

void solve(int tc = 1)
{
	// int n;
	// cin >> n;
	// set<int> S;
	// for (int i = 1; i <= n; i++)
	// {
	// 	int v;
	// 	cin >> v;
	// 	if (S.count(v))
	// 	{
	// 		v++;
	// 	}
	// 	S.insert(v);
	// }
	// cout << S.size() << '\n';
	int n;
	cin >> n;
	int ans = 0, prev = 0, ok = 0;
	for (int i = 1; i <= n; i++)
	{
		int a;
		cin >> a;
		ans += (a != prev);
		ok |= (a == prev);
		if (prev + 1 < a)
		{
			ans += ok;
			ok = 0;
		}
		prev = a;
	}
	ans += ok;
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
