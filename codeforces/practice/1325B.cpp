#include <bits/stdc++.h>
using namespace std;

void solve(int tc = 1)
{
	int n;
	cin >> n;
	vector<int64_t> v(n);
	for (int i = 0; i < n; i++)
		cin >> v[i];
	sort(v.begin(), v.end());
	int res = 0;
	for (int i = 0; i < n; i++)
	{
		while (i < n - 1 && v[i] == v[i + 1])
			i++;

		res++;
	}
	if (res <= n)
	{
		cout << res << '\n';
		return;
	}
	else
	{
		cout << n << '\n';
		return;
	}
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
