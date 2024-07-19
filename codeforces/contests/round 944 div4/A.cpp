#include <bits/stdc++.h>
using namespace std;

void solve(int tc = 1)
{
	int x, y;
	cin >> x >> y;
	cout << min(x, y) << " " << max(x, y) << '\n';
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
