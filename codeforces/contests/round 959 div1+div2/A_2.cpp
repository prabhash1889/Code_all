#include <bits/stdc++.h>
using namespace std;

void solve(int tc = 1)
{
	int n, m;
	cin >> n >> m;
	vector<vector<int>> a(n, vector<int>(m));
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			cin >> a[i][j];
		}
	}
	if (n * m == 1)
		cout << -1 << '\n';
	else
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				cout << a[(i + 1) % n][(j + 1) % m] << " ";
			}
		}
		cout << '\n';
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
