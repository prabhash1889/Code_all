#include <bits/stdc++.h>
using namespace std;

void solve(int tc = 1)
{
	int n, m;
	cin >> n >> m;
	int a[n][m];
	// vector<vector<int>> a;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			cin >> a[i][j];
		}
	}
	if (n == 1 and m == 1)
	{
		cout << -1 << '\n';
		return;
	}
	// vector<vector<int>> b;
	int b[n][m];
	if (n == 1)
	{
		for (int i = 0; i < m - 1; i++)
		{
			int temp = a[0][i + 1];
			a[0][i + 1] = a[0][i];
			a[0][i] = temp;
		}
	}
	else if (m == 1)
	{
		for (int i = 0; i < n - 1; i++)
		{
			int temp = a[i + 1][0];
			a[i + 1][0] = a[i][0];
			a[i][0] = temp;
		}
	}
	else
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m - 1; j++)
			{
				int temp = a[i][j + 1];
				a[i][j + 1] = a[i][j];
				a[i][j] = temp;
			}
		}
	}
	// cout << n * m << '\n';
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			cout << a[i][j] << " ";
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
