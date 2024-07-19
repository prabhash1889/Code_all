#include <bits/stdc++.h>
using namespace std;

void solve(int tc = 1)
{
	int M[5][5];
	int a, b;
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			cin >> M[i][j];

			if (M[i][j] == 1)
			{
				a = i + 1;
				b = j + 1;
			}
		}
	}
	cout << abs(a - 3) + abs(b - 3) << '\n';
}

int main()
{
	ios::sync_with_stdio(false);
	cin.tie(0);

	int tc = 1;
	// cin >> tc;
	for (int t = 0; t < tc; t++)
	{
		solve(t);
	}
}
