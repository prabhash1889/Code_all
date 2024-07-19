#include <bits/stdc++.h>
using namespace std;

void solve(int tc = 1)
{
	int n;
	cin >> n;
	int x = 0;
	for (int i = 0; i < n; i++)
	{
		string s;
		cin >> s;
		for (int j = 0; j < 3; j++)
		{
			if (s[j] == '+')
			{
				x++;
				break;
			}
			else if (s[j] == '-')
			{
				x--;
				break;
			}
		}
	}
	cout << x << '\n';
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
