#include <bits/stdc++.h>
using namespace std;

void solve(int tc = 1)
{
	int n;
	cin >> n;
	string s;
	cin >> s;
	int k1 = 0;
	int k0 = 0;
	for (int i = 0; i < n; i++)
	{
		if (s[i] == '1')
		{
			k1 += 1;
		}
		else
		{
			if (i == 0 || s[i - 1] != '0')
			{
				k0 += 1;
			}
		}
	}
	cout << (k1 > k0 ? "Yes" : "No") << '\n';
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
