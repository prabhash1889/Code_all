#include <bits/stdc++.h>
using namespace std;

void solve(int tc = 1)
{
	string s;
	cin >> s;
	int ans = 1;
	for (int i = 1; i < s.size(); i++)
	{
		if (s[i - 1] != s[i])
		{
			ans++;
		}
	}
	if (ans == 1)
	{
		cout << ans << '\n';
	}
	else if (ans > 2)
	{
		cout << ans - 1 << '\n';
	}
	else
	{
		if (s[0] == '1')
		{
			cout << ans << '\n';
		}
		else
		{
			cout << ans - 1 << '\n';
		}
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
