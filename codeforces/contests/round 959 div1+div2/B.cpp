#include <bits/stdc++.h>
using namespace std;

void solve(int tc = 1)
{
	int n;
	cin >> n;
	string s, t;
	cin >> s >> t;
	for (int i = 0; i < s.size() && s[i] == '0'; i++)
	{
		if (t[i] != '0')
		{
			cout << "NO" << '\n';
			return;
		}
	}
	cout << "YES" << '\n';
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
