#include <bits/stdc++.h>
using namespace std;

void solve(int tc = 1)
{
	int n;
	cin >> n;
	string s, t;
	cin >> s >> t;
	int ks = 0, kt = 0;
	for (int i = 0; i < s.size(); i++)
	{
		if (s[i] == '0')
		{
			ks++;
		}
		else
			break;
	}
	for (int j = 0; j < t.size(); j++)
	{
		if (t[j] == '0')
		{
			kt++;
		}
		else
			break;
	}
	cout << (ks > kt ? "NO" : "YES") << '\n';
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
