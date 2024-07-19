#include <bits/stdc++.h>
using namespace std;

void solve(int tc = 1)
{
	int r, g, b;
	cin >> r >> g >> b;
	string s;
	cin >> s;
	if (s == "Red")
	{
		cout << min(g, b) << '\n';
		return;
	}
	if (s == "Green")
	{
		cout << min(r, b) << '\n';
		return;
	}
	if (s == "Blue")
	{
		cout << min(r, g) << '\n';
		return;
	}
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
