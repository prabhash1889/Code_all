#include <bits/stdc++.h>
using namespace std;

void solve(int tc = 1)
{
	string s;
	cin >> s;
	string R = s;
	sort(R.begin(), R.end());
	if (R != s)
	{
		cout << "YES" << '\n';
		cout << R << '\n';
		return;
	}
	reverse(R.begin(), R.end());

	if (R != s)
	{
		cout << "YES" << '\n';
		cout << R << '\n';
		return;
	}

	cout << "NO" << '\n';
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
