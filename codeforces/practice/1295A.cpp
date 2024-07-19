#include <bits/stdc++.h>
using namespace std;

void solve(int tc = 1)
{
	int n;
	cin >> n;
	if (n % 2 != 0)
	{
		cout << "7";
		n -= 3;
	}

	while (n > 0)
	{
		cout << "1";
		n -= 2;
	}

	cout << '\n';
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
