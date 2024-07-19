#include <bits/stdc++.h>
using namespace std;

void solve(int tc = 1)
{
	int n;
	cin >> n;
	vector<int> a(n);
	for (auto &A : a)
		cin >> A;
	sort(a.begin(), a.end());
	if (a[0] == a[n - 1])
		cout << n << '\n';
	else
		cout << 1 << '\n';
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
