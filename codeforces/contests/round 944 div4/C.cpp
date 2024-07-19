#include <bits/stdc++.h>
using namespace std;

bool ordered(int x, int y, int z)
{
	return (x < y && y < z) || (z < x && x < y) || (y < z && z < x);
}

void solve(int tc = 1)
{
	int A, B, C, D;
	cin >> A >> B >> C >> D;
	cout << (ordered(A, C, B) ^ (ordered(A, D, B)) ? "YES" : "NO") << '\n';
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
