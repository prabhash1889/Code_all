#include <bits/stdc++.h>
using namespace std;
#include <iostream>

typedef long long ll;
typedef vector<int> vi;
ll n, m, k, q, l, r, x, y, z;
const ll template_array_size = 1e6 + 585;
ll a[template_array_size];
ll b[template_array_size];
ll c[template_array_size];
string s, t;
ll ans = 0;

void solve(int tc = 0)
{
	int x, y, x1, y1;
	cin >> x >> y >> x1 >> y1;
	if (abs(x - x1) == abs(y - y1))
	{
		if (abs(x - x1) == 0 and abs(y - y1) == 0)
		{
			cout << "YES" << '\n';
		}
		else
			cout << "NO" << '\n';
	}
	else
		cout << "YES" << '\n';
}

int main()
{

	int tc = 1;
	cin >> tc;

	// cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}