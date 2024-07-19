#include <bits/stdc++.h>
using namespace std;

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
	int n;
	cin >> n;
	int ans = 0;
	map<int, bool> mp;
	while (!mp.count(n))
	{
		mp[n] = 1;
		ans++;
		n++;
		while (n % 10 == 0)
		{
			n /= 10;
		}
	}
	cout << ans << '\n';
}

int main()
{

	int tc = 1;
	// cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}