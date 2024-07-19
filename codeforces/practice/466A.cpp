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
	int n, m, a, b;
	cin >> n >> m >> a >> b;
	int ans = 0;
	int k;
	if (b <= a * m)
	{
		k = n % m;
		int tt = (n - k) / m;
		ans += (tt * b);
		if (k * a > b)
		{
			ans += (b);
		}
		else
			ans += (k * a);
		cout << ans << '\n';
		return;
	}
	else if (b > a * m)
	{
		ans += n * a;
		cout << ans << '\n';
		return;
	}
}

int main()
{

	int tc = 1;
	// cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}