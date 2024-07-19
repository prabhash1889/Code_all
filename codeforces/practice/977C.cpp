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
	int n, k;
	cin >> n >> k;
	vector<int64_t> a(n);
	for (int i = 0; i < n; i++)
	{
		cin >> a[i];
	}
	sort(a.begin(), a.end());
	if (k == 0)
	{
		ans = a[0] - 1;
	}
	else
		ans = a[k - 1];
	int cnt = 0;

	for (int i = 0; i < n; ++i)
		if (a[i] <= ans)
			++cnt;

	if (cnt != k || !(1 <= ans && ans <= 1000 * 1000 * 1000))
	{
		cout << -1 << '\n';
		return;
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