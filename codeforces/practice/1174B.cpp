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
	vector<int> a(n);
	int even = 0, odd = 0;
	for (int i = 0; i < n; i++)
	{
		cin >> a[i];
		if (a[i] % 2 == 0)
			even++;
		else
			odd++;
	}
	if (even == 0 or odd == 0)
	{
		for (int i = 0; i < n; i++)
			cout << a[i] << " ";
	}
	else
	{
		sort(a.begin(), a.end());
		for (int i = 0; i < n; i++)
		{
			cout << a[i] << " ";
		}
	}
	cout << '\n';
}

int main()
{

	int tc = 1;
	// cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}