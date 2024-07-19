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
	long long n;

	cin >> n;
	int g = 0;
	while (n > 0)
	{
		g += (n % 10 == 4 || n % 10 == 7);
		n /= 10;
	}

	cout << ((g == 4 || g == 7) ? "YES" : "NO");
}

int main()
{

	int tc = 1;
	// cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}