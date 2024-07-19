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
	for (int i = 0; i < k; i++)
	{
		if (n % 10 == 0)
		{
			n /= 10;
		}
		else if (n % 10 != 10)
		{
			n--;
		}
	}
	cout << n << '\n';
}

int main()
{

	int tc = 1;
	// cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}