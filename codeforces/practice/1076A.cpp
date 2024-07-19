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
	string s;
	cin >> n >> s;
	int pos = n - 1;
	for (int i = 0; i < n - 1; ++i)
	{
		if (s[i] > s[i + 1])
		{
			pos = i;
			break;
		}
	}
	cout << s.substr(0, pos) + s.substr(pos + 1) << endl;
}

int main()
{

	int tc = 1;
	// cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}