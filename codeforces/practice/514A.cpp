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
	string s;
	cin >> s;
	for (int i = 0; i < s.size(); i++)
	{
		if (!(i == 0 && s[i] == '9') && s[i] - '0' > 9 - (s[i] - '0'))
		{
			s[i] = 9 - (s[i] - '0') + '0';
		}
	}

	cout << s << '\n';
}

int main()
{

	int tc = 1;
	// cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}