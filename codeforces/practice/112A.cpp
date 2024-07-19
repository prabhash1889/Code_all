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
	string s1, s2;
	cin >> s1 >> s2;
	for (int i = 0; i < s1.size(); i++)
	{
		s1[i] = toupper(s1[i]);
		s2[i] = toupper(s2[i]);
	}
	if (s1 < s2)
		cout << -1 << '\n';
	else if (s1 == s2)
		cout << 0 << '\n';
	else
		cout << 1 << '\n';
}

int main()
{

	int tc = 1;
	// cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}