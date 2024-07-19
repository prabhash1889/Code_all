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
	vector<int> v;
	for (int i = 0; i < s.size(); i++)
	{
		if (s[i] == '+')
		{
			continue;
		}
		else
		{
			v.push_back(s[i] - '0');
		}
	}
	// for (int i = 0; i < v.size(); i++)
	// 	cout << v[i] << " ";
	sort(v.begin(), v.end());
	for (int i = 0; i < v.size() - 1; i++)
	{
		cout << v[i];
		cout << '+';
	}
	cout << v[v.size() - 1] << '\n';
}

int main()
{

	int tc = 1;
	// cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}