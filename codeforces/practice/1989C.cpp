#include <bits/stdc++.h>
using namespace std;
#include <iostream>

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
	int a[n], b[n];
	for (int i = 0; i < n; i++)
	{
		cin >> a[i];
	}

	for (int i = 0; i < n; i++)
		cin >> b[i];
	int v1 = 0, v2 = 0, sp = 0, sm = 0;
	for (int i = 0; i < n; i++)
	{
		if (a[i] > b[i])
			v1 += a[i];
		else if (b[i] > a[i])
			v2 += b[i];
		else if (a[i] == 1)
			sp++;
		else if (a[i] == -1)
			sm++;
	}
	while (sp--)
	{
		if (v1 < v2)
			v1++;
		else
			v2++;
	}
	while (sm--)
	{
		if (v1 > v2)
			v1--;
		else
			v2--;
	}
	cout << min(v1, v2) << '\n';
}

int main()
{

	int tc = 1;
	cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}