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
	int x;
	cin >> x;
	vector<int> a;

	while (x > 0)
	{
		if ((x % 2) == 0)
		{
			a.push_back(0);
		}
		else
		{
			if ((x - 1) / 2 % 2 == 0)
			{
				a.push_back(1);
				x -= 1;
			}
			else
			{
				a.push_back(-1);
				x += 1;
			}
		}
		x /= 2;
	}
	cout << a.size() << '\n';
	for (int i = 0; i < int(a.size()); i++)
	{
		cout << a[i] << " \n"[i == int(a.size()) - 1];
	}
}

int main()
{

	int tc = 1;
	cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}