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

ll length(int xa, int xb, int ya, int yb)
{
	return ((xa - xb) * (xa - xb)) + ((ya - yb) * (ya - yb));
}

void solve(int tc = 0)
{
	int xa, ya, xb, yb, xc, yc;
	cin >> xa >> ya >> xb >> yb >> xc >> yc;

	if (length(xa, xb, ya, yb) + length(xa, xc, ya, yc) == length(xb, xc, yb, yc) || length(xa, xb, ya, yb) + length(xb, xc, yb, yc) == length(xa, xc, ya, yc) ||
		length(xb, xc, yb, yc) + length(xa, xc, ya, yc) == length(xa, xb, ya, yb))
	{
		cout << "YES";
	}
	else
		cout << "NO";
}

int main()
{

	int tc = 1;
	// cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}