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
	int n, p;
	cin >> n >> p;
	pair<int, int> A[n];
	for (int i = 0; i < n; i++)
	{
		cin >> A[i].first;
	}
	for (int i = 0; i < n; i++)
		cin >> A[i].second;
	sort(A, A + n, [&](pair<int, int> x, pair<int, int> y)
		 { return x.second < y.second; });
	ll ans = 0;
	queue<int> Q;
	for (int i = 0; i < n; i++)
	{
		if (Q.empty())
			ans += p;
		else
			ans += Q.front();
		Q.pop();
		if (A[i].second < p)
		{
			for (int j = 0; j < A[i].first; j++)
			{
				if (Q.size() >= n)
				{
					break;
				}
				Q.push(A[i].second);
			}
		}
	}
	cout << ans << '\n';
}

int main()
{

	int tc = 1;
	cin >> tc;
	for (int t = 0; t < tc; t++)
		solve(t);
}