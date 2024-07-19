#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <tuple>
#include <queue>
#include <numeric>
#define pii pair<int, int>
#define pll pair<long long, long long>
#define piii pair<int, pii>
#define plll pair<long long, pll>
#define tiii tuple<int, int, int>
#define tiiii tuple<int, int, int, int>
#define ff first
#define ss second
#define ee ss.ff
#define rr ss.ss
const int INF = (int)1e9 + 7;

using namespace std;

int main()
{
	ios_base::sync_with_stdio(false);
	cin.tie(nullptr);

	int T;
	cin >> T;
	while (T--)
	{
		int n, p;
		cin >> n >> p;
		pii A[n];
		for (int i = 0; i < n; ++i)
			cin >> A[i].ff;
		for (int i = 0; i < n; ++i)
			cin >> A[i].ss;
		sort(A, A + n, [&](pii x, pii y)
			 { return x.ss < y.ss; });

		long long ans = 0;

		queue<int> Q;
		for (int i = 0; i < n; ++i)
		{
			if (Q.empty())
				ans += p;
			else
				ans += Q.front(), Q.pop();
			if (A[i].ss < p)
				for (int j = 0; j < A[i].ff; ++j)
				{
					if (Q.size() >= n)
						break;
					Q.push(A[i].ss);
				}
		}

		cout << ans << '\n';
	}
}