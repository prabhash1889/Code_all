#include <iostream>
#include <cstdio>
using namespace std;
const int INF = 1e5 + 5;
int a[INF], n, b[INF];
signed main()
{
	ios::sync_with_stdio(false);
	cin >> n;
	int sum = 0, sum1 = 0;
	for (int i = 1; i <= n; i++)
		cin >> a[i];
	for (int i = 1; i <= n; i++)
		cin >> b[i];
	for (int i = 1; i <= n; i++)
		if (a[i] && !b[i])
			sum++;
		else if (!a[i] && b[i])
			sum1++;
	if (sum)
		cout << (sum1 + 1 + sum - 1) / sum << "\n";
	else
		cout << "-1\n";
	return 0;
}