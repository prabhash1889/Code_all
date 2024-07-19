#include <bits/stdc++.h>
#include <iostream>
using namespace std;

int main()
{
	int t;
	cin >> t;
	while (t--)
	{
		string a, b;
		cin >> a >> b;
		int n = a.size(), m = b.size();
		int ans = n + m;
		for (int i = 0; i < m; ++i)
		{
			int j = i;
			for (int k = 0; k < a.size(); k++)
			{
				if (j < m && a[k] == b[j])
					++j;
			}
			ans = min(ans, n + m - (j - i));
		}
		cout << ans << '\n';
	}
}