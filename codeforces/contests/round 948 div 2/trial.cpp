#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
// when n! <1e18 and you need raw
long long int nFactorial(int n)
{
	long long int ans = 1;
	for (int i = 2; i <= n; i++)
		ans *= i;
	return ans;
}
int main()
{
	ll a = nFactorial(20);
	cout << a << '\n';
}