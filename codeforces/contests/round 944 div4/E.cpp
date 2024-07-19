#include <bits/stdc++.h>
using namespace std;

void solve(int tc = 1)
{
	int N, K, Q;
	cin >> N >> K >> Q;
	vector<int> A(K), B(K);

	for (auto &a : A)
		cin >> a;

	for (auto &b : B)
		cin >> b;

	A.insert(A.begin(), 0);
	B.insert(B.begin(), 0);

	auto query = [&](int d) -> int64_t
	{
		if (d == 0)
			return 0;

		int index = int(lower_bound(A.begin(), A.end(), d) - A.begin());
		assert(0 < index && index <= K);
		int a_diff = A[index] - A[index - 1];
		int b_diff = B[index] - B[index - 1];
		int d_diff = d - A[index - 1];
		return B[index - 1] + int64_t(b_diff) * d_diff / a_diff;
	};

	for (int q = 0; q < Q; q++)
	{
		int d;
		cin >> d;
		cout << query(d) << (q < Q - 1 ? ' ' : '\n');
	}
}

int main()
{
	ios::sync_with_stdio(false);
	cin.tie(0);

	int tc = 1;
	cin >> tc;
	for (int t = 0; t < tc; t++)
	{
		solve(t);
	}
}
