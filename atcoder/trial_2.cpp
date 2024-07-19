#undef _GLIBCXX_DEBUG

// https://judge.yosupo.jp/submission/159163
// idk...

#pragma GCC optimize("O3,unroll-loops")
// #pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")
// #pragma GCC target("avx512f")

using i32 = int;
using u32 = unsigned;
using i64 = long long;
using u64 = unsigned long long;

#include <bits/stdc++.h>
using namespace std;

/*
 * Solution based on https://codeforces.com/blog/entry/117947
 * Iterative and in-place version.
 * Uses signed montgomery
 * Optimized to minimize memory usage
 * Trying out 4 radix recursive DFS NTT
 */

// IO from https://judge.yosupo.jp/submission/142782
// #include <sys/mman.h>
#include <sys/stat.h>
#include <cstring>

// Solution based on https://codeforces.com/blog/entry/117947

const int MOD = 998244353;
const long long MOD2 = (long long)MOD * MOD;
const int root = 3;
const int alim = 64; // Bound for using O(n^2) polynomial mult

int modpow(int b, int e)
{
	int ans = 1;
	for (; e; b = (long long)b * b % MOD, e /= 2)
		if (e & 1)
			ans = (long long)ans * b % MOD;
	return ans;
}

const int MODinv = 2 - MOD; // pow(-MOD, -1, 2**32)
inline int m_reduce(long long x)
{
	int m = x * MODinv;
	return (x >> 32) - (((long long)m * MOD) >> 32);
}

const int r2 = modpow(2, 64);
inline int m_transform(int x)
{
	return m_reduce((long long)x * r2);
}

inline int m_add(int x, int y)
{
	int z = x + y;
	return z < 0 ? z + MOD : z - MOD;
}

inline int m_sub(int x, int y)
{
	int z = x - y;
	return z < 0 ? z + MOD : z - MOD;
}

inline int m_mult(int x, int y)
{
	return m_reduce((long long)x * y);
}

vector<int> rt = {1};
vector<int> transformed_rt;
vector<int> transformed_rt2;

template <int a>
void transform(vector<int> &P)
{
	int m = P.size();
	int n = m / a;

	int size = rt.size();
	while (2 * size < n)
	{
		rt.resize(n / 2);
		int r = modpow(root, MOD / (4 * size));
		for (int i = 0; i < size; ++i)
			rt[i + size] = (long long)r * rt[i] % MOD;
		size *= 2;
	}

	// For montgomery
	for (int i = transformed_rt.size(); i < rt.size(); ++i)
	{
		transformed_rt.resize(rt.size());
		transformed_rt[i] = m_transform(rt[i]);
		transformed_rt2.resize(rt.size());
		transformed_rt2[i] = (unsigned int)MODinv * transformed_rt[i];
	}

	// Radix 4 recursive NTT
	auto dfs = [&](auto &&self, int i, int k) -> void
	{
		if (k == 1)
			return;
		int step = k * a;
		int quarter_step = step / 4;

		int R20 = transformed_rt2[2 * i];
		int RR0 = transformed_rt[2 * i];

		int R21 = transformed_rt2[2 * i + 1];
		int RR1 = transformed_rt[2 * i + 1];

		int R2 = transformed_rt2[i];
		int RR = transformed_rt[i];

		int *P1 = &P[i * step];
		int *P2 = P1 + quarter_step;
		int *P3 = P2 + quarter_step;
		int *P4 = P3 + quarter_step;

#pragma GCC ivdep
		for (int j = 0; j < quarter_step; ++j)
		{
			int z0;
			{
				int z = P3[j];
				int m = (unsigned int)R2 * z;
				z0 = ((long long)z * RR - (long long)m * MOD) >> 32;
			}

			int z1;
			{
				int z = P4[j];
				int m = (unsigned int)R2 * z;
				z1 = ((long long)z * RR - (long long)m * MOD) >> 32;
			}

			int sum0 = m_add(P1[j], z0);
			int diff0 = m_sub(P1[j], z0);
			int sum1 = P2[j] + z1;
			int diff1 = P2[j] - z1;

			// [sum0, sum1, diff0, diff1]

			int zz0;
			{
				int z = sum1;
				int m = (unsigned int)R20 * z;
				zz0 = ((long long)z * RR0 - (long long)m * MOD) >> 32;
			}

			int zz1;
			{
				int z = diff1;
				int m = (unsigned int)R21 * z;
				zz1 = ((long long)z * RR1 - (long long)m * MOD) >> 32;
			}

			P1[j] = m_add(sum0, zz0);
			P2[j] = m_sub(sum0, zz0);
			P3[j] = m_add(diff0, zz1);
			P4[j] = m_sub(diff0, zz1);
		}

		self(self, 4 * i + 0, k / 4);
		self(self, 4 * i + 1, k / 4);
		self(self, 4 * i + 2, k / 4);
		self(self, 4 * i + 3, k / 4);
	};

	int k = n;
	while (k >= 4)
		k /= 4;

	if (k == 2)
	{
		int step = n * a;
		int half_step = step / 2;
		for (int j1 = 0; j1 < half_step; ++j1)
		{
			int j2 = j1 + half_step;

			int diff = m_sub(P[j1], P[j2]);
			P[j1] = m_add(P[j1], P[j2]);
			P[j2] = diff;
		}
		k = n / 2;
		dfs(dfs, 0, k);
		dfs(dfs, 1, k);
	}
	else
	{
		k = n;
		dfs(dfs, 0, k);
	}

	for (int i = 0; i < m; ++i)
		if (P[i] < 0)
			P[i] += MOD;
}

template <int a>
void inverse_transform(vector<int> &P)
{
	int m = P.size();
	int n = m / a;
	int n_inv = m_transform(modpow(n, MOD - 2));

	vector<int> rev(n);
	for (int i = 1; i < n; ++i)
	{
		rev[i] = rev[i / 2] / 2 + (i & 1) * n / 2;
	}

	// P = [p * n_inv for p in P]
	for (int i = 0; i < m; ++i)
		P[i] = m_mult(n_inv, P[i]);

	// P = [P[a * rev[i // a] + (i % a)] for i in range(m)]
	for (int i = 1; i < n; ++i)
		if (i < rev[i])
			swap_ranges(P.begin() + a * i, P.begin() + a * i + a, P.begin() + a * rev[i]);

	// P = [P[-a * (i // a) + (i % a)] for i in range(m)]
	for (int i = 1; i < n / 2; ++i)
		swap_ranges(P.begin() + a * i, P.begin() + a * i + a, P.begin() + a * (n - i));

	transform<a>(P);

	// P = [P[a * rev[i // a] + (i % a)] for i in range(m)]
	for (int i = 1; i < n; ++i)
		if (i < rev[i])
			swap_ranges(P.begin() + a * i, P.begin() + a * i + a, P.begin() + a * rev[i]);
}

template <int a>
void fast_polymult_mod(vector<int> &P, vector<int> &Q)
{
	int m = P.size();
	int n = m / a;

	transform<a>(P);
	transform<a>(Q);

	vector<int> &PQ = P;
	for (int i = 0; i < n; ++i)
	{
		vector<unsigned long long> res(2 * a);
		for (int j = 0; j < a; ++j)
		{
			if (j >= 10 && j % 9 == 8)
				for (int k = j; k < j + a - 10; ++k)
					res[k] -= (res[k] >> 63) * 9 * MOD2;
			for (int k = 0; k < a; ++k)
				res[j + k] += (long long)P[i * a + j] * Q[i * a + k];
		}

		int c = rt[i / 2];
		if (i & 1)
			c = MOD - c;
		for (int j = 0; j < a; ++j)
			PQ[i * a + j] = (res[j] + c * (res[j + a] % MOD)) % MOD;
	}

	inverse_transform<a>(PQ);
}

template <size_t... N>
void work(std::index_sequence<N...>, int x, std::vector<int> &a, std::vector<int> &b)
{
	static void (*ptrs[])(std::vector<int> &, std::vector<int> &) = {&fast_polymult_mod<N + 1>...};
	ptrs[x - 1](a, b);
}

void fast_polymult(vector<int> &P, vector<int> &Q)
{
	int m1 = P.size();
	int m2 = Q.size();
	int res_len = m1 + m2 - 1;

	int b = 1;
	while ((alim << b) < res_len)
		++b;
	int a = ((res_len - 1) >> b) + 1;
	int m = a << b;

	P.resize(m);
	Q.resize(m);

	// Call fast_polymult_mod<a>(P, Q);
	work(std::make_index_sequence<alim>{}, a, P, Q);

	P.resize(res_len);
}

#include <bits/stdc++.h>

using namespace std;

#ifdef LOCAL
#include "algo/debug.h"
#else
#define debug(...) 42
#endif

template <typename T>
T inverse(T a, T m)
{
	T u = 0, v = 1;
	while (a != 0)
	{
		T t = m / a;
		m -= t * a;
		swap(a, m);
		u -= t * v;
		swap(u, v);
	}
	assert(m == 1);
	return u;
}

template <typename T>
class Modular
{
public:
	using Type = typename decay<decltype(T::value)>::type;

	constexpr Modular() : value() {}
	template <typename U>
	Modular(const U &x)
	{
		value = normalize(x);
	}

	template <typename U>
	static Type normalize(const U &x)
	{
		Type v;
		if (-mod() <= x && x < mod())
			v = static_cast<Type>(x);
		else
			v = static_cast<Type>(x % mod());
		if (v < 0)
			v += mod();
		return v;
	}

	const Type &operator()() const { return value; }
	template <typename U>
	explicit operator U() const { return static_cast<U>(value); }
	constexpr static Type mod() { return T::value; }

	Modular &operator+=(const Modular &other)
	{
		if ((value += other.value) >= mod())
			value -= mod();
		return *this;
	}
	Modular &operator-=(const Modular &other)
	{
		if ((value -= other.value) < 0)
			value += mod();
		return *this;
	}
	template <typename U>
	Modular &operator+=(const U &other) { return *this += Modular(other); }
	template <typename U>
	Modular &operator-=(const U &other) { return *this -= Modular(other); }
	Modular &operator++() { return *this += 1; }
	Modular &operator--() { return *this -= 1; }
	Modular operator++(int)
	{
		Modular result(*this);
		*this += 1;
		return result;
	}
	Modular operator--(int)
	{
		Modular result(*this);
		*this -= 1;
		return result;
	}
	Modular operator-() const { return Modular(-value); }

	template <typename U = T>
	typename enable_if<is_same<typename Modular<U>::Type, int>::value, Modular>::type &operator*=(const Modular &rhs)
	{
		value = normalize(static_cast<int64_t>(value) * static_cast<int64_t>(rhs.value));
		return *this;
	}
	template <typename U = T>
	typename enable_if<is_same<typename Modular<U>::Type, long long>::value, Modular>::type &operator*=(const Modular &rhs)
	{
		long long q = static_cast<long long>(static_cast<long double>(value) * rhs.value / mod());
		value = normalize(value * rhs.value - q * mod());
		return *this;
	}
	// template <typename U = T>
	// typename enable_if<!is_integral<typename Modular<U>::Type>::value, Modular>::type& operator*=(const Modular& rhs) {
	//   value = normalize(value * rhs.value);
	//   return *this;
	// }

	Modular &operator/=(const Modular &other) { return *this *= Modular(inverse(other.value, mod())); }

	friend const Type &abs(const Modular &x) { return x.value; }

	template <typename U>
	friend bool operator==(const Modular<U> &lhs, const Modular<U> &rhs);

	template <typename U>
	friend bool operator<(const Modular<U> &lhs, const Modular<U> &rhs);

	template <typename V, typename U>
	friend V &operator>>(V &stream, Modular<U> &number);

private:
	Type value;
};

template <typename T>
bool operator==(const Modular<T> &lhs, const Modular<T> &rhs) { return lhs.value == rhs.value; }
template <typename T, typename U>
bool operator==(const Modular<T> &lhs, U rhs) { return lhs == Modular<T>(rhs); }
template <typename T, typename U>
bool operator==(U lhs, const Modular<T> &rhs) { return Modular<T>(lhs) == rhs; }

template <typename T>
bool operator!=(const Modular<T> &lhs, const Modular<T> &rhs) { return !(lhs == rhs); }
template <typename T, typename U>
bool operator!=(const Modular<T> &lhs, U rhs) { return !(lhs == rhs); }
template <typename T, typename U>
bool operator!=(U lhs, const Modular<T> &rhs) { return !(lhs == rhs); }

template <typename T>
bool operator<(const Modular<T> &lhs, const Modular<T> &rhs) { return lhs.value < rhs.value; }

template <typename T>
Modular<T> operator+(const Modular<T> &lhs, const Modular<T> &rhs) { return Modular<T>(lhs) += rhs; }
template <typename T, typename U>
Modular<T> operator+(const Modular<T> &lhs, U rhs) { return Modular<T>(lhs) += rhs; }
template <typename T, typename U>
Modular<T> operator+(U lhs, const Modular<T> &rhs) { return Modular<T>(lhs) += rhs; }

template <typename T>
Modular<T> operator-(const Modular<T> &lhs, const Modular<T> &rhs) { return Modular<T>(lhs) -= rhs; }
template <typename T, typename U>
Modular<T> operator-(const Modular<T> &lhs, U rhs) { return Modular<T>(lhs) -= rhs; }
template <typename T, typename U>
Modular<T> operator-(U lhs, const Modular<T> &rhs) { return Modular<T>(lhs) -= rhs; }

template <typename T>
Modular<T> operator*(const Modular<T> &lhs, const Modular<T> &rhs) { return Modular<T>(lhs) *= rhs; }
template <typename T, typename U>
Modular<T> operator*(const Modular<T> &lhs, U rhs) { return Modular<T>(lhs) *= rhs; }
template <typename T, typename U>
Modular<T> operator*(U lhs, const Modular<T> &rhs) { return Modular<T>(lhs) *= rhs; }

template <typename T>
Modular<T> operator/(const Modular<T> &lhs, const Modular<T> &rhs) { return Modular<T>(lhs) /= rhs; }
template <typename T, typename U>
Modular<T> operator/(const Modular<T> &lhs, U rhs) { return Modular<T>(lhs) /= rhs; }
template <typename T, typename U>
Modular<T> operator/(U lhs, const Modular<T> &rhs) { return Modular<T>(lhs) /= rhs; }

template <typename T, typename U>
Modular<T> power(const Modular<T> &a, const U &b)
{
	assert(b >= 0);
	Modular<T> x = a, res = 1;
	U p = b;
	while (p > 0)
	{
		if (p & 1)
			res *= x;
		x *= x;
		p >>= 1;
	}
	return res;
}

template <typename T>
bool IsZero(const Modular<T> &number)
{
	return number() == 0;
}

template <typename T>
string to_string(const Modular<T> &number)
{
	return to_string(number());
}

// U == std::ostream? but done this way because of fastoutput
template <typename U, typename T>
U &operator<<(U &stream, const Modular<T> &number)
{
	return stream << number();
}

// U == std::istream? but done this way because of fastinput
template <typename U, typename T>
U &operator>>(U &stream, Modular<T> &number)
{
	typename common_type<typename Modular<T>::Type, long long>::type x;
	stream >> x;
	number.value = Modular<T>::normalize(x);
	return stream;
}

/*
using ModType = int;

struct VarMod { static ModType value; };
ModType VarMod::value;
ModType& md = VarMod::value;
using Mint = Modular<VarMod>;
*/

constexpr int md = 998244353;
using Mint = Modular<std::integral_constant<decay<decltype(md)>::type, md>>;

vector<Mint> fact(1, 1);
vector<Mint> inv_fact(1, 1);

Mint C(int n, int k)
{
	if (k < 0 || k > n)
	{
		return 0;
	}
	while ((int)fact.size() < n + 1)
	{
		fact.push_back(fact.back() * (int)fact.size());
		inv_fact.push_back(1 / fact.back());
	}
	return fact[n] * inv_fact[k] * inv_fact[n - k];
}

template <typename T>
class NTT
{
public:
	using Type = typename decay<decltype(T::value)>::type;

	static Type md;
	static Modular<T> root;
	static int base;
	static int max_base;
	static vector<Modular<T>> roots;
	static vector<int> rev;

	static void clear()
	{
		root = 0;
		base = 0;
		max_base = 0;
		roots.clear();
		rev.clear();
	}

	static void init()
	{
		md = T::value;
		assert(md >= 3 && md % 2 == 1);
		auto tmp = md - 1;
		max_base = 0;
		while (tmp % 2 == 0)
		{
			tmp /= 2;
			max_base++;
		}
		root = 2;
		while (power(root, (md - 1) >> 1) == 1)
		{
			root++;
		}
		assert(power(root, md - 1) == 1);
		root = power(root, (md - 1) >> max_base);
		base = 1;
		rev = {0, 1};
		roots = {0, 1};
	}

	static void ensure_base(int nbase)
	{
		if (md != T::value)
		{
			clear();
		}
		if (roots.empty())
		{
			init();
		}
		if (nbase <= base)
		{
			return;
		}
		assert(nbase <= max_base);
		rev.resize(1 << nbase);
		for (int i = 0; i < (1 << nbase); i++)
		{
			rev[i] = (rev[i >> 1] >> 1) + ((i & 1) << (nbase - 1));
		}
		roots.resize(1 << nbase);
		while (base < nbase)
		{
			Modular<T> z = power(root, 1 << (max_base - 1 - base));
			for (int i = 1 << (base - 1); i < (1 << base); i++)
			{
				roots[i << 1] = roots[i];
				roots[(i << 1) + 1] = roots[i] * z;
			}
			base++;
		}
	}

	static void fft(vector<Modular<T>> &a)
	{
		int n = (int)a.size();
		assert((n & (n - 1)) == 0);
		int zeros = __builtin_ctz(n);
		ensure_base(zeros);
		int shift = base - zeros;
		for (int i = 0; i < n; i++)
		{
			if (i < (rev[i] >> shift))
			{
				swap(a[i], a[rev[i] >> shift]);
			}
		}
		for (int k = 1; k < n; k <<= 1)
		{
			for (int i = 0; i < n; i += 2 * k)
			{
				for (int j = 0; j < k; j++)
				{
					Modular<T> x = a[i + j];
					Modular<T> y = a[i + j + k] * roots[j + k];
					a[i + j] = x + y;
					a[i + j + k] = x - y;
				}
			}
		}
	}

	static vector<Modular<T>> multiply(vector<Modular<T>> a, vector<Modular<T>> b)
	{
		if (a.empty() || b.empty())
		{
			return {};
		}
		int eq = (a == b);
		int need = (int)a.size() + (int)b.size() - 1;
		int nbase = 0;
		while ((1 << nbase) < need)
			nbase++;
		ensure_base(nbase);
		int sz = 1 << nbase;
		a.resize(sz);
		b.resize(sz);
		fft(a);
		if (eq)
			b = a;
		else
			fft(b);
		Modular<T> inv_sz = 1 / static_cast<Modular<T>>(sz);
		for (int i = 0; i < sz; i++)
		{
			a[i] *= b[i] * inv_sz;
		}
		reverse(a.begin() + 1, a.end());
		fft(a);
		a.resize(need);
		return a;
	}
};

template <typename T>
typename NTT<T>::Type NTT<T>::md;
template <typename T>
Modular<T> NTT<T>::root;
template <typename T>
int NTT<T>::base;
template <typename T>
int NTT<T>::max_base;
template <typename T>
vector<Modular<T>> NTT<T>::roots;
template <typename T>
vector<int> NTT<T>::rev;

template <typename T>
vector<Modular<T>> inverse(const vector<Modular<T>> &a)
{
	assert(!a.empty());
	int n = (int)a.size();
	vector<Modular<T>> b = {1 / a[0]};
	while ((int)b.size() < n)
	{
		vector<Modular<T>> x(a.begin(), a.begin() + min(a.size(), b.size() << 1));
		x.resize(b.size() << 1);
		b.resize(b.size() << 1);
		vector<Modular<T>> c = b;
		NTT<T>::fft(c);
		NTT<T>::fft(x);
		Modular<T> inv = 1 / static_cast<Modular<T>>((int)x.size());
		for (int i = 0; i < (int)x.size(); i++)
		{
			x[i] *= c[i] * inv;
		}
		reverse(x.begin() + 1, x.end());
		NTT<T>::fft(x);
		rotate(x.begin(), x.begin() + (x.size() >> 1), x.end());
		fill(x.begin() + (x.size() >> 1), x.end(), 0);
		NTT<T>::fft(x);
		for (int i = 0; i < (int)x.size(); i++)
		{
			x[i] *= c[i] * inv;
		}
		reverse(x.begin() + 1, x.end());
		NTT<T>::fft(x);
		for (int i = 0; i < ((int)x.size() >> 1); i++)
		{
			b[i + ((int)x.size() >> 1)] = -x[i];
		}
	}
	b.resize(n);
	return b;
}

template <typename T>
vector<Modular<T>> inverse_old(vector<Modular<T>> a)
{
	assert(!a.empty());
	int n = (int)a.size();
	if (n == 1)
	{
		return {1 / a[0]};
	}
	int m = (n + 1) >> 1;
	vector<Modular<T>> b = inverse(vector<Modular<T>>(a.begin(), a.begin() + m));
	int need = n << 1;
	int nbase = 0;
	while ((1 << nbase) < need)
	{
		++nbase;
	}
	NTT<T>::ensure_base(nbase);
	int size = 1 << nbase;
	a.resize(size);
	b.resize(size);
	NTT<T>::fft(a);
	NTT<T>::fft(b);
	Modular<T> inv = 1 / static_cast<Modular<T>>(size);
	for (int i = 0; i < size; ++i)
	{
		a[i] = (2 - a[i] * b[i]) * b[i] * inv;
	}
	reverse(a.begin() + 1, a.end());
	NTT<T>::fft(a);
	a.resize(n);
	return a;
}

template <typename T>
vector<Modular<T>> operator*(const vector<Modular<T>> &a, const vector<Modular<T>> &b)
{
	if (a.empty() || b.empty())
	{
		return {};
	}
	if (min(a.size(), b.size()) < 150)
	{
		vector<Modular<T>> c(a.size() + b.size() - 1, 0);
		for (int i = 0; i < (int)a.size(); i++)
		{
			for (int j = 0; j < (int)b.size(); j++)
			{
				c[i + j] += a[i] * b[j];
			}
		}
		return c;
	}
	return NTT<T>::multiply(a, b);
}

template <typename T>
vector<Modular<T>> &operator*=(vector<Modular<T>> &a, const vector<Modular<T>> &b)
{
	return a = a * b;
}

int main()
{
	ios::sync_with_stdio(false);
	cin.tie(0);
	int tt;
	cin >> tt;
	while (tt--)
	{
		int n, m, b0;
		cin >> n >> m >> b0;
		if (b0 >= m || b0 >= n)
		{
			cout << power(Mint(m), n) << '\n';
			continue;
		}
		if (m == 1)
		{
			cout << 0 << '\n';
			continue;
		}
		auto orig_m = m;
		auto inv = (orig_m == 1 ? 0 : 1 / Mint(orig_m - 1));
		vector<Mint> pw(n + 1);
		vector<Mint> pw_m(n + 1);
		vector<Mint> inv_pw(n + 1);
		pw[0] = 1;
		pw_m[0] = 1;
		inv_pw[0] = 1;
		for (int i = 1; i <= n; i++)
		{
			pw[i] = pw[i - 1] * (orig_m - 1);
			pw_m[i] = pw_m[i - 1] * orig_m;
			inv_pw[i] = inv_pw[i - 1] * inv;
		}
		m = min(m, b0 + n + 1);
		Mint ans = 0;
		int par = b0 & 1;
		vector<Mint> dp((m + 1 - par) >> 1);
		dp[b0 >> 1] = 1;
		while (n > 0)
		{
			int steps = min(n, m);
			int new_par = par ^ (steps & 1);
			vector<int> p(steps + 1);
			for (int i = 0; i <= steps; i++)
			{
				p[i] = (C(steps, i) * pw[i])();
			}
			vector<int> dp_int(dp.size());
			for (int i = 0; i < int(dp.size()); i++)
			{
				dp_int[i] = dp[i]();
			}
			Mint collected = 0;
			fast_polymult(dp_int, p);
			// vector<Mint> dp_int(dp_int_int.size());
			// for (int i = 0; i < int(dp_int.size()); i++) {
			//   dp_int[i] = dp_int_int[i];
			// }
			vector<Mint> new_dp((m + 1 - new_par) >> 1);
			for (int i = 0; i < int(dp_int.size()); i++)
			{
				int at = (2 * i + par) - steps;
				if (at < 0)
				{
					continue;
				}
				if (at >= m)
				{
					collected += dp_int[i];
					continue;
				}
				new_dp[at >> 1] += dp_int[i];
				{
					int fake = 2 * m - at;
					if (fake + steps >= 0 && ((fake + steps) >> 1) < int(dp_int.size()))
					{
						Mint sub = dp_int[(fake + steps) >> 1];
						sub *= inv_pw[(fake - at) >> 1];
						collected += sub;
						new_dp[at / 2] -= sub;
					}
				}
				{
					int fake = -2 - at;
					if (fake + steps >= 0 && ((fake + steps) >> 1) < int(dp_int.size()))
					{
						Mint sub = dp_int[(fake + steps) >> 1];
						sub *= pw[(at - fake) >> 1];
						new_dp[at / 2] -= sub;
					}
				}
			}
			swap(dp, new_dp);
			ans += collected * pw_m[n - steps];
			n -= steps;
			par = new_par;
		}
		ans += accumulate(dp.begin(), dp.end(), Mint(0));
		cout << ans << '\n';
	}
	debug(clock());
	return 0;
}
