#include<bits/stdc++.h>
using namespace std;
#include<iostream>

typedef long long ll;
typedef vector<int> vi;
ll n, m, k, q, l, r, x, y, z;
const ll template_array_size = 1e6 + 585;
ll a[template_array_size];
ll b[template_array_size];
ll c[template_array_size];
string s, t;
ll ans = 0;


bool t(vector<int> a){
	for(int i=0;i<a.size();i++){
		if(a[i]==i){
			return true;
		}
	}
	return false;
}

void solve(int tc = 0) {
	int n;
	cin>>n;
	vector<int> a;
	vector<int> b;

	for(int i=0;i<n;i++){
		int x;
		cin>>x;
		a.push_back(x);
		if(x==i){
			b.push_back(x);
		}
	}
	int ans=0;
	for(int i=0;i<100;i++){
		
	}
	

}

int main() {
	
 
	int tc = 1;
	// cin >> tc;
	for (int t = 0; t < tc; t++) solve(t);
	
	
} 