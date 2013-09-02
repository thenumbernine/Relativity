#include <iostream>
using namespace std;

#include "../tensor.h"

//I know, it's not a legitimate unit test.
// but maybe some day it will become one ...

int test_tensors() {
	typedef double real;
	typedef tensor<real,upper<3>> vector;
	vector v;
	cout << "size " << v.size() << endl;
	cout << "rank " << v.rank << endl;
	for (int i = 0; i < 3; ++i) {
		cout << "v^" << i << " = " << v(i) << endl;
	}

	typedef tensor<real,lower<3>> oneform;
	oneform w;
	cout << "size " << w.size() << endl;
	cout << "rank " << w.rank << endl;
	for (int i = 0; i < 3; ++i) {
		cout << "w_" << i << " = " << w(i) << endl;
	}

	typedef tensor<real,symmetric<lower<3>,lower<3>>> metric;
	metric g;
	for (int i = 0; i < 3; ++i) {
		g.body(i,i) = 1;
	}
	cout << "size " << g.size() << endl;
	cout << "rank " << g.rank << endl;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			cout << "g_" << i << j << " = " << g.body(i,j) << endl;
		}
	}

	typedef tensor<real,upper<3>,upper<3>> matrix;
	matrix h;
	int index = 0;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			h(i,j) = ++index;
		}
	}
	cout << "size " << h.size() << endl;
	cout << "rank " << h.rank << endl;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			cout << "h^" << i << j << " = " << h(i,j) << endl;
		}
		matrix h2;
		h2(i) = h(i);
		cout << "h^"<<i<<" = "<< h2 << endl;
	}

	//subtensor access
	tensor<real, upper<3>, upper<3>, upper<3>> ta;
	for (tensor<real, upper<3>, upper<3>, upper<3>>::iterator i = ta.begin(); i != ta.end(); ++i) *i = 1.;
	cout << "ta = " << ta << endl;
	tensor<real, upper<3>, upper<3>> tb;
	for (tensor<real, upper<3>, upper<3>>::iterator i = tb.begin(); i != tb.end(); ++i) *i = 2.;
	cout << "tb = " << tb << endl;
	ta(0) = tb;
	cout << "ta = " << ta << endl;
	tensor<real, upper<3>> tc;
	for (tensor<real, upper<3>>::iterator i = tc.begin(); i != tc.end(); ++i) *i = 3.;
	cout << "tc = " << tc << endl;
	ta(0,0) = tc;
	cout << "ta = " << ta << endl;


#if 0	//not yet working
	typedef tensor<real, 
		antisymmetric<lower<2>, lower<2>>,
		antisymmetric<lower<2>, lower<2>>
	> RiemannTensor;
	typedef RiemannTensor::TensorStats RiemannTensorStats;
	typedef RiemannTensorStats::InnerType RiemannTensorStatsInnerType;
	typedef RiemannTensorStatsInnerType::BodyType RiemannTensorStatsInnerBodyType;
	RiemannTensor r;
	cout << "size " << r.size() << endl;
	cout << "rank " << r.rank << endl;
	int e = 0;
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < i; ++j) {
			for (int k = 0; k < 2; ++k) {
				for (int l = 0; l < k; ++l) {
					((RiemannTensorStatsInnerBodyType&)(r.body(i,j)))(k,l) = ++e;
				}
			}
		}
	}
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			for (int k = 0; k < 2; ++k) {
				for (int l = 0; l < 2; ++l) {
					cout << "r_" << i << j << k << l << " = " << ((RiemannTensorStatsInnerBodyType&)(r.body(i,j)))(k,l) << endl;
				}
			}
		}
	}
#endif
}

int main() {
	test_tensors();
}

