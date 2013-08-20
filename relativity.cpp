#include <assert.h>

/*
Starting again from scratch now that I've gone through the ADM chapter of Gravitation and Numerical Relativity 
*/

template<int dim_, typename type_>
struct vec {
	enum { dim = dim_ };
	typedef type_ type;

	type v[dim];

	vec() {
		for (int i = 0; i < dim; ++i) {
			v[i] = type();
		}
	}

	vec(const type &x) {
		assert(dim >= 1);
		v[0] = x;
	}

	vec(const type &x, const type &y) {
		assert(dim >= 2);
		v[0] = x;
		v[1] = y;
	}

	vec(const type &x, const type &y, const type &z) {
		assert(dim >= 3);
		v[0] = x;
		v[1] = y;
		v[2] = z;
	}

	type &operator()(int i) { return v[i]; }

	type volume() const {
		type vol = type(1);
		for (int i = 0; i < dim; ++i) {
			vol *= v[i];
		}
		return vol;
	}

	static type dot(const vec &a, const vec &b) {
		type d = type(0);
		for (int i = 0; i < dim; ++i) {
			d += a(i) * b(i);
		}
		return d;
	}
};

/*
symmat(i,j) == symmat(j,k)
*/
template<int dim_, typename type_>
struct symmat {
	enum { dim = dim_ };
	enum { size = dim * (dim+1) / 2 };
	typedef type_ type;

	type v[size];

	/*
	initialize to identity or zero?
	identity
	*/
	symmat() {
		for (int i = 0; i < size; ++i) {
			v[i] = type();
		}
		for (int i = 0; i < dim; ++i) {
			(*this)(i,i) = type(1);
		}
	}

	/*
	math-index: i is the row, j is the column
	row-major: i is nested inner-most
	upper triangular: i <= j
	*/
	static int index(int i, int j) {
		if (i > j) return index(j,i);
		//j == 0: return 0
		//j == 1: return 1 + i
		//j == 2: return 1 + 2 + i
		//j == j: return j * (j+1)/2 + i
		return j * (j + 1) / 2 + i;
	}
	
	type &operator()(int i, int j) { return v[index(i,j)];}
};

//rank is templated, but dim is not as it varies per-rank
//so this is dynamically-sized
template<int rank_, typename type_>
struct tensor {
	enum { rank = rank_ };
	typedef type_ type;

	typedef ::vec<rank,int> veci;

	type *v;
	veci size;
	
	//cached for quick access by dot with index vector
	//step[0] = 1, step[1] = size[0], step[j] = product(i=1,j-1) size[i]
	veci step;

	tensor(const veci &size_) : size(size_) {
		v = new type[size.volume()];
		step[0] = 1;
		for (int i = 1; i < rank; ++i) {
			step[i] = step[i-1] * size[i-1];
		}
	}

	struct iterator {
		veci index;
		
		iterator() {}
		iterator(const iterator &iter) : index(iter.index) {}
		
		bool operator==(const iterator &b) const {
			return index == b.index;
		}
		iterator &operator++(int) {
			for (int i = 0; i < rank-1; ++i) {	//allow the last index to overflow for sake of comparing it to end
				++index[i];
				if (index[i] < size[i]) break;
				index[i] = 0;
			}
			return *this;
		}
	};

	iterator begin() {
		iterator i;
		--i.index[0];
		return i;
	}
	iterator end() {
		iterator i;
		i.index[rank-1] = size[rank-1];
		return i;
	}

	type &operator()(const veci &index) {
		return v[veci::dot(index, step)];
	}
};

template<int dim_, typename real_>
struct Cell {
	typedef real_ real;
	enum { dim = dim_ };

	typedef ::vec<dim,real> vec;
	typedef ::symmat<dim,real> symmat;

	//lapse
	real alpha;

	//shift
	//beta^i
	//beta^t = 0
	vec beta_u;

	//spatial metric
	//gamma_ij
	//gamma_it = gamma_tj = 0
	symmat gamma_ll;
	
	//extrinsic curvature
	//K_ij
	//K_it = K_tj = 0
	symmat K_ll;
};

template<int dim_, typename real_>
struct Simulation {
	enum { dim = dim_ };
	typedef real_ real;

	typedef ::vec<dim,real> vec;
	typedef ::vec<dim,int> veci;
	typedef ::symmat<dim,int> symmat;
	typedef ::Cell<dim,real> Cell;
	typedef ::tensor<dim,Cell> CellTensor;

	CellTensor readCells, writeCells;

	Simulation(const veci &size) 
	: readCells(size), writeCells(size)
	{
	}

	void update() {
		typename CellTensor::iterator iter = writeCells.begin();
		for (; iter != writeCells.end(); ++iter) {
			Cell &writeCell = *iter;
			const Cell &readCell = readCells(iter.index);
			
			const real &alpha = readCell.alpha;
			const vec &beta_u = readCell.beta_u;
			const symmat &K_ll = readCell.K_ll;

			//beta_u = g_uv beta^v
			vec beta_l;

			//diff_i beta_j = 
			symmat D_beta_ll;
			
			Cell partial_t;
		
			//partial_t gamma_ij = -2 alpha K_ij + D_i beta_j + D_j beta_i
			for (int j = 0; j < dim; ++j) {
				for (int i = 0; i <= j; ++i) {
					partial_t.gamma_ll(i,j) = -2 * alpha * K_ll(i,j) + D_beta_ll(i,j) + D_beta_ll(j,i);
				}
			}
			
			//partial_t K_ij = alpha (R_ij - 2 K_ik K^k_j + K K_ij) - D_i D_j alpha - 8 pi alpha (S_ij - 1/2 gamma_ij (S - rho)) + beta^k partial_k K_ij + K_ik partial_j beta^k + K_kj partial_i beta^k

		}
	}
};

int main() {
	Simulation<2,double> sim(vec<2,int>(10, 10));

}

