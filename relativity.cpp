#include <assert.h>
#include <algorithm>	//min/max

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
		for (int i = 0; i < dim; ++i) {
			v[i] = x;
		}
	}

	vec(const type &x, const type &y) {
		assert(dim >= 2);
		v[0] = x;
		v[1] = y;
		for (int i = 2; i < dim; ++i) v[i] = type();
	}

	vec(const type &x, const type &y, const type &z) {
		assert(dim >= 3);
		v[0] = x;
		v[1] = y;
		v[2] = z;
		for (int i = 3; i < dim; ++i) v[i] = type();
	}

	//index access
	type &operator()(int i) { return v[i]; }
	const type &operator()(int i) const { return v[i]; }

	//cast operation
	template<typename type2>
	operator vec<dim,type2>() const {
		vec<dim,type2> result;
		for (int i = 0; i < dim; ++i) {
			result(i) = (type2)(*this)(i);
		}
		return result;
	}

	//product of elements / flat-space volume operator
	type volume() const {
		type vol = type(1);
		for (int i = 0; i < dim; ++i) {
			vol *= v[i];
		}
		return vol;
	}

	//inner product / flat-space dot product
	static type dot(const vec &a, const vec &b) {
		type d = type(0);
		for (int i = 0; i < dim; ++i) {
			d += a(i) * b(i);
		}
		return d;
	}

	//bounds
	static vec clamp(const vec &a, const vec &min, const vec &max) {
		vec b;
		for (int i = 0; i < dim; ++i) {
			b(i) = std::max(min(i), std::min(max(i), a(i)));
		}
		return b;
	}

	//math operations

	bool operator==(const vec &b) const {
		const vec &a = *this;
		for (int i = 0; i < dim; ++i) {
			if (a(i) != b(i)) return false;
		}
		return true;
	}

	bool operator!=(const vec &b) const { return ! this->operator==(b); }

	vec operator+(const vec &b) const {
		const vec &a = *this;
		vec c;
		for (int i = 0; i < dim; ++i) {
			c(i) = a(i) + b(i);
		}
		return c;
	}
	
	vec operator-(const vec &b) const {
		const vec &a = *this;
		vec c;
		for (int i = 0; i < dim; ++i) {
			c(i) = a(i) - b(i);
		}
		return c;
	}

	vec operator*(const type &b) const {
		const vec &a = *this;
		vec c;
		for (int i = 0; i < dim; ++i) {
			c(i) = a(i) * b; 
		}
		return c;
	}

	vec operator*(const vec &b) const {
		const vec &a = *this;
		vec c;
		for (int i = 0; i < dim; ++i) {
			c(i) = a(i) * b(i); 
		}
		return c;
	}
	
	vec operator/(const type &b) const {
		const vec &a = *this;
		vec c;
		for (int i = 0; i < dim; ++i) {
			c(i) = a(i) / b; 
		}
		return c;
	}

	vec operator/(const vec &b) const {
		const vec &a = *this;
		vec c;
		for (int i = 0; i < dim; ++i) {
			c(i) = a(i) / b(i); 
		}
		return c;
	}
};

/*
child being the whatever curious whatever thing that returns its child
child is used for
	vector dereference index calculation
	operator return type
*/
template<int dim_, typename type_, typename child, int size_>
struct generic_dense_matrix {
	enum { dim = dim_ };
	enum { size = size_ };
	typedef type_ type;

	type v[size];

	/*
	initialize to identity or zero?
	identity
	*/
	generic_dense_matrix() {
		for (int i = 0; i < size; ++i) {
			v[i] = type();
		}
		for (int i = 0; i < dim; ++i) {
			(*this)(i,i) = type(1);
		}
	}

	generic_dense_matrix(const child &a) {
		for (int i = 0; i < size; ++i) {
			v[i] = a.v[i];
		}
	}

	//index access
	type &operator()(int i, int j) { return v[child::index(i,j)];}
	const type &operator()(int i, int j) const { return v[child::index(i,j)];}

	//math operations

	bool operator==(const child &b) const {
		const child &a = *this;
		for (int i = 0; i < size; ++i) {
			if (a.v[i] != b.v[i]) return false;
		}
		return true;
	}

	bool operator!=(const child &b) const { return ! this->operator==(b); }

	child operator+(const child &b) const {
		const generic_dense_matrix &a = *this;
		child c;
		for (int i = 0; i < size; ++i) {
			c.v[i] = a.v[i] + b.v[i];
		}
		return c;
	}
	
	child operator-(const child &b) const {
		const generic_dense_matrix &a = *this;
		child c;
		for (int i = 0; i < size; ++i) {
			c.v[i] = a.v[i] - b.v[i];
		}
		return c;
	}

	child operator*(const type &b) const {
		const generic_dense_matrix &a = *this;
		child c;
		for (int i = 0; i < dim; ++i) {
			c.v[i] = a.v[i] * b;
		}
		return c;
	}

	child operator/(const type &b) const {
		const generic_dense_matrix &a = *this;
		child c;
		for (int i = 0; i < dim; ++i) {
			c.v[i] = a.v[i] / b;
		}
		return c;
	}
};

template<int dim_, typename type_>
struct matrix : public generic_dense_matrix<dim_, type_, matrix<dim_, type_>, dim_ * dim_ > {
	typedef generic_dense_matrix<dim_, type_, matrix<dim_, type_>, dim_ * dim_> parent;
	
	enum { dim = parent::dim };
	typedef typename parent::type type;

	matrix() : parent() {} 

	/*
	math-index: i is the row, j is the column
	row-major: i is nested inner-most
	*/
	static int index(int i, int j) { return i + dim * j; }
	
};

template<int dim_, typename type_>
struct symmat;

template<int dim, typename type>
struct symmat_det {
	type operator()(const symmat<dim,type> &a);
};

template<int dim, typename type>
struct symmat_invert {
	symmat<dim,type> operator()(const symmat<dim,type> &a);
};

/*
symmat(i,j) == symmat(j,k)
*/
template<int dim_, typename type_>
struct symmat : public generic_dense_matrix<dim_, type_, symmat<dim_, type_>, dim_ * (dim_ + 1) / 2> {
	typedef generic_dense_matrix<dim_, type_, symmat<dim_, type_>, dim_ * (dim_ + 1) / 2> parent;
	
	enum { dim = parent::dim };
	typedef typename parent::type type;

	symmat() : parent() {}
	symmat(const symmat &a) : parent(a) {}

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

	//determinant
	static type det(const symmat &a) {
		return symmat_det<dim,type>()(a);
	}
	
	//inverse
	static symmat invert(const symmat &a) {
		return symmat_invert<dim,type>()(a);
	}
};

template<typename type>
struct symmat_det<2,type> {
	type operator()(const symmat<2,type> &a) const {
		return a(0,0) * a(1,1) - a(0,1) * a(0,1);
	}
};

template<typename type>
struct symmat_invert<2,type> {
	symmat<2,type> operator()(const symmat<2,type> &a) const {
		type det = symmat<2,type>::det(a);
		symmat<2,type> inv;
		inv(0,0) = a(1,1) / det;
		inv(1,1) = a(0,0) / det;
		inv(0,1) = -a(0,1) / det;
		return inv;
	}
};

//rank is templated, but dim is not as it varies per-rank
//so this is dynamically-sized
template<int rank_, typename type_>
struct Grid {
	enum { rank = rank_ };
	typedef type_ type;

	typedef ::vec<rank,int> veci;

	type *v;
	veci size;
	
	//cached for quick access by dot with index vector
	//step[0] = 1, step[1] = size[0], step[j] = product(i=1,j-1) size[i]
	veci step;

	Grid(const veci &size_) : size(size_) {
		v = new type[size.volume()];
		step(0) = 1;
		for (int i = 1; i < rank; ++i) {
			step(i) = step(i-1) * size(i-1);
		}
	}

	struct iterator {
		Grid *parent;
		veci index;
		
		iterator() : parent(NULL) {}
		iterator(Grid *parent_) : parent(parent_) {}
		iterator(const iterator &iter) : parent(iter.parent), index(iter.index) {}
		
		bool operator==(const iterator &b) const { return index == b.index; }
		bool operator!=(const iterator &b) const { return index != b.index; }
		
		iterator &operator++() {
			for (int i = 0; i < rank-1; ++i) {	//allow the last index to overflow for sake of comparing it to end
				++index(i);
				if (index(i) < parent->size(i)) break;
				index(i) = 0;
			}
			return *this;
		}

		Grid::type &operator*() const { return (*parent)(index); }
	};

	iterator begin() {
		iterator i(this);
		--i.index(0);
		return i;
	}
	iterator end() {
		iterator i(this);
		i.index(rank-1) = size(rank-1);
		return i;
	}

	type &operator()(const veci &index) { return v[veci::dot(veci::clamp(index, veci(0), size-1), step)]; }
	const type &operator()(const veci &index) const { return v[veci::dot(veci::clamp(index, veci(0), size-1), step)]; }
};

template<int dim_, typename real_>
struct Cell {
	typedef real_ real;
	enum { dim = dim_ };

	typedef ::vec<dim,real> vec;
	typedef ::symmat<dim,real> symmat;

	//	state variables

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

	// aux values computed and stored for partials

	//beta_i
	//beta_t = beta^k beta_k
	vec beta_l;
};

template<int dim_, typename real_>
struct Simulation {
	enum { dim = dim_ };
	typedef real_ real;

	//statically-sized mathematical types 
	typedef ::vec<dim,real> vec;		//3D hypersurface vector
	typedef ::vec<dim,int> veci;		//integer vector, used for indexes
	typedef ::vec<dim+1,real> tvec;		//4D hypersurface + timelike component vector
	typedef ::matrix<dim,real> matrix;	//3D matrix / rank-2 tensor of some sort
	typedef ::symmat<dim,real> symmat;	//3D symmetric matrix / rank-2 tensor of some sort

	typedef ::Cell<dim,real> Cell;
	typedef ::Grid<dim,Cell> CellGrid;
	CellGrid readCells, writeCells;

	//kronecher delta for integer vector / indexes, such that dxi(i)(j) == i == j
	::vec<dim,veci> dxi;

	//resolution of our grids, stored here as well as in each grid for convenience
	veci size;

	//x(0) == min
	vec min;
	
	//x(size-1) == max
	vec max;
	
	//range = max - min
	vec range;
	
	//dx = range / size
	vec dx;

	Simulation(const vec &min_, const vec &max_, const veci &size_)
	:	size(size_),
		min(min_),
		max(max_),
		range(max_ - min_),
		dx((max_ - min_) / vec(size_)),
		readCells(size_),
		writeCells(size_)
	{
		//basis vectors as integers for indexes
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j < dim; ++j) {
				dxi(i)(j) == i == j;
			}
		}
	}

	/*
	partial derivative operator
	for now let's use 3-point everywhere: d/dx^i f(x) ~= (f(x + dx^i) - f(x - dx^i)) / (2 * |dx^i|)
	*/
	template<typename T>
	T partial(T Cell::*field, const veci &index, int partialIndex) {
		return (readCells(index + dxi(dim)).*field - readCells(index - dxi(dim)).*field) / (2. * dx(partialIndex));
	}

	void update(real dt) {
		typedef typename CellGrid::iterator CellGridIter;
		CellGridIter iter;
		
		//first compute and store aux values that will be subsequently used for partial differentiation
		//during this process read and write to the same cell
		for (iter = writeCells.begin(); iter != writeCells.end(); ++iter) {
			Cell &readCell = readCells(iter.index);
		
			//beta_u(k) := beta^k
			const vec &beta_u = readCell.beta_u;
			
			//gamma_ll(i,j) := gamma_ij
			const symmat &gamma_ll = readCell.gamma_ll;

			//beta_l(i) := g_ij beta^j
			//exclude sum of beta^t = 0
			//not storing beta_t here, since beta_t = beta^k beta_k
			vec &beta_l = readCell.beta_l;
			for (int i = 0; i < dim; ++i) {
				beta_l(i) = 0;
				for (int j = 0; j < dim; ++j) {
					beta_l(i) += gamma_ll(i,j) * beta_u(j);
				}
			}
		}


		for (iter = writeCells.begin(); iter != writeCells.end(); ++iter) {
			const Cell &readCell = readCells(iter.index);
			Cell &writeCell = writeCells(iter.index);
			
			const real &alpha = readCell.alpha;
			
			//K_ll(i,j) := K_ij
			const symmat &K_ll = readCell.K_ll;

			//gamma_ll(i,j) := gamma_ij
			const symmat &gamma_ll = readCell.gamma_ll;

			//gamma_uu(i,j) := gamma^ij = inverse of gamma_ij
			// I could write the tensor formula for matrix inverses out (see "Gravitation", exercise 5.5e)
			symmat gamma_uu = symmat::invert(gamma_ll);

			//partial_gamma_lll[k](i,j) := partial_k gamma_ij
			symmat partial_gamma_lll[dim];
			for (int k = 0; k < dim; ++k) {
				partial_gamma_lll[k] = partial(&Cell::gamma_ll, iter.index, k);
			}

			//3D hypersurface connection coefficients
			//only need spatial coefficients (since gamma^at = 0, so conn^t_ab = 0.  see "Numerical Relativity", p.48)
			//conn_lll[i](j,k) := conn_ijk = 1/2 (partial_k gamma_ij + partial_j gamma_ik - partial_i gamma_jk)
			symmat conn_lll[dim];
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					for (int k = 0; k <= j; ++k) {
						conn_lll[i](j,k) = .5 * (partial_gamma_lll[k](i,j) + partial_gamma_lll[j](i,k) - partial_gamma_lll[i](j,k));
					}
				}
			}

			//conn_ull[i](j,k) := conn^i_jk = gamma^il conn_ljk
			symmat conn_ull[dim];
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					for (int k = 0; k <= j; ++k) {
						conn_ull[i](j,k) = 0;
						for (int l = 0; l < dim; ++l) {
							conn_ull[i](j,k) += gamma_uu(i,l) * conn_lll[l](j,k);
						}
					}
				}
			}
			
			//beta_l(k) := beta_k
			const vec &beta_l = readCell.beta_l;

			//partial_beta_ll[j](i) := partial_j beta_i
			vec partial_beta_ll[dim];
			for (int j = 0; j < dim; ++j) {
				partial_beta_ll[j] = partial(&Cell::beta_l, iter.index, j);
			}
			

			//diff_beta_ll(j,i) := diff_j beta_i = partial_j beta_i - conn^k_ij beta_k
			matrix diff_beta_ll;
			for (int j = 0; j < dim; ++j) {
				for (int i = 0; i < dim; ++i) {
					diff_beta_ll(j,i) = partial_beta_ll[j](i);
					for (int k = 0; k < dim; ++k) {
						diff_beta_ll(j,i) -= conn_ull[k](i,j) * beta_l(k);
					}
				}
			}
			
			//D_i beta_j = diff_i beta_j 
			//-- but only for lower indexes.
			//doesn't work if either is upper (you get extra terms for the timelike component we're dropping)
			//(see "Numerical Relativity", exercise 2.30)
			const matrix &D_beta_ll = diff_beta_ll;
			
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
	typedef double real;
	enum { dim = 2 };
	typedef ::vec<dim,real> vec;
	typedef ::vec<dim,int> veci;
	typedef ::Simulation<dim,real> Simulation;

	real dt = .01;
	real dist = 2e+3;	//2km, schwarzschild radius of our sun
	Simulation sim(vec(-dist, -dist), vec(dist, dist), veci(10, 10));
	sim.update(dt);
}

