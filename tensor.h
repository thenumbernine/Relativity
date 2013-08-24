#pragma once

/*
new experimental template-driven tensor class
made to unify all of vector, one-form, symmetric, antisymmetric, and just regular kind matrixes
and any subsequently nested numeric structures.

so far so good except 
*/

#include "generic_vector.h"
#include "symmat.h"			//used by symmetric
#include "antisymmat.h"		//used by antisymmetric

/*
vector is used by generic_rank1 
tensor is supposed to replace vectors
sooo
how about we rename 'vector' to 'tuple' or 'fixed_array' or something?
*/
#include "vector.h"

template<int dim_>
struct generic_rank1 {
	enum { dim = dim_ };
	enum { rank = 1 };

	//now i need a specialization of this for when 'type' is a primitive
	// in such cases, accessing nested members results in compiler errors
	template<typename inner_type, typename scalar_type>
	struct body : public generic_vector<dim, inner_type, scalar_type, body<inner_type, scalar_type>> {
		typedef generic_vector<dim, inner_type, scalar_type, body<inner_type, scalar_type>> parent;
		body() : parent() {}
		body(const body &b) : parent(b) {}
		body(const inner_type &t) : parent(t) {}

		inner_type &operator()(const ::vector<1,int> &deref) { return parent::v[deref(0)]; }
		const inner_type &operator()(const ::vector<1,int> &deref) const { return parent::v[deref(0)]; }
	};
};

/*
dim = dimension of this upper index
used like so: 
	typedef tensor<real, upper<3>> vector3;
*/
template<int dim>
struct upper : public generic_rank1<dim> {};

/*
dim = dimension of this lower index
used like so:
	typedef tensor<real, lower<3>> oneform3;
*/
template<int dim>
struct lower : public generic_rank1<dim> {};

/*
used like so: 
	typedef tensor<real, symmetric<lower<3>, lower<3>>> metric3;
	typedef tensor<real, symmetric<lower<3>, lower<3>>, symmetric<lower<3>, lower<3>>> riemann_metric3;
*/
template<typename index1, typename index2>
struct symmetric {
	static_assert(index1::dim == index2::dim, "symmetric can only accept two indexes of equal dimension");
	enum { dim = index1::dim };
	enum { rank = index1::rank + index2::rank };
	
	template<typename inner_type, typename scalar_type>
	struct body : public generic_symmat<index1::dim, inner_type, scalar_type, body<inner_type, scalar_type>> {
		typedef generic_symmat<index1::dim, inner_type, scalar_type, body<inner_type, scalar_type>> parent;
		body() : parent() {}
		body(const body &b) : parent(b) {}
		body(const inner_type &t) : parent(t) {}
	};
};

template<typename index1, typename index2>
struct antisymmetric {
	static_assert(index1::dim == index2::dim, "symmetric can only accept two indexes of equal dimension");
	enum { dim = index1::dim };
	enum { rank = index1::rank + index2::rank };

	template<typename inner_type, typename scalar_type>
	struct body : public generic_antisymmat<index1::dim, inner_type, scalar_type, body<inner_type, scalar_type>> {
		typedef generic_antisymmat<index1::dim, inner_type, scalar_type, body<inner_type, scalar_type>> parent;
		body() : parent() {}
		body(const body &b) : parent(b) {}
		body(const inner_type &t) : parent(t) {}
	};
};

/*
support class for metaprogram specializations of tensor

rank 		= total rank of the structure
				= the depth of the body_type
inner_type 	= the next-inner-most tensor_stats type
body_type 	= the body_type (the generic_vector subclass) associated with this nesting level
deref_type 	= the dereference type that would be needed to dereference this body_type 
				= int vector with dimension equal to rank
*/
template<typename scalar_type, typename... args>
struct tensor_stats;

template<typename scalar_type_, typename index_, typename... args>
struct tensor_stats<scalar_type_, index_, args...> {
	typedef scalar_type_ scalar_type;
	typedef index_ index;
	enum { rank = index::rank + tensor_stats<scalar_type, args...>::rank };

	typedef tensor_stats<scalar_type, args...> inner_type;

	typedef typename index::template body<typename inner_type::body_type, scalar_type> body_type;
	
	typedef vector<rank,int> deref_type;	//rank-n tensors (n>1) use vecs for dereferencing

		//inductive cases
	
	//get an element in a tensor
	template<int total_rank, int offset>
	static scalar_type &get(body_type &body, const ::vector<total_rank,int> &deref) {
		::vector<index::rank,int> subderef;
		for (int i = 0; i < index::rank; ++i) {
			subderef(i) = deref(i + offset);
		}
		return inner_type::template get<total_rank, offset + index::rank>(body(subderef), deref);
	}

	template<int total_rank, int offset>
	static const scalar_type &get_const(const body_type &body, const ::vector<total_rank,int> &deref) {
		::vector<index::rank,int> subderef;
		for (int i = 0; i < index::rank; ++i) {
			subderef(i) = deref(i + offset);
		}
		return inner_type::template get_const<total_rank, offset + index::rank>(body(subderef), deref);
	}

	//get the size of the tensor
	template<int total_rank, int offset>
	static void assign_size(::vector<total_rank,int> &s) {
		for (int i = offset; i < offset + index::rank; ++i) {
			s.v[i] = index::dim;
		}
		inner_type::template assign_size<total_rank, offset + index::rank>(s);
	}
};

template<typename scalar_type_, typename index_>
struct tensor_stats<scalar_type_, index_> {
	typedef scalar_type_ scalar_type;
	typedef index_ index;
	enum { rank = index::rank };

	typedef tensor_stats<scalar_type> inner_type;
	
	typedef typename index::template body<scalar_type, scalar_type> body_type;

	typedef ::vector<rank,int> deref_type;	//rank-1 tensors use ints for dereferencing

		//second-to-base (could be base) case
	
	//get an element in a tensor
	template<int total_rank, int offset>
	static scalar_type &get(body_type &body, const ::vector<total_rank,int> &deref) {
		::vector<index::rank,int> subderef;
		for (int i = 0; i < index::rank; ++i) {
			subderef(i) = deref(i + offset);
		}
		return body(subderef);
	}

	template<int total_rank, int offset>
	static const scalar_type &get_const(const body_type &body, const ::vector<total_rank,int> &deref) {
		::vector<index::rank,int> subderef;
		for (int i = 0; i < index::rank; ++i) {
			subderef(i) = deref(i + offset);
		}
		return body(subderef);
	}

	//get the size of the tensor
	template<int total_rank, int offset>
	static void assign_size(::vector<total_rank,int> &s) {
		for (int i = offset; i < offset + index::rank; ++i) {
			s.v[i] = index::dim;
		}
	}
};

//appease the vararg specialization recursive reference
//(it doesn't seem to recognize the single-entry base case)
template<typename scalar_type_>
struct tensor_stats<scalar_type_> {
	typedef scalar_type_ scalar_type;
	enum { rank = 0 };

	typedef scalar_type body_type;

		//base case: do nothing
	
		//get an element in a tensor
	
	template<int total_rank, int offset>
	static scalar_type &get(body_type &body, const ::vector<total_rank,int> &deref) {}
	
	template<int total_rank, int offset>
	static const scalar_type &get_const(const body_type &body, const ::vector<total_rank,int> &deref) {}
	
	template<int total_rank, int offset>
	static void assign_size(::vector<total_rank,int> &s) {}
};

/*
3-element vector:					tensor<double,upper<3> >
4-element one-form:					tensor<double,lower<3> >
matrix (contravariant matrix):		tensor<double, upper<3>, upper<3> >
metric tensor (symmetric covariant matrix):	tensor<double, symmetric<lower,3> >
*/
template<typename type_, typename... args>
struct tensor {
	typedef type_ type;

	//tensor_stats metaprogram calculates rank
	//it pulls individual entries from the index args
	//I could have the index args themselves do the calculation
	// but that would mean making base-case specializations for each index class 
	typedef ::tensor_stats<type_, args...> tensor_stats;
	typedef typename tensor_stats::body_type body_type;
	typedef typename tensor_stats::deref_type deref_type;
	
	enum { rank = tensor_stats::rank };

	tensor() {}
	tensor(const body_type &body_) : body(body_) {}
	tensor(const tensor &t) : body(t.body) {}
	tensor(const type &v) : body(v) {}
	
	/*
	tensor<real, upper<3>> v;
	v.body(0) will return of type real

	tensor<real, upper<3>, upper<4>> v;
	v.body(0) will return of type upper<4>::body<real, real>
	*/
	body_type body;

	//no way to specify a typed list of arguments
	// (solving that problem would give us arbitrary parameter constructors for vector classes)
	//so here's the special case instances for up to N=4
	typename tensor_stats::inner_type::body_type &operator()(int i) { return body(i); }
	const typename tensor_stats::inner_type::body_type &operator()(int i) const { return body(i); }
	//...and I haven't implemented these yet ...
	type &operator()(int i, int j) { return tensor_stats::template get<2,0>(body, ::vector<2,int>(i, j)); }
	const type &operator()(int i, int j) const { return tensor_stats::template get_const<2,0>(body, ::vector<2,int>(i, j)); }
	type &operator()(int i, int j, int k) { return tensor_stats::template get<3,0>(body, ::vector<3,int>(i, j, k)); }
	const type &operator()(int i, int j, int k) const { return tensor_stats::template get_const<3,0>(body, ::vector<3,int>(i, j, k)); }
	type &operator()(int i, int j, int k, int l) { return tensor_stats::template get<4,0>(body, ::vector<4,int>(i, j, k, l)); }
	const type &operator()(int i, int j, int k, int l) const { return tensor_stats::template get_const<4,0>(body, ::vector<4,int>(i, j, k, l)); }
	
	type &operator()(const deref_type &deref) { return tensor_stats::template get<rank,0>(body, deref); }
	const type &operator()(const deref_type &deref) const { return tensor_stats::template get_const<rank,0>(body, deref); }

	tensor operator-() const { return tensor(-body); }
	tensor operator+(const tensor &b) const { return tensor(body + b.body); }
	tensor operator-(const tensor &b) const { return tensor(body - b.body); }
	tensor operator*(const type &b) const { return tensor(body * b); }
	tensor operator/(const type &b) const { return tensor(body / b); }

	//tempted to make this a generic ctor,
	// but that would probably screw up the rest of our ctors
	template<typename type2, typename... args2>
	tensor &operator=(const ::tensor<type2, args2...> &t) {
		for (tensor::iterator i = begin(); i != end(); ++i) {
			*i = t(i.index);
		}
	}

	deref_type size() const {
		deref_type s;
		tensor_stats::template assign_size<rank, 0>(s);
		return s;
	};

	struct iterator {
		tensor *parent;
		deref_type index;
		
		iterator() : parent(NULL) {}
		iterator(tensor *parent_) : parent(parent_) {}
		iterator(const iterator &iter) : parent(iter.parent), index(iter.index) {}
		
		bool operator==(const iterator &b) const { return index == b.index; }
		bool operator!=(const iterator &b) const { return index != b.index; }
		
		//NOTICE this doesn't take symmetric indexes into account
		// we could shave off a few assignments with that
		iterator &operator++() {
			for (int i = 0; i < rank-1; ++i) {	//allow the last index to overflow for sake of comparing it to end
				++index(i);
				if (index(i) < parent->size()(i)) break;
				index(i) = 0;
			}
			return *this;
		}

		type &operator*() const { return (*parent)(index); }
	};

	iterator begin() {
		iterator i(this);
		--(i.index(0));
		return i;
	}
	iterator end() {
		iterator i(this);
		i.index(rank-1) = size()(rank-1);
		return i;
	}

	operator body_type&() { return body; }
	operator const body_type&() const { return body; }
};

