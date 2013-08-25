#pragma once

#include "vector.h"

//rank is templated, but dim is not as it varies per-rank
//so this is dynamically-sized
template<typename type_, int rank_>
struct Grid {
	typedef type_ type;
	enum { rank = rank_ };

	typedef vector<int,rank> deref_type;

	type *v;
	deref_type size;
	
	//cached for quick access by dot with index vector
	//step[0] = 1, step[1] = size[0], step[j] = product(i=1,j-1) size[i]
	deref_type step;

	Grid(const deref_type &size_) : size(size_) {
		v = new type[size.volume()];
		step(0) = 1;
		for (int i = 1; i < rank; ++i) {
			step(i) = step(i-1) * size(i-1);
		}
	}

	struct iterator {
		Grid *parent;
		deref_type index;
		
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

	type &operator()(const deref_type &index) { return v[deref_type::dot(index, step)]; }
	const type &operator()(const deref_type &index) const { return v[deref_type::dot(index, step)]; }
};

