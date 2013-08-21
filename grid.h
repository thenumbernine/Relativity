#pragma once

#include "vec.h"

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

