#pragma once

/*
the 'parent' curious pattern whatever for generic_vector and generic_matrix
	size_ is the size of the dense vector
	type_ is the type of each element
	child is the curious reoccurring chlid class that uses this
*/
template<int size_, typename type_, typename child>
struct generic_array {
	enum { size = size_ };
	typedef type_ type;

	type v[size];

	//default ctor: init all to zero
	generic_array() {
		for (int i = 0; i < size; ++i) {
			v[i] = type();
		}
	}

	//default copy ctor
	generic_array(const child &a) {
		for (int i = 0; i < size; ++i) {
			v[i] = a.v[i];
		}
	}

	//default fill ctor
	generic_array(const type &x) {
		for (int i = 0; i < size; ++i) {
			v[i] = x;
		}
	}

	//bounds
	static child clamp(const child &a, const child &min, const child &max) {
		child b;
		for (int i = 0; i < size; ++i) {
			b.v[i] = std::max(min.v[i], std::min(max.v[i], a.v[i]));
		}
		return b;
	}

	//boolean operations

	bool operator==(const child &b) const {
		const generic_array &a = *this;
		for (int i = 0; i < size; ++i) {
			if (a.v[i] != b.v[i]) return false;
		}
		return true;
	}

	bool operator!=(const child &b) const { return ! this->operator==(b); }

	//math operations

	child operator+(const child &b) const {
		const generic_array &a = *this;
		child c;
		for (int i = 0; i < size; ++i) {
			c.v[i] = a.v[i] + b.v[i];
		}
		return c;
	}
	
	child operator-(const child &b) const {
		const generic_array &a = *this;
		child c;
		for (int i = 0; i < size; ++i) {
			c.v[i] = a.v[i] - b.v[i];
		}
		return c;
	}

	child operator*(const type &b) const {
		const generic_array &a = *this;
		child c;
		for (int i = 0; i < size; ++i) {
			c.v[i] = a.v[i] * b;
		}
		return c;
	}


	child operator/(const type &b) const {
		const generic_array &a = *this;
		child c;
		for (int i = 0; i < size; ++i) {
			c.v[i] = a.v[i] / b;
		}
		return c;
	}
};

