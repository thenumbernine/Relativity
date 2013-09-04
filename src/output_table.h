#pragma once

/*
trying to organize my output stuff instead of typing headers into two places and hoping they match up
*/
#include <ostream>
#include "grid.h"
#include "cell.h"
#include "admformalism.h"
#include "tensor.h"

struct CoordNames {
	static const char *names[];
};
const char *CoordNames::names[] = {"x", "y", "z"};

template<typename FieldType>
struct OutputField {
	static void outputName(std::ostream &o, const char *name) {
		o << name << "\t";
	}

	static void outputField(std::ostream &o,const FieldType &field) {
		o << field << "\t";
	}
};

template<typename real, int dim> 
struct OutputFieldVector {
	static void outputName(std::ostream &o, const char *name, const char *indexSymbol) {
		for (int i = 0; i < dim; ++i) {
			o << name << indexSymbol << CoordNames::names[dim] << "\t";
		}
	}
	static void outputField(std::ostream &o, const typename generic_rank1<dim>::template body<real, real> &body) {
		for (int i = 0; i < dim; ++i) {
			o << body.v[i] << "\t";
		}
	}
};

template<typename real, int dim>
struct OutputField<tensor<real, lower<dim>>> : public OutputFieldVector<real, dim> {
	static void outputName(std::ostream &o, const char *name) {
		//TODO put the "_" vs "^" within the upper and lower classes
		//then write that mpl method of extracting the upper/lower for a particular index
		//	as well as the nested body type for a particular index
		//	as well as the dimension for a particular index
		//	etc...
		OutputFieldVector<real, dim>::outputName(o, name, "_");
	}
	static void outputField(std::ostream &o, const tensor<real, lower<dim>> &t) {
		OutputFieldVector<real, dim>::outputField(o, t.body);
	}
};

template<typename real, int dim>
struct OutputField<tensor<real, upper<dim>>> : public OutputFieldVector<real, dim> {
	static void outputName(std::ostream &o, const char *name) {
		OutputFieldVector<real, dim>::outputName(o, name, "^");
	}
	static void outputField(std::ostream &o, const tensor<real, upper<dim>> &t) {
		OutputFieldVector<real, dim>::outputField(o, t.body);
	}
};


//TODO write iterators could simplify this a lot
template<typename real, int dim>
struct OutputField<tensor<real, symmetric<lower<dim>, lower<dim>>>> {
	static void outputName(std::ostream &o, const char *name) {
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j <= i; ++j) {
				o << name << "_" << CoordNames::names[i] << CoordNames::names[j] << "\t";
			}
		}
	}
	static void outputField(std::ostream &o, const tensor<real, symmetric<lower<dim>, lower<dim>>> &t) {
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j <= i; ++j) {
				o << t(i,j) << "\t";
			}
		}
	}
};

template<typename real, int dim>
struct OutputField<tensor<real, symmetric<upper<dim>, upper<dim>>>> {
	static void outputName(std::ostream &o, const char *name) {
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j <= i; ++j) {
				o << name << "^" << CoordNames::names[i] << CoordNames::names[j] << "\t";
			}
		}
	}
	static void outputField(std::ostream &o, const tensor<real, symmetric<upper<dim>, upper<dim>>> &t) {
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j <= i; ++j) {
				o << t(i,j) << "\t";
			}
		}
	}
};


template<typename CellType>
struct ICellField {
	virtual const char *getFullName() const = 0;
	virtual void outputName(std::ostream &o) const = 0;
	virtual void outputValues(std::ostream &o, const CellType &cell) const = 0;
};

template<typename FieldType, typename CellType>
struct CellField : public ICellField<CellType> {
	const char *name;
	const char *fullname;
	FieldType CellType::*field;
	CellField(const char *name_, const char *fullname_, FieldType CellType::*field_)
	: name(name_), fullname(fullname_), field(field_) {}

	virtual const char *getFullName() const {
		return fullname; 
	}

	virtual void outputName(std::ostream &o) const {
		OutputField<FieldType>::outputName(o, name);
	}

	virtual void outputValues(std::ostream &o, const CellType &cell) const {
		OutputField<FieldType>::outputField(o, cell.*field);
	}

};
template<typename FieldType, typename CellType>
CellField<FieldType, CellType> *makeField(const char *name, const char *fullname, FieldType CellType::*field) {
	return new CellField<FieldType, CellType>(name, fullname, field);
}

template<typename real, int dim, typename CellType>
struct OutputCellFields;

//I've tried to avoid template-macro combos up until now ...
//stupid C++ ... you couldn't think of a more convoluted way to communicate compile-time information
//Boost isn't the *answer* to this.  It is the *problem* with this.
#define BEGIN_CELL_FIELDS(CELL)	\
template<typename real, int dim>	\
struct OutputCellFields<real, dim, CELL<real, dim>> {	\
	typedef CELL<real, dim> CellType;	\
	static ICellField<CellType> **fields() {	\
		static ICellField<CellType> *fields[] = {

#define CELL_FIELD(FIELD, SUFFIX...)	makeField(#FIELD, #FIELD #SUFFIX, &CellType::FIELD##SUFFIX),

#define END_CELL_FIELDS()	\
			NULL	\
		};	\
		return fields;	\
	}	\
};

BEGIN_CELL_FIELDS(GeomCell)
	CELL_FIELD(alpha)
	CELL_FIELD(beta, _u)
	CELL_FIELD(phi)
	CELL_FIELD(gamma, _ll)
	CELL_FIELD(K, _ll)
	CELL_FIELD(K)
END_CELL_FIELDS()

BEGIN_CELL_FIELDS(MatterCell)
	CELL_FIELD(rho)
	CELL_FIELD(S, _u)
	CELL_FIELD(S, _ll)
END_CELL_FIELDS()

BEGIN_CELL_FIELDS(AuxCell)
	CELL_FIELD(beta, _l)
	CELL_FIELD(psi)
	CELL_FIELD(gammaBar, _ll)
	CELL_FIELD(gammaBar, _uu)
	//CELL_FIELD(connBar, _ull)	//haven't written a OutputField for usl yet
	CELL_FIELD(RBar, _ll)
	CELL_FIELD(RBar)
	CELL_FIELD(gamma)
	CELL_FIELD(gamma, _uu)
	//CELL_FIELD(conn, _ull)
	CELL_FIELD(R, _ll)
	CELL_FIELD(R)
	CELL_FIELD(H)
	CELL_FIELD(M, _u)
END_CELL_FIELDS()

template<typename real, int dim>
struct OutputTable {
	typedef ::ADMFormalism<real, dim> ADMFormalism;
	typedef ::GeomCell<real, dim> GeomCell;
	typedef ::AuxCell<real, dim> AuxCell;
	typedef ::MatterCell<real, dim> MatterCell;

	template<typename CellType> 
	static void outputHeadersForCellType(std::ostream &o, const std::vector<bool> &columns, int &columnIndex) {
		typedef ::ICellField<CellType> ICellField;
		ICellField **fields = OutputCellFields<real, dim, CellType>::fields();
		for (int i = 0; ; ++i) {
			ICellField *field = fields[i];
			if (!field) break;
			if (columns[columnIndex++]) field->outputName(o);
		}
	}

	template<typename CellType> 
	static void outputValuesForCellType(std::ostream &o, const CellType &cell, const std::vector<bool> &columns, int &columnIndex) {
		typedef ::ICellField<CellType> ICellField;
		ICellField **fields = OutputCellFields<real, dim, CellType>::fields();
		for (int i = 0; ; ++i) {
			ICellField *field = fields[i];
			if (!field) break;
			if (columns[columnIndex++]) field->outputValues(o, cell);
		}
	}
	
	template<typename CellType> 
	static bool findColumnName(const std::string &fullname, int &columnIndex) {
		typedef ::ICellField<CellType> ICellField;
		ICellField **fields = OutputCellFields<real, dim, CellType>::fields();
		for (int i = 0; ; ++i) {
			ICellField *field = fields[i];
			if (!field) break;
			if (fullname == field->getFullName()) return true;
			++columnIndex;
		}
		return false;
	}
	
	template<typename CellType> 
	static void countNumColumns(int &columnIndex) {
		typedef ::ICellField<CellType> ICellField;
		ICellField **fields = OutputCellFields<real, dim, CellType>::fields();
		for (int i = 0; ; ++i) {
			if (!fields[i]) break;
			++columnIndex;
		}
	}
	
	static void header(std::ostream &o, std::vector<bool> &columns) {
		o << "#";
		o << "t\t";
		for (int i = 0; i < dim; ++i) {
			o << CoordNames::names[i] << "\t";
		}
		int columnIndex = 0;
		outputHeadersForCellType<GeomCell>(o, columns, columnIndex);
		outputHeadersForCellType<MatterCell>(o, columns, columnIndex);
		outputHeadersForCellType<AuxCell>(o, columns, columnIndex);
		o << std::endl;
	}

	static void state(std::ostream &o, const ADMFormalism &sim, std::vector<bool> &columns) {
		for (typename Grid<AuxCell, dim>::const_iterator iter = sim.auxGrid.begin(); iter != sim.auxGrid.end(); ++iter) {
			const AuxCell &auxCell = *iter;
			const MatterCell &matterCell = sim.matterGrid(iter.index);
			const GeomCell &geomCell = (*sim.geomGridReadCurrent)(iter.index);
			
			o << sim.time << "\t";

			::vector<real, dim> x = sim.coordForIndex(iter.index);
			for (int i = 0; i < dim; ++i) {
				o << x(i) << "\t";
			}

			int columnIndex = 0;
			outputValuesForCellType<GeomCell>(o, geomCell, columns, columnIndex);
			outputValuesForCellType<MatterCell>(o, matterCell, columns, columnIndex);
			outputValuesForCellType<AuxCell>(o, auxCell, columns, columnIndex);
			o << std::endl;
		}
		o << std::endl;
	}

	//only handles fullnames.  i.e. handles 'gamma_ll' rather than just 'gamma' for gamma_ij
	static int getColumnIndex(const std::string &columnName) {
		int columnIndex = 0;
		if (findColumnName<GeomCell>(columnName, columnIndex)) return columnIndex;
		if (findColumnName<MatterCell>(columnName, columnIndex)) return columnIndex;
		if (findColumnName<AuxCell>(columnName, columnIndex)) return columnIndex;
		throw Exception() << "failed to find column named " << columnName;
	}

	static int getNumColumns() {
		int columnIndex = 0;
		countNumColumns<GeomCell>(columnIndex);
		countNumColumns<MatterCell>(columnIndex);
		countNumColumns<AuxCell>(columnIndex);
		return columnIndex;
	}
};

