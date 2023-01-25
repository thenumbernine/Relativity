#pragma once

/*
trying to organize my output stuff instead of typing headers into two places and hoping they match up
*/
#include "cell.h"
#include "admformalism.h"
#include "parallel.h"
#include "Tensor/Grid.h"
#include "Tensor/Tensor.h"
#include "Common/Exception.h"
#include <ostream>

struct CoordNames {
	static char const * names[];
};
char const * CoordNames::names[] = {"x", "y", "z"};

template<typename FieldType>
struct OutputField {
	static void outputName(std::ostream &o, char const * name) {
		o << name << "\t";
	}

	static void outputField(std::ostream &o, FieldType const & field) {
		o << field << "\t";
	}
};

template<typename Real, int dim> 
struct OutputFieldVector {
	static void outputName(std::ostream &o, char const * name, char const * indexSymbol) {
		for (int i = 0; i < dim; ++i) {
			o << name << indexSymbol << CoordNames::names[i] << "\t";
		}
	}
	static void outputField(std::ostream &o, Tensor::vec<Real, dim> const & t) {
		for (int i = 0; i < dim; ++i) {
			o << t[i] << "\t";
		}
	}
};

template<typename Real, int dim>
struct OutputField<Tensor::vec<Real, dim>> : public OutputFieldVector<Real, dim> {
	static void outputName(std::ostream &o, char const * name) {
		OutputFieldVector<Real, dim>::outputName(o, name, "_");	//"^");//TODO can't determine anymore, no longer using Upper and Lower
	}
	static void outputField(std::ostream &o, Tensor::vec<Real, dim> const & t) {
		OutputFieldVector<Real, dim>::outputField(o, t);
	}
};

//TODO write iterators could simplify this a lot
template<typename Real, int dim>
struct OutputField<Tensor::sym<Real, dim>> {
	static void outputName(std::ostream &o, char const * name) {
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j <= i; ++j) {
				o << name << "_" << CoordNames::names[i] << CoordNames::names[j] << "\t";	// TODO Upper vs Lower
			}
		}
	}
	static void outputField(std::ostream &o, Tensor::sym<Real, dim> const & t) {
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j <= i; ++j) {
				o << t(i,j) << "\t";
			}
		}
	}
};

template<typename CellType>
struct ICellField {
	virtual char const * getFullName() const = 0;
	virtual void outputName(std::ostream &o) const = 0;
	virtual void outputValues(std::ostream &o, CellType const & cell) const = 0;
};

template<typename FieldType, typename CellType>
struct CellField : public ICellField<CellType> {
	char const * name;
	char const * fullname;
	FieldType CellType::*field;
	CellField(char const * name_, char const * fullname_, FieldType CellType::*field_)
	: name(name_), fullname(fullname_), field(field_) {}

	virtual char const * getFullName() const {
		return fullname; 
	}

	virtual void outputName(std::ostream &o) const {
		OutputField<FieldType>::outputName(o, name);
	}

	virtual void outputValues(std::ostream &o, CellType const & cell) const {
		OutputField<FieldType>::outputField(o, cell.*field);
	}

};
template<typename FieldType, typename CellType>
CellField<FieldType, CellType> *makeField(char const * name, char const * fullname, FieldType CellType::*field) {
	return new CellField<FieldType, CellType>(name, fullname, field);
}

template<typename Real, int dim, typename CellType>
struct OutputCellFields;

//I've tried to avoid template-macro combos up until now ...
//stupid C++ ... you couldn't think of a more convoluted way to communicate compile-time information
//Boost isn't the *answer* to this.  It is the *problem* with this.
#define BEGIN_CELL_FIELDS(CELL)	\
template<typename Real, int dim>	\
struct OutputCellFields<Real, dim, CELL<Real, dim>> {	\
	using CellType = CELL<Real, dim>;	\
	static ICellField<CellType> **fields() {	\
		static ICellField<CellType> *fields[] = {

#define CELL_FIELD(FIELD, SUFFIX...)	makeField(#FIELD, #FIELD #SUFFIX, &CellType::FIELD##SUFFIX),

#define END_CELL_FIELDS()	\
			nullptr	\
		};	\
		return fields;	\
	}	\
};

BEGIN_CELL_FIELDS(GeomCell)
	CELL_FIELD(alpha)
	CELL_FIELD(beta, _u)
	CELL_FIELD(phi)
	CELL_FIELD(gammaBar, _ll)
	CELL_FIELD(ATilde, _ll)
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
	CELL_FIELD(gammaBar, _uu)
	//CELL_FIELD(connBar, _ull)	//haven't written a OutputField for usl yet
	CELL_FIELD(RBar, _ll)
	CELL_FIELD(RBar)
	CELL_FIELD(gamma)
	CELL_FIELD(gamma, _ll)
	CELL_FIELD(gamma, _uu)
	//CELL_FIELD(conn, _ull)
	CELL_FIELD(RBar, _ll)
	CELL_FIELD(RPhi, _ll)
	CELL_FIELD(R)
	CELL_FIELD(H)
	CELL_FIELD(M, _u)
END_CELL_FIELDS()

template<typename Real, int dim>
struct OutputTable {
	using ADMFormalism = ::ADMFormalism<Real, dim>;
	using GeomCell = ::GeomCell<Real, dim>;
	using AuxCell = ::AuxCell<Real, dim>;
	using MatterCell = ::MatterCell<Real, dim>;

	template<typename CellType> 
	static void outputHeadersForCellType(std::ostream &o, std::vector<bool> const & columns, int & columnIndex) {
		using ICellField = ::ICellField<CellType>;
		ICellField **fields = OutputCellFields<Real, dim, CellType>::fields();
		for (int i = 0; ; ++i) {
			ICellField *field = fields[i];
			if (!field) break;
			if (columns[columnIndex++]) field->outputName(o);
		}
	}

	template<typename CellType> 
	static void outputValuesForCellType(std::ostream &o, CellType const & cell, std::vector<bool> const & columns, int &columnIndex) {
		using ICellField = ::ICellField<CellType>;
		ICellField **fields = OutputCellFields<Real, dim, CellType>::fields();
		for (int i = 0; ; ++i) {
			ICellField *field = fields[i];
			if (!field) break;
			if (columns[columnIndex++]) field->outputValues(o, cell);
		}
	}
	
	template<typename CellType> 
	static bool findColumnName(std::string const & fullname, int &columnIndex) {
		using ICellField = ::ICellField<CellType>;
		ICellField **fields = OutputCellFields<Real, dim, CellType>::fields();
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
		using ICellField = ::ICellField<CellType>;
		ICellField **fields = OutputCellFields<Real, dim, CellType>::fields();
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

	static void state(std::ostream &o, ADMFormalism const & sim, std::vector<bool> &columns) {
		Tensor::RangeObj<dim> range = sim.auxGrid.range();
		parallel.foreach(range.begin(), range.end(), [&](Tensor::intN<dim> index) {
			AuxCell const & auxCell = sim.auxGrid(index);
			MatterCell const & matterCell = sim.matterGrid(index);
			GeomCell const & geomCell = (*sim.geomGridReadCurrent)(index);
			
			o << sim.time << "\t";

			auto x = sim.coordForIndex(index);
			for (auto const & xi : x) {
				o << xi << "\t";
			}

			int columnIndex = 0;
			outputValuesForCellType<GeomCell>(o, geomCell, columns, columnIndex);
			outputValuesForCellType<MatterCell>(o, matterCell, columns, columnIndex);
			outputValuesForCellType<AuxCell>(o, auxCell, columns, columnIndex);
			o << std::endl;
		});
		o << std::endl;
	}

	//only handles fullnames.  i.e. handles 'gamma_ll' rather than just 'gamma' for gamma_ij
	static int getColumnIndex(std::string const & columnName) {
		int columnIndex = 0;
		if (findColumnName<GeomCell>(columnName, columnIndex)) return columnIndex;
		if (findColumnName<MatterCell>(columnName, columnIndex)) return columnIndex;
		if (findColumnName<AuxCell>(columnName, columnIndex)) return columnIndex;
		throw Common::Exception() << "failed to find column named " << columnName;
	}

	static int getNumColumns() {
		int columnIndex = 0;
		countNumColumns<GeomCell>(columnIndex);
		countNumColumns<MatterCell>(columnIndex);
		countNumColumns<AuxCell>(columnIndex);
		return columnIndex;
	}
};
