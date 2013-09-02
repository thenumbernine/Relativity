#pragma once

/*
trying to organize my output stuff instead of typing headers into two places and hoping they match up
*/
#include <ostream>
#include "grid.h"
#include "cell.h"
#include "admformalism.h"

template<typename CellType>
struct ICellField {
	virtual void output(CellType &cell, std::ostream &o) = 0;
	virtual const char *getName() = 0;
};

template<typename FieldType, typename CellType>
struct CellField : public ICellField<CellType> {
	const char *name;
	FieldType CellType::*field;
	CellField(const char *name_, FieldType CellType::*field_)
	: name(name_), field(field_) {}
	virtual void output(CellType &cell, std::ostream &o) {
		return o << cell.*field;
	}
	virtual const char *getName() const { return name; }

};
template<typename FieldType, typename CellType>
CellField<FieldType, CellType> *makeField(const char *name, FieldType CellType::*field) {
	return new CellField<FieldType, CellType>(name, field);
}

template<typename real, int dim, typename CellType>
struct OutputCellFields;

//I've tried to avoid template-macro combos up until now ...
//stupid C++ ... you couldn't think of a more convoluted way to communicate compile-time information
//Boost isn't the *answer* to this.  It is the *problem* with this.
#define BEGIN_CELL_FIELDS(CELL)	\
template<typename real, int dim>	\
struct OutputCellFields<real, dim, CELL<real, dim>> {	\
	typedef GeomCell<real, dim> CellType;	\
	static ICellField<CellType> **fields() {	\
		static ICellField<CellType> *fields[] = {

#define CELL_FIELD(FIELD)	\
			makeField(#FIELD, &CellType::FIELD),

#define END_CELL_FIELDS()	\
		};	\
		return fields;	\
	}	\
};

BEGIN_CELL_FIELDS(GeomCell)
	CELL_FIELD(alpha)
	CELL_FIELD(beta_u)
	CELL_FIELD(ln_sqrt_gamma)
	CELL_FIELD(gamma_ll)
	CELL_FIELD(K_ll)
	CELL_FIELD(K)
END_CELL_FIELDS()

template<typename real, int dim>
struct OutputTable {
	typedef ::ADMFormalism<real, dim> ADMFormalism;
	typedef ::GeomCell<real, dim> GeomCell;
	typedef ::AuxCell<real, dim> AuxCell;
	typedef ::MatterCell<real, dim> MatterCell;

	typedef OutputCellFields<real, dim, GeomCell> OutputGeomCellFields;
	typedef OutputCellFields<real, dim, AuxCell> OutputAuxCellFields;
	typedef OutputCellFields<real, dim, MatterCell> OutputMatterCellFields;
	
	void header(std::ostream &o) {
		typedef ICellField<GeomCell> GeomField;
		GeomField **geomFields = OutputGeomCellFields::fields();
		for (int i = 0; i < numberof(geomFields); ++i) {
			GeomField *field = geomFields[i];
			o << field->getName() << "\t";
		}

		o << std::endl;
	}

	void state(const ADMFormalism &sim) {
	}
};


