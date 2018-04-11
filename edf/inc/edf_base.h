#ifndef EDF_BASEH
#define EDF_BASEH
// ===========================================================================
#include "utils.h"
#include "potential.h"
#include "tables_aa.h"
// ===========================================================================
// Base class for EDF
class edf{
protected:
public:
	// Potential and Action calculator
	Potential_JS *pot;
	Actions_AxisymmetricFudge_InterpTables *ActionCalculator;

	edf(Potential_JS *Pot, Actions_AxisymmetricFudge_InterpTables *ActionCalculator):pot(Pot),ActionCalculator(ActionCalculator){}
};

struct DF_st{
	edf *EDF; VecDoub acts;
	DF_st(edf *EDF,VecDoub acts)
		:EDF(EDF), acts(acts){}
};
// ===========================================================================
#endif
