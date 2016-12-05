package methods

import (
	"testing"

	gcf "github.com/NumberXNumbers/types/gc/functions"
	gcfargs "github.com/NumberXNumbers/types/gc/functions/arguments"
	v "github.com/NumberXNumbers/types/gc/vectors"
)

func TestEvalV(t *testing.T) {
	x := gcfargs.NewVar(gcfargs.Vector)
	regVars := []gcfargs.Var{x}
	testFunction := gcf.MakeFuncPanic(regVars, x)

	_, err := evalV(testFunction, v.NewVector(v.RowSpace, 3))
	if err == nil {
		t.Fail()
	}
}
