package methods

import (
	"errors"

	gcf "github.com/NumberXNumbers/types/gc/functions"
	"github.com/NumberXNumbers/types/gc/functions/arguments"
	gcv "github.com/NumberXNumbers/types/gc/values"
)

// evalV will evaluates a gcf function at inputs and expects a gcv Value, else returns error
func evalV(f *gcf.Function, inputs ...interface{}) (gcv.Value, error) {
	solution, err := f.Eval(inputs...)
	if err != nil {
		return nil, err
	}

	if solution.Type() != arguments.Value {
		return nil, errors.New("Incorrect type")
	}

	return solution.Value(), nil
}
