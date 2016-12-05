package methods

import (
	"errors"
	"math"

	gcf "github.com/NumberXNumbers/types/gc/functions"
	gcv "github.com/NumberXNumbers/types/gc/values"
	gcvops "github.com/NumberXNumbers/types/gc/values/ops"
)

// Bisection1D is for solving the 1D root finding bisection method
func Bisection1D(intervalBegin float64, intervalEnd float64, TOL float64, maxIteration int, f *gcf.Function) (gcv.Value, error) {
	fOfA, errfA := evalV(f, intervalBegin)
	if errfA != nil {
		return nil, errfA
	}

	currentX := intervalBegin + (intervalEnd-intervalBegin)/float64(2)
	fOfCurrentX, errCurrentX := evalV(f, currentX)
	if errCurrentX != nil {
		return nil, errCurrentX
	}

	var root gcv.Value
	solutionFound := false

	for i := 0; i < maxIteration; i++ {
		if fOfCurrentX.Real() == 0 || (intervalEnd-intervalBegin)/float64(2) < TOL {
			root = gcv.MakeValue(currentX)
			solutionFound = true
			break
		}

		if gcvops.Mult(fOfA, fOfCurrentX).Real() > 0 {
			intervalBegin = currentX
			fOfA = fOfCurrentX
		} else {
			intervalEnd = currentX
		}

		currentX = intervalBegin + (intervalEnd-intervalBegin)/float64(2)
		fOfCurrentX, errCurrentX = evalV(f, currentX)
		if errCurrentX != nil {
			return nil, errCurrentX
		}
	}

	if solutionFound {
		return root, nil
	}

	return nil, errors.New("Unable to find root of given function")
}

// FixedPointIteration1D is for solving the 1D root finding fixed point iteration method
func FixedPointIteration1D(initialApprox float64, TOL float64, maxIteration int, f *gcf.Function) (gcv.Value, error) {
	previousApprox := gcv.MakeValue(initialApprox)
	currentApprox, errCurrentApprox := evalV(f, previousApprox)
	if errCurrentApprox != nil {
		return nil, errCurrentApprox
	}

	var root gcv.Value
	solutionFound := false

	for i := 0; i < maxIteration; i++ {
		if gcvops.Abs(gcvops.Sub(currentApprox, previousApprox)).Real() < TOL {
			root = gcv.MakeValue(currentApprox)
			solutionFound = true
			break
		}

		previousApprox = currentApprox
		currentApprox, errCurrentApprox = evalV(f, previousApprox)
		if errCurrentApprox != nil {
			return nil, errCurrentApprox
		}
	}

	if solutionFound {
		return root, nil
	}

	return nil, errors.New("Unable to find root of given function")
}

// Newton1D is for solving the 1D  root finding newton's method
func Newton1D(initialApprox float64, TOL float64, maxIteration int, f *gcf.Function, df *gcf.Function) (gcv.Value, error) {
	previousApprox := gcv.MakeValue(initialApprox)
	fPA, errfPA := evalV(f, previousApprox)
	if errfPA != nil {
		return nil, errfPA
	}

	dfPA, errdfPA := evalV(df, previousApprox)
	if errdfPA != nil {
		return nil, errdfPA
	}

	currentApprox := gcvops.Sub(previousApprox, gcvops.Div(fPA, dfPA))
	var root gcv.Value
	solutionFound := false

	for i := 0; i < maxIteration; i++ {
		if gcvops.Abs(gcvops.Sub(currentApprox, previousApprox)).Real() < TOL {
			root = currentApprox
			solutionFound = true
			break
		}

		previousApprox = currentApprox

		fPA, errfPA = evalV(f, previousApprox)
		if errfPA != nil {
			return nil, errfPA
		}

		dfPA, errdfPA = evalV(df, previousApprox)
		if errdfPA != nil {
			return nil, errdfPA
		}

		currentApprox = gcvops.Sub(previousApprox, gcvops.Div(fPA, dfPA))
	}

	if solutionFound {
		return root, nil
	}

	return nil, errors.New("Unable to find root of given function")
}

// ModifiedNewton1D is a modification for solving the 1D  root finding newton's method
func ModifiedNewton1D(initialApprox float64, TOL float64, maxIteration int, f *gcf.Function,
	df *gcf.Function, ddf *gcf.Function) (gcv.Value, error) {
	previousApprox := gcv.MakeValue(initialApprox)
	two := gcv.MakeValue(2)

	fPA, errfPA := evalV(f, previousApprox)
	if errfPA != nil {
		return nil, errfPA
	}

	dfPA, errdfPA := evalV(df, previousApprox)
	if errdfPA != nil {
		return nil, errdfPA
	}

	ddfPA, errddfPA := evalV(ddf, previousApprox)
	if errddfPA != nil {
		return nil, errddfPA
	}

	ratioA := gcvops.Mult(fPA, dfPA)
	ratioB := gcvops.Mult(fPA, ddfPA)
	currentApprox := gcvops.Sub(previousApprox, gcvops.Div(ratioA, gcvops.Sub(gcvops.Pow(dfPA, two), ratioB)))
	var root gcv.Value
	solutionFound := false

	for i := 0; i < maxIteration; i++ {
		if gcvops.Abs(gcvops.Sub(currentApprox, previousApprox)).Real() < TOL {
			root = currentApprox
			solutionFound = true
			break
		}

		previousApprox = currentApprox
		fPA, errfPA = evalV(f, previousApprox)
		if errfPA != nil {
			return nil, errfPA
		}

		dfPA, errdfPA = evalV(df, previousApprox)
		if errdfPA != nil {
			return nil, errdfPA
		}

		ddfPA, errddfPA = evalV(ddf, previousApprox)
		if errddfPA != nil {
			return nil, errddfPA
		}

		ratioA := gcvops.Mult(fPA, dfPA)
		ratioB := gcvops.Mult(fPA, ddfPA)
		currentApprox = gcvops.Sub(previousApprox, gcvops.Div(ratioA, gcvops.Sub(gcvops.Pow(dfPA, two), ratioB)))
	}

	if solutionFound {
		return root, nil
	}

	return nil, errors.New("Unable to find root of given function")
}

// Secant1D is for solving the 1D root finding secant method
func Secant1D(initialApprox1 float64, initialApprox2 float64, TOL float64, maxIteration int, f *gcf.Function) (gcv.Value, error) {
	previousApprox1 := gcv.MakeValue(initialApprox1)
	previousApprox2 := gcv.MakeValue(initialApprox2)

	fPA1, errfPA1 := evalV(f, previousApprox1)
	if errfPA1 != nil {
		return nil, errfPA1
	}

	fPA2, errfPA2 := evalV(f, previousApprox2)
	if errfPA2 != nil {
		return nil, errfPA2
	}

	ratioA := gcvops.Div(gcvops.Sub(previousApprox2, previousApprox1), gcvops.Sub(fPA2, fPA1))
	currentApprox := gcvops.Sub(previousApprox2, gcvops.Mult(fPA2, ratioA))
	var root gcv.Value
	solutionFound := false

	for i := 1; i < maxIteration; i++ {
		if gcvops.Abs(gcvops.Sub(currentApprox, previousApprox2)).Real() < TOL {
			root = currentApprox
			solutionFound = true
			break
		}

		previousApprox1 = previousApprox2
		previousApprox2 = currentApprox

		fPA1, errfPA1 = evalV(f, previousApprox1)
		if errfPA1 != nil {
			return nil, errfPA1
		}

		fPA2, errfPA2 = evalV(f, previousApprox2)
		if errfPA2 != nil {
			return nil, errfPA2
		}

		ratioA := gcvops.Div(gcvops.Sub(previousApprox2, previousApprox1), gcvops.Sub(fPA2, fPA1))
		currentApprox = gcvops.Sub(previousApprox2, gcvops.Mult(fPA2, ratioA))
	}

	if solutionFound {
		return root, nil
	}

	return nil, errors.New("Unable to find root of given function")
}

// FalsePosition1D is for solving the 1D root finding false position method
func FalsePosition1D(initialApprox1 float64, initialApprox2 float64, TOL float64, maxIteration int, f *gcf.Function) (gcv.Value, error) {
	previousApprox1 := gcv.MakeValue(initialApprox1)
	previousApprox2 := gcv.MakeValue(initialApprox2)

	fPA1, errfPA1 := evalV(f, previousApprox1)
	if errfPA1 != nil {
		return nil, errfPA1
	}

	fPA2, errfPA2 := evalV(f, previousApprox2)
	if errfPA2 != nil {
		return nil, errfPA2
	}

	ratioA := gcvops.Div(gcvops.Sub(previousApprox2, previousApprox1), gcvops.Sub(fPA2, fPA1))
	currentApprox := gcvops.Sub(previousApprox2, gcvops.Mult(fPA2, ratioA))

	fCA, errfCA := evalV(f, currentApprox)
	if errfCA != nil {
		return nil, errfCA
	}

	var root gcv.Value
	solutionFound := false

	for i := 1; i < maxIteration; i++ {
		if gcvops.Abs(gcvops.Sub(currentApprox, previousApprox2)).Real() < TOL {
			root = currentApprox
			solutionFound = true
			break
		}

		if gcvops.Mult(fCA, fPA2).Real() < 0 {
			previousApprox1 = previousApprox2
			fPA1 = fPA2
		}

		previousApprox2 = currentApprox
		fPA2 = fCA

		ratioA = gcvops.Div(gcvops.Sub(previousApprox2, previousApprox1), gcvops.Sub(fPA2, fPA1))
		currentApprox = gcvops.Sub(previousApprox2, gcvops.Mult(fPA2, ratioA))

		fCA, errfCA = evalV(f, currentApprox)
		if errfCA != nil {
			return nil, errfCA
		}
	}

	if solutionFound {
		return root, nil
	}

	return nil, errors.New("Unable to find root of given function")
}

// Steffensen1D is for solving the 1D root finding Steffensen's mehtod
func Steffensen1D(initialApprox float64, TOL float64, maxIteration int, f *gcf.Function) (gcv.Value, error) {
	previousApprox1 := initialApprox
	previousApprox2 := f.MustEval(previousApprox1).Value().Real()
	previousApprox3 := f.MustEval(previousApprox2).Value().Real()
	currentApprox := previousApprox1 - math.Pow((previousApprox2-previousApprox1), 2)/(previousApprox3-2*previousApprox2+previousApprox1)
	root := float64(0)
	solutionFound := false

	for i := 0; i < maxIteration; i++ {
		if math.Abs(currentApprox-previousApprox1) < TOL {
			root = currentApprox
			solutionFound = true
			break
		}

		previousApprox1 = currentApprox
		previousApprox2 = f.MustEval(previousApprox1).Value().Real()
		previousApprox3 = f.MustEval(previousApprox2).Value().Real()

		currentApprox = previousApprox1 - math.Pow((previousApprox2-previousApprox1), 2)/(previousApprox3-2*previousApprox2+previousApprox1)
	}

	if solutionFound {
		return gcv.MakeValue(root), nil
	}

	return nil, errors.New("Unable to find root of given function")
}
