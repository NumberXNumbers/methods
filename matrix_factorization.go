package methods

import (
	"errors"

	m "github.com/NumberXNumbers/types/gc/matrices"
	gcv "github.com/NumberXNumbers/types/gc/values"
	gcvops "github.com/NumberXNumbers/types/gc/values/ops"
)

// LU will return an L U decomposition of matrix A as well as any permutation matrix P, else error
func LU(A m.Matrix) (L, U, P m.Matrix, err error) {
	if !A.IsSquare() {
		return nil, nil, nil, errors.New("Matrix is not square")
	}
	degree, _ := A.Dim()
	matrixCopy := A.Copy()
	p := m.NewIdentityMatrix(degree)
	l := m.NewMatrix(degree, degree)
	u := m.NewMatrix(degree, degree)
	for i := 0; i < degree; i++ {
		sumUII := gcv.NewValue()
		for k := 0; k < i; k++ {
			sumUII = gcvops.Add(sumUII, gcvops.Mult(l.Get(i, k), u.Get(k, i)))
		}
		valueU := gcvops.Sub(matrixCopy.Get(i, i), sumUII)
		if i != degree-1 && valueU.Complex() == 0 {
			count := i + 1
			for count < degree {
				sumUII := gcv.NewValue()
				for k := 0; k < i; k++ {
					sumUII = gcvops.Add(sumUII, gcvops.Mult(l.Get(count, k), u.Get(k, i)))
				}
				valueAtCount := gcvops.Sub(matrixCopy.Get(count, i), sumUII)
				if valueAtCount.Complex() != 0 {
					valueU = valueAtCount

					matrixCopy.Swap(i, count)
					p.Swap(i, count)
					l.Swap(i, count)

					break
				}
				count++
				if count == degree {
					return nil, nil, nil, errors.New("Unable to factorize matrix")
				}
			}
		}

		u.Set(i, i, valueU)
		l.Set(i, i, 1)

		for j := i + 1; j < degree; j++ {
			sumU := gcv.NewValue()
			sumL := gcv.NewValue()
			for k := 0; k < i; k++ {
				sumU = gcvops.Add(sumU, gcvops.Mult(l.Get(i, k), u.Get(k, j)))
				sumL = gcvops.Add(sumL, gcvops.Mult(l.Get(j, k), u.Get(k, i)))
			}
			u.Set(i, j, gcvops.Sub(matrixCopy.Get(i, j), sumU))
			l.Set(j, i, gcvops.Div(gcvops.Sub(matrixCopy.Get(j, i), sumL), u.Get(i, i)))
		}
	}
	return l, u, p, nil
}
