package methods

import (
	"errors"

	m "github.com/NumberXNumbers/types/gc/matrices"
	v "github.com/NumberXNumbers/types/gc/vectors"
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
		sumUII := gcv.Zero()
		for k := 0; k < i; k++ {
			sumUII = gcvops.Add(sumUII, gcvops.Mult(l.Get(i, k), u.Get(k, i)))
		}
		valueU := gcvops.Sub(matrixCopy.Get(i, i), sumUII)
		if i != degree-1 && valueU.Complex() == 0 {
			count := i + 1
			for count < degree {
				sumUII := gcv.Zero()
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
			sumU := gcv.Zero()
			sumL := gcv.Zero()
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

// LDLt will return the L D mattrices of the L D Lt facorization of matrix A if it is
// hermitian positive-definite else error
// TODO: add in checks for positive-definite and hermitian forms
func LDLt(A m.Matrix) (L, D  m.Matrix, err error) {
	if !A.IsSquare() {
		return nil, nil, errors.New("Matrix is not square")
	}
	degree, _ := A.Dim()
	l := m.NewMatrix(degree, degree)
	d := m.NewMatrix(degree, degree)

	for i := 0; i < degree; i++ {
		sumLi := gcv.Zero()
		diff := gcv.Zero()
		sumLj := gcv.Zero()
		counter := 0
		tempVector := v.NewVector(v.RowSpace, i+1)
		for j := 0; j < degree; j++ {
			if j <= i - 1 {
				tempValue := gcvops.Mult(l.Get(i, j), d.Get(j, j))
				tempVector.Set(j, tempValue)
				sumLi = gcvops.Add(sumLi, gcvops.Mult(l.Get(i, j), tempValue))
			} else if j == i {
				diff = gcvops.Sub(A.Get(i, i), sumLi)
				if diff.Complex() == 0 {
					return nil, nil, errors.New("Matrix can not be factorized")
				}
				d.Set(i, i, diff)
				l.Set(i, i, 1)
			} else {
				if j > tempVector.Len() {
					sumLi = gcv.Zero()
				} else {
					sumLj = gcvops.Add(sumLj, gcvops.Mult(l.Get(j, counter), tempVector.Get(counter)))
				}
				div := gcvops.Div(gcvops.Sub(A.Get(j, i), sumLj), diff)
				l.Set(j, i, div)
				counter++
			}
		}
	}
	return l, d, nil
}
