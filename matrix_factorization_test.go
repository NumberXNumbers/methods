package methods

import (
	"reflect"
	"testing"

	m "github.com/NumberXNumbers/types/gc/matrices"
	mops "github.com/NumberXNumbers/types/gc/matrices/ops"
	v "github.com/NumberXNumbers/types/gc/vectors"
)

func TestLU(t *testing.T) {
	testVectorAa := v.MakeVector(v.RowSpace, 0, 0, -1, 1)
	testVectorAb := v.MakeVector(v.RowSpace, 1, 1, -1, 2)
	testVectorAc := v.MakeVector(v.RowSpace, -1, -1, 2, 0)
	testVectorAd := v.MakeVector(v.RowSpace, 1, 2, 0, 2)
	testVectorsA := v.MakeVectors(v.RowSpace, testVectorAa, testVectorAb, testVectorAc, testVectorAd)
	testMatrixA := m.MakeMatrixAlt(testVectorsA)

	L, U, P, errA := LU(testMatrixA)

	if errA != nil {
		t.Error("Unexpected Error")
	}

	P.Trans()
	A := mops.MustMultSimple(mops.MustMultSimple(P, L), U)
	if !reflect.DeepEqual(testMatrixA, A) {
		t.Errorf("Expected %+v, received %+v", testMatrixA, A)
	}

	testVectorsB := v.MakeVectors(v.RowSpace, testVectorAa, testVectorAb, testVectorAc)
	testMatrixB := m.MakeMatrixAlt(testVectorsB)

	_, _, _, errB := LU(testMatrixB)
	if errB == nil {
		t.Error("Expected Error")
	}

	testVectorCa := v.MakeVector(v.RowSpace, 0, 0, -1)
	testVectorCb := v.MakeVector(v.RowSpace, 0, 1, -1)
	testVectorCc := v.MakeVector(v.RowSpace, 0, -1, 2)
	testVectorsC := v.MakeVectors(v.RowSpace, testVectorCa, testVectorCb, testVectorCc)
	testMatrixC := m.MakeMatrixAlt(testVectorsC)

	_, _, _, errC := LU(testMatrixC)
	if errC == nil {
		t.Error("Expected Error")
	}
}

func TestLDLtrans(t *testing.T) {
	testVectorAa := v.MakeVector(v.RowSpace, 4, -1, 1)
	testVectorAb := v.MakeVector(v.RowSpace, -1, 4.25, 2.75)
	testVectorAc := v.MakeVector(v.RowSpace, 1, 2.75, 3.5)
	testVectorsA := v.MakeVectors(v.RowSpace, testVectorAa, testVectorAb, testVectorAc)
	testMatrixA := m.MakeMatrixAlt(testVectorsA)

	L, D, errA := LDLt(testMatrixA)
	if errA != nil {
		t.Error("Unexpected Error")
	}

	A := mops.MustMultSimple(mops.MustMultSimple(L, D), m.MakeTransMatrix(L))
	if !reflect.DeepEqual(testMatrixA, A) {
		t.Errorf("Expected %+v, received %+v", testMatrixA, A)
	}

	testVectorsB := v.MakeVectors(v.RowSpace, testVectorAa, testVectorAb)
	testMatrixB := m.MakeMatrixAlt(testVectorsB)

	_, _, errB := LDLt(testMatrixB)
	if errB == nil {
		t.Error("Expected Error")
	}

	testVectorCa := v.MakeVector(v.RowSpace, 0, -1, 0)
	testVectorCb := v.MakeVector(v.RowSpace, -1, 0, 2.75)
	testVectorCc := v.MakeVector(v.RowSpace, 0, 2.75, 0)
	testVectorsC := v.MakeVectors(v.RowSpace, testVectorCa, testVectorCb, testVectorCc)
	testMatrixC := m.MakeMatrixAlt(testVectorsC)

	_, _, errC := LDLt(testMatrixC)
	if errC == nil {
		t.Error("Expected Error")
	}
}
