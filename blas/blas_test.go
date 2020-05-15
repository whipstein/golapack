package blas

import (
	"strconv"
	"testing"
)

func TestLsame(t *testing.T) {
	testvals :=[]struct {
		x        byte
		y        byte
		exemplar bool
	}{
		{
			'A',
			'A',
			true,
		},
		{
			'A',
			'a',
			true,
		},
		{
			'A',
			'B',
			false,
		},
		{
			'A',
			'b',
			false,
		},
		{
			'B',
			'A',
			false,
		},
		{
			'b',
			'A',
			false,
		},
	}
	for _, val := range testvals {
		if out := Lsame(&val.x, &val.y); out != val.exemplar {
			t.Errorf("Unexpected result for The Answer. Input: %c and %c  Got: %s  want: %s", val.x, val.y, strconv.FormatBool(out), strconv.FormatBool(val.exemplar))
		}
	}
}
